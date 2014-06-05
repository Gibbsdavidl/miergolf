

''' 
Run this like:
ipython flowSim.py config.txt graphfile.txt 100 100 6 
where we specify (in the config file):

   config file name, 
   the graph file name,
   the number of nodes in the graph, 
   the number of edges, 
   the power law degree,
   timesteps
   dissipation rate (in decimal)
   fullrun

1.) The graph is randomly generated, with random edge weights, and written as a file
2.) The graph is transformed to a line graph, and that's written to a file
3.) A simulation is run where:
        * Info generation & flow decisions *
        A. First each node generates an info-block, tagged with source and time step
        B. Then each node randomly chooses to dissipate a block, or transfer
                If transferring, this is done probabilistically based on normalized edge weights
        * syncronous update of blocks *
        C. Then each block moves blocks, either dissipating or moving to another node
                For each node, we keep track of
                    where the info originated from
                    the path the info took
                    ...

                    
'''



def main():
    # parse command line options
    try:
        opts, args = getopt.getopt(sys.argv[1:], "h", ["help"])
    except getopt.error, msg:
        print msg
        print "for help use --help"
        sys.exit(2)
    # process options
    for o, a in opts:
        if o in ("-h", "--help"):
            print __doc__
            sys.exit(0)

    # process arguments
    state = initProgramState(args)
    state = initSimState(args,state)

    # generate a random graph .... 
    nodenames = randomGraph(state)

    # convert to the linegraph and print it out for viz
    (nodes, sparseMat) = buildGraph (state)
    printLineGraph(state, nodes, sparseMat)
    
    # run the explicit flow
    (rxhistory,nstore) = flowtron2(state, sparseMat, nodes)
     
    # some statistics on the flow
    counts = processSim(rxhistory, nstore)
    txhistory = counts[6]
    rxhistory = counts[7]
    
    # output results ... did we find the solution? OR how much of it?
    # we should plot out graphs, and such.
    printResults(state, nodes, counts)

    # then search for the solution.
    (score,soln,rxed,txed) = search(nodes, state, txhistory, rxhistory)

    printRXTXTables(state, state["timesteps"], soln, rxhistory, txhistory)

    wts = weightsum(nodes,soln)
    
    if state["full"] == "fullrun":
        # process arguments
        sys.stderr.write("\nRunning optimization..\n")
        cpus = state["cpus"]
        pool = mp.Pool(cpus)
        s2 = optimize (pool, state, nodes, sparseMat)
        pool.close()

        print ("Config\tGraph\tType\tOptimTo\tMode\tSteps\tNodes\tAnts\tTx\tRx\tDamp\tLocal\tSimScore\tSimSoln\tAntScore\tAntSoln\tWt")
        print ( str(state["config"]) +"\t"+
                str(state["graphfile"]) +"\t"+
                str(state["lineGraph"]) +"\t"+
                str(state["opton"]) +"\t"+
                str(state["mode"]) +"\t"+
                str(state["timesteps"]) + "\t"+
                str(len(nodes)) + "\t"+
                str(state["ants"]) + "\t"+
                str(state["tx"]) +"\t"+
                str(state["rx"]) +"\t"+
                str(state["damp"]) +"\t"+
                str(state["local"]) +"\t"+
                str(score) +"\t"+
                str(soln) +"\t"+
                str(state["bestEver"][1]) +"\t"+
                str(state["bestEver"][2]) +"\t"+
                str(wts))
    else:
        print "Done:"
        print score
        print soln
        print rxed
        print txed



def printLineGraph(state, nodes, sparseMat):
    filename = state["graphfile"]+".linegraph.txt"
    fout = open(filename,'w')
    n = len(nodes)
    for i in xrange(n):
        for j in xrange(n):
            p = sparseMat.getrow(i).getcol(j).toarray().flatten()[0]
            if p > 0.0:
                fout.write(str(i) +"\t" + str(j) +"\t" + str(p) +"\n") 
    fout.close()
    sys.stderr.write("Printed line graph.\n")


def randomGraph(state):
    # this function is going to generate a graph, and write it
    # to a file, as specified in the config file.
    # ns: number of nodes, es: number of edges, deg: power law degree, s: state 
    ns = state["nodes"]
    es = state["edges"]
    deg = state["degpow"]
    #g1 = Graph.Static_Power_Law(ns,es,deg)
    #gedgesOrdered = g1.get_edgelist()
    #gedgesOrdered = HairBall(ns)
    #gedgesOrdered = BAModel(1, ns, 0.1, 0.2)
    gedgesOrdered = Bollobas(0.4, 0.2, 0.4, 0.1, 0.0, 0.2, ns)
    gedges = []
    for ge in gedgesOrdered:  # give the edges a random direction
        if np.random.random() > 0.5:
            gedges.append( (ge[1],ge[0]) )
        else:
            gedges.append(ge)
    nodenames = set()
    for (e1,e2) in gedges: nodenames.add(e1); nodenames.add(e2);
    nodenames = list(nodenames)
    wts = [np.random.random() for i in xrange(len(gedges))]
    fout = open(state["graphfile"],'w')
    for i in xrange(len(gedges)):
        fout.write(str(gedges[i][0]) +"\t"+ str(gedges[i][1]) +"\t"+ str(wts[i]) +"\n")
    print("Done writing random graph.")
    return(nodenames)


def partitionGraph(nodes, state):
    n = len(nodes)
    s = np.random.choice(range(0,n), int(state["k"]))
    t = [i for i in range(0,n) if i not in s]
    return(s,t)


def newInfoBlock(i,step):
    # the source and time should uniquely identify the block #
    d = dict()
    d["source"] = i   # the node that created the info-block
    d["path"] = [i]   # the path of nodes it followed 
    d["time"] = step  # when it was created
    return(d)


def initNodeList(n):
    d = dict()
    for ni in xrange(n):
        d[ni] = []
    return(d)

    
def initNodeSet(n):
    d = dict()
    for ni in xrange(n):
        d[ni] = set()
    return(d)


def flowstatus(liner,step):
    if step % 97 == 0:
        sys.stderr.write(".")
        if liner < 0:
            sys.stderr.write("\n")
            liner = 80
        liner -= 1
    return(liner)


def nextState(cps, r):
    idx = 0
    for c in cps:
        if r < c:
            return idx
        else:
            idx += 1
    return (len(cps)-1)


def flowtron(state, sparseMat, nodes):
    nids = range(0,len(nodes))
    disp = state["disp"]
    steps = state["timesteps"]
    n = len(nodes)                # number of nodes in the line graph
    nodestore   = initNodeList(n) # holds the infoblocks
    nodehistory = initNodeList(n) # holds the list of where info came from
    sys.stderr.write("Running simulation for " +str(steps)+ " timesteps.\n")
    liner = 80
    # first produce new information and determine transitions #
    for step in xrange(steps):
        liner = flowstatus(liner,step)
        movement = [0 for i in xrange(n)]                 # record where each node will transfer to..
        whichblock = [0 for i in xrange(n)]
        for ni in xrange(n):                              # for each node
            nodestore[ni].append(newInfoBlock(ni, step))      # generate new information
            ps = sparseMat.getrow(ni).toarray().flatten()     # the transition probs
            ri = np.random.random()    
            if all(ps == 0.0) or (ri < disp):                 # dissipate?
                movement[ni] = -1
            else:
                cps = np.cumsum(ps)
                rnd = np.random.random() * cps[-1]
                nj = np.searchsorted(cps, rnd)                # choose an edge
                movement[ni] = nj
                whichblock[ni] = np.random.randint(0,len(nodestore[ni]),1)  # choose a block

        # then do a syncronous info-block move or dissipation. #
        for ni in xrange(n):
            nj = movement[ni]                          # the move destination
            ib = nodestore[whichblock[ni]]
            if nj == -1:                               # dissipate a random block, put it in the dustbin
                nodestore[ni].remove(ib)
                #dustbin.append(ib)               # throw it in the dustbin  
            else:
                ib["path"].append(nj)            # record where this block has been
                nodehistory[nj].append(ib["source"])  # node nj has received a block with source ib["source"]
                nodestore[ni].remove(ib)
                nodestore[nj].append(ib)
    sys.stderr.write("\n")
    return((nodehistory, nodestore, dustbin)) # nodehistory is rxhist


def updateNodeStore(nodestore, node, block, storesize):
    if len(nodestore[node]) > storesize:
        i = np.random.randint(0,storesize)
        nodestore[node].remove(nodestore[node][i])
    nodestore[node].append(block)
    return(nodestore)

# there's an issue where some nodes, have a massive store of
# info-blocks, dominated by a small set of nodes, and the
# minority member blocks, are passed up for selection for passage.
# A solution might be to hold a set of blocks rather than a list.
# or to support a limited store size.
def flowtron2(state, sparseMat, nodes):
    nids = range(0,len(nodes))
    disp = state["disp"]
    steps = state["timesteps"]
    n = len(nodes)                # number of nodes in the line graph
    nodestore   = initNodeList(n) # holds the infoblocks
    nodehistory = initNodeSet(n) # holds the list of where info came from
    # nodestore = np.zeros( (n,n) )   # how much info contained in the store? Just an int at the right location
    # nodehistory = np.zeros( (n,n) ) # make updating the store and array an O(1) operation
    # movement = np.zeros(n)          # where is each going? make sure each gets written to at each iter
    sys.stderr.write("Running simulation for " +str(steps)+ " timesteps.\n")
    liner = 80
    # first produce new information and determine transitions #
    for step in xrange(steps):
        liner = flowstatus(liner,step)
        movement = [0 for i in xrange(n)]                 # record where each node will transfer to..
        # elim the above line
        for ni in xrange(n):                              # for each node
            nodestore = updateNodeStore(nodestore, ni, newInfoBlock(ni,step), state["storesize"])
            # nodestore[ni]
            ps = sparseMat.getrow(ni).toarray().flatten()     # the transition probs
            ri = np.random.random()    
            if all(ps == 0.0) or (ri < disp):                 # dissipate
                movement[ni] = -1
            else:
                cps = np.cumsum(ps)
                rnd = np.random.random() * cps[-1]
                nj = np.searchsorted(cps, rnd)
                movement[ni] = nj
        # then do a syncronous info-block move or dissipation. #
        for ni in xrange(n):
            nj = movement[ni]                          # the move destination
            randblocki = np.random.randint(0, len(nodestore[ni]))
            ib = nodestore[ni][randblocki]  # we just generated one, so there has to be at least 1!
            if nj == -1:                               # dissipate a random block, put it in the dustbin
                nodestore[ni].remove(ib)
            else:
                nodehistory[nj].add(ib["source"])  # node nj has received a block with source ib["source"]
                nodestore[ni].remove(ib)
                nodestore = updateNodeStore(nodestore, nj, ib, state["storesize"])
    sys.stderr.write("\n")
    return((nodehistory, nodestore)) # nodehistory is rxhist

    
    
def processSim(nhist, nstore):
    sys.stderr.write("Processing Simulation...\n")
    rxhist = dict()
    n = len(nhist)
    storG = np.zeros(n)  # the amt of info generated by self
    storR = np.zeros(n)  # the amt of into received by others
    totalR = np.zeros(n) # total amt received
    uniqR = np.zeros(n)  # the amt of unique received info
    uniqT = np.zeros(n)  # the unique number of nodes sent to
    totalT = np.zeros(n) # total amt sent out
    txhist = []          # where the info *ultimately* ended up
    for ni in xrange(n):
        rxhist[ni] = list(nhist[ni])
        storG[ni] = sum(map(lambda y: y == ni, map(lambda x: x["source"],nstore[ni])))
        storR[ni] = sum(map(lambda y: y != ni, map(lambda x: x["source"],nstore[ni])))
        uniqR[ni] = len(set(nhist[ni]))
        totalR[ni] = len(nhist[ni])
        x = 0; y = 0; l = [];
        for nj in xrange(n):
            if ni != nj:
                x += int(ni in nhist[nj])
                yi = sum(map(lambda z: z == ni, nhist[nj]))
                l = l + list(np.repeat(nj,yi)) # ni transmitted to nj yi times.
                y += yi
        uniqT[ni] = x
        totalT[ni] = y
        txhist.append(l)
    return( (storG, storR, totalR, uniqR, uniqT, totalT, txhist, rxhist) )
            

def nuniq(x):
    return(len(set(x)))
        

def search(nodes, state, txhist, rxhist):
    sys.stderr.write("Searching Solutions...\n")
    cpus = int(state["cpus"])
    pool = mp.Pool(cpus)
    n = len(nodes)
    idx = range(0,n)
    bestScore = 0; bestSoln = []; bestList = [];
    rxed = []; txed = []
    step = 1
    liner = 80
    alltups = []

    # get all the possible solutions together
    for tup in it.combinations(idx, int(state["k"])):  # for each combination of k nodes
        alltups.append(tup)
       
    # split up the alltups into cpus number of lists
    tuplist = np.array_split(np.array(alltups), cpus)

    # build up a list of tuples with rxhist, txhist, and tups
    dat = it.izip( it.repeat(state, cpus),
                   it.repeat(nodes, cpus),
                   it.repeat(rxhist,cpus), 
                   it.repeat(txhist,cpus), tuplist)
    allscores = pool.map(subSearch, dat)

    maxi = 0; bestScore = 0.0
    for score in allscores:
        if score[0] > bestScore:
            bestScore = score[0]
            bestSoln = score[1]
            rxed = score[2]
            txed = score[3]

    pool.close()
    return( (bestScore, bestSoln, rxed, txed) )


def subSearch( (state, nodes, rxhist, txhist, tups) ):
    # score each tup
    sys.stderr.write("Working on " + str(len(tups)) + " number of solutions...\n")
    bestScore = 0.0; bestSoln = []; bestList = []
    for tup in tups:
        subscore = []; tx = []; rx = [];               # gather all nodes that we TXed to, or RXed from 
        for t in tup:
            if state["mode"] == "tx":
                subscore += txhist[t]
                tx += txhist[t]
            elif state["mode"] == "rx":
                subscore += rxhist[t]
                rx += rxhist[t]
            else:
                subscore += txhist[t] + rxhist[t]
                tx += txhist[t]
                rx += rxhist[t]
        subscore = [si for si in subscore if si not in tup] # don't score what's in tup!
        score = nuniq(subscore) # number of unique members .. PLUS ..
        score = score + weightsum(nodes,tup)
        if score > bestScore:
            bestScore = score
            bestSoln = tup
            bestList = subscore
            rxed = rx; txed = tx
    return( (bestScore, bestSoln, rxed, txed) )


def weightsum(nodes,tup):
    totwt = 0.0
    for ti in tup:
        totwt += nodes[ti][2]
    return(totwt)


def printResults(state, nodes, (storG, storR, totalR, uniqR, uniqT, totalT, txhist, rxhist)):
    fout = open((state["graphfile"]+"simResultsTable.txt"),'w')
    n = len(nodes)
    fout.write("From\tTo\tWt\tStoreGen\tStoreRx\ttotalRX\tuniqRX\ttotalTX\tuniqTX\n")
    for ni in xrange(n):
        line = (nodes[ni][0] +"\t"+ nodes[ni][1] +"\t"+
                str(round(nodes[ni][2],3)) +"\t"+
                str(storG[ni]) +"\t"+ str(storR[ni]) +"\t"+
                str(totalR[ni]) +"\t"+ str(uniqR[ni]) +"\t"+
                str(totalT[ni]) +"\t"+ str(uniqT[ni]) +"\n")
        fout.write(line)
    fout.close()


def printRXTXTables(state, steps, soln, rx, tx):
    fout = open((state["graphfile"]+"simFtables.txt"),'w')
    n = len(rx)
    soln = list(soln)
    ts = [i for i in xrange(n) if i not in soln]
    for ti in ts:
        fout.write(str(ti)+"\t")
    fout.write("\n")
    for i in soln:
        fout.write(str(i)+"\t")
        rs = np.array(rx[i])
        for j in ts:
            if j in rs:
                numjs = float(sum(rs == j))/float(steps)
                fout.write(str(numjs) + "\t")
            else:
                fout.write("0.0\t")
        fout.write("\n")
    fout.close()
    fout = open((state["graphfile"]+"simHtables.txt"),'w')
    for ti in ts:
        fout.write(str(ti)+"\t") # column name
    fout.write("\n")
    for i in soln:
        fout.write(str(i)+"\t")  # row name
        rs = np.array(tx[i])
        for j in ts:
            if j in rs:
                numjs = float(sum(rs == j))/float(steps)
                fout.write(str(numjs) + "\t")
            else:
                fout.write("0.0\t")
        fout.write("\n")
    fout.close()


def printColElem(sm, ts, i):
    xs = sm.getcol(i).toarray().flatten()
    for k in xrange(len(ts)):
    	 print str(k) +"\t"+ str(ts[k]) +"\t"+ str(xs[k])


def printRowElem(sm, ts, i):
    xs = sm.getrow(i).toarray().flatten()
    for k in xrange(len(ts)):
    	 print str(k) +"\t"+ str(ts[k]) +"\t"+ str(xs[k])


def writeMatrix(mat, fileout):
    fout = open(fileout,'w')
    (m,n) = mat.get_shape()  # m - rows, n - cols
    for i in xrange(m):
        for j in xrange(n):
            fout.write(str(mat[i,j]) + "\t")
        fout.write("\n")
    fout.close()


if __name__ == "__main__":
    import sys
    import getopt
    from math import sqrt, ceil
    import multiprocessing as mp
    import numpy as np
    import itertools as it
    from programState import *
    from graph import *
    from antopt import *
    from randomNets import *
    main()


# args = ["../../Sim/config_test.txt", "../../Sim/simgraph.txt", "10","10","4","10", "0.1"]
