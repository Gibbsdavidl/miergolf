

''' 
Run this like:
ipython flowSim.py config.txt graphfile.txt 100 100 6 
where we specify the:

   config file name, 
   the graph file name,
   the number of nodes in the graph, 
   the number of edges, 
   the power law degree,
   timesteps
   dissipation rate (in decimal)
   --fullrun  or --notfullrun

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
    state = initProgramState (args)

    # generate a random graph .... 
    nodenames = randomGraph(int(args[2]), int(args[3]), int(args[4]), state)

    # convert to the linegraph and print it out for viz
    (nodes, sparseMat) = buildGraph (state)
    printLineGraph(state, nodes, sparseMat)
    
    # run the explicit flow
    timesteps = int(args[5])
    dissipation = float(args[6])
    (txhist,nstore,nbin) = flowtron(state, sparseMat, nodes, nodenames, dissipation, timesteps)

    # some statistics on the flow
    counts = processSim(txhist, nstore)
    rxhist = counts[6]
    
    # output results ... did we find the solution? OR how much of it?
    # we should plot out graphs, and such.
    printResults(nodes, counts)

    # then search for the solution.
    (score,soln,rxed,txed) = search(nodes, state, counts, txhist)

    if args[7] == "fullrun":
        # process arguments
        sys.stderr.write("\nRunning optimization..\n")
        s = initProgramState (args)
        cpus = s["cpus"]
        pool = mp.Pool(cpus)
        (nodes, sparseMat) = buildGraph (s)
        s2 = optimize (pool, s, nodes, sparseMat)
        pool.close()

        print ("Config\tGraph\tType\tOptimTo\tMode\tAnts\tTx\tRx\tDamp\tLocal\tSimScore\tSimSoln\tAntScore\tAntSoln")
        print ( str(s["config"]) +"\t"+
                str(s["graphfile"]) +"\t"+
                str(s["lineGraph"]) +"\t"+
                str(s["opton"]) +"\t"+
                str(s["mode"]) +"\t"+
                str(s["ants"]) + "\t"+
                str(s["tx"]) +"\t"+
                str(s["rx"]) +"\t"+
                str(s["damp"]) +"\t"+
                str(s["local"]) +"\t"+
                str(score) +"\t"+
                str(soln) +"\t"+
                str(s["bestEver"][1]) +"\t"+
                str(s["bestEver"][2]) + "\t")

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


    
def randomGraph(ns, es, deg, state):
    # this function is going to generate a graph, and write it
    # to a file, as specified in the config file.
    # ns: number of nodes, es: number of edges, deg: power law degree, s: state 
    g1 = Graph.Static_Power_Law(ns,es,deg)
    gedgesOrdered = g1.get_edgelist()
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


def flowtron(state, sparseMat, nodes, nodenames, disp, steps):    
    n = len(nodes)                # number of nodes in the line graph
    nodestore   = initNodeList(n) # holds the infoblocks
    nodehistory = initNodeList(n) # holds the list of where info came from
    dustbin = []                  # the dust bin of history, where info-blocks go to retire
    sys.stderr.write("Running simulation for " +str(steps)+ " timesteps.\n")

    # first produce new information and determine transitions #
    for step in xrange(steps):
        movement = [0 for i in xrange(n)]                 # record where each node will transfer to..
        for ni in xrange(n):                              # for each node
            nodestore[ni].append(newInfoBlock(ni, step))      # generate new information
            ps = sparseMat.getrow(ni).toarray().flatten()     # the transition probs
            ri = np.random.random()    
            if all(ps == 0) or (ri < disp):                   # dissipate
                movement[ni] = -1
            else:                                             # or transfer the block
                psi = np.random.random()
                tps = np.where(ps > 0)                           # possible transitions, small vec
                cps = ps.cumsum()[ps > 0]                        # the cumulative probabilities
                nj = tps[0][(cps >= psi).argmax()]               # where we're going
                movement[ni] = nj

        # then do a syncronous info-block move or dissipation. #
        for ni in xrange(n):
            nj = movement[ni]                # the move destination
            ib = np.random.choice(nodestore[ni],1)[0]  # we just generated one, so there has to be at least 1!
            if nj == -1:                     # dissipate a random block, put it in the dustbin
                nodestore[ni].remove(ib)
                dustbin.append(ib)               # throw it in the dustbin  
            else:
                ib["path"].append(nj)            # record where this block has been
                nodehistory[nj].append(ib["source"])  # node nj has received a block with source ib["source"]
                nodestore[ni].remove(ib)
                nodestore[nj].append(ib)

    return((nodehistory, nodestore, dustbin))


    
def processSim(nhist, nstore):
    sys.stderr.write("Processing Simulation...\n")
    n = len(nhist)
    storG = np.zeros(n)  # the amt of info generated by self
    storR = np.zeros(n)  # the amt of into received by others
    totalR = np.zeros(n) # total amt received
    uniqR = np.zeros(n)  # the amt of unique received info
    uniqT = np.zeros(n)  # the unique number of nodes sent to
    totalT = np.zeros(n) # total amt sent out
    rxhist = []          # where the info *ultimately* ended up
    for ni in xrange(n):
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
        rxhist.append(l)
    return( (storG, storR, totalR, uniqR, uniqT, totalT, rxhist) )
            

def printResults(nodes, (storG, storR, totalR, uniqR, uniqT, totalT, rxhist)):
    fout = open("simResultsTable.txt",'w')
    n = len(nodes)
    fout.write("From\tTo\tWt\tStoreGen\tStoreRx\ttotalRX\tuniqRX\ttotalTX\tuniqTX\n")
    for ni in xrange(n):
        line = (nodes[ni][0] +"\t"+ nodes[ni][1] +"\t"+ str(round(nodes[ni][2],3)) +"\t"+
                str(storG[ni]) +"\t"+ str(storR[ni]) +"\t"+
                str(totalR[ni]) +"\t"+ str(uniqR[ni]) +"\t"+
                str(totalT[ni]) +"\t"+ str(uniqT[ni]) +"\n")
        fout.write(line)


def nuniq(x):
    return(len(set(x)))
        

def search(nodes, state, counts, txhist):
    rxhist = counts[6]
    n = len(nodes)
    idx = range(0,n)
    bestScore = 0; bestSoln = []; bestList = [];
    rxed = []; txed = []
    for tup in it.combinations(idx, int(state["k"])):  # for each combination of k nodes
        subscore = []; tx = []; rx = [];               #    gather all nodes that we TXed to, or RXed from 
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
        score = nuniq(subscore)
        if score > bestScore:
            bestScore = score
            bestSoln = tup
            bestList = subscore
            rxed = rx; txed = tx
    return( (bestScore, bestSoln, rxed, txed) )

        
if __name__ == "__main__":
    import sys
    import getopt
    import multiprocessing as mp
    from igraph import Graph
    import numpy as np
    import itertools as it
    from programState import *
    from graph import *
    from antopt import *
    main()


# args = ["../../Sim/config_test.txt", "../../Sim/simgraph.txt", "10","10","4","10", "0.1"]
