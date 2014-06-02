
from math import ceil, sqrt, pow
import numpy as np
from igraph import *


def randSelect(ps):
    cps = np.cumsum(ps)
    rnd = np.random.random() * cps[-1]
    return(np.searchsorted(cps, rnd))

    
def uncorScaleFreeNetwork(nodes,degpow,mindeg):
    network = []
    degmax = int(ceil(sqrt(nodes)))
    probDist = [pow(x, (-degpow)) for x in range(mindeg,degmax)]
    nodeDegs = [min(randSelect(probDist)+mindeg,degmax) for xi in xrange(nodes)]
    totalEdges = sum(nodeDegs)/2
    if totalEdges%2 != 0:
        idx = np.random.choice(range(0,nodes),1)
        nodeDegs[idx] +=1
        totalEdges = sum(nodeDegs)/2
    for e in xrange(totalEdges):
        i = randSelect(nodeDegs) # pick a node
        nodeDegs[i] -= 1         # decrement the degree on node i
        j = randSelect(nodeDegs)
        while i == j:
            j = randSelect(nodeDegs)
        nodeDegs[j] -= 1
        network.append( (i,j) )
    for i in xrange(len(network)):
        network[i] = ( network[i][0]+1 , network[i][1]+1)
    return(network)


def Bollobas(alpha, beta, gamma, deltaIn, deltaOut, randDir, timesteps):
    #http://research.microsoft.com/en-us/um/people/borgs/Papers/dirSCgrph.pdf
    urand = np.random.random

    # alpha+beta+gamma == 1, alpha+gamma > 0, deltaOut == 0, deltaIn > 0
    if alpha+beta+gamma < 1:
        print "alpha+beta+gamma < 1, alpha+gamma > 0, deltaOut == 0, deltaIn > 0"
        sys.exit(1)
    
    # start with a two nodes, network is G
    nextNode = 2;
    net    = [(1,0)]
    degIn  = [1,0]
    degOut = [0,1]
    
    for t in range(1,timesteps):
        if urand() < alpha:
            # add vertex, connect to vj given by deltaIn+degIn
            pdenom = t + deltaIn*nextNode
            ps = np.array([(d+deltaIn)/pdenom for d in degIn])
            nj = randSelect(ps[0:nextNode])
            degOut.append(0); degIn.append(0);
            degIn[nj] += 1
            degOut[nextNode] += 1
            net.append( (nextNode, nj) )
            nextNode +=1

        if urand() < beta:
            # connect two existing nodes, first the in coming
            pdenom = t + deltaIn*nextNode
            ps = np.array([(d+deltaIn)/pdenom for d in degIn])
            ni = randSelect(ps[0:nextNode])
            degIn[ni] += 1
            # out going
            pdenom = t + deltaOut*nextNode
            ps = np.array([(d+deltaOut)/pdenom for d in degOut])
            nj = randSelect(ps[0:nextNode])
            degOut[nj] += 1
            net.append( (nj, ni) )
                        
        if urand() < gamma:
            # add vertex, add edge from nj to ni, choose by degOut+deltaOut
            pdenom = t + deltaOut*nextNode
            ps = np.array([(d+deltaOut)/pdenom for d in degOut])
            nj = randSelect(ps[0:nextNode])
            degOut.append(0); degIn.append(0);
            degIn[nextNode] += 1
            degOut[nj] += 1
            net.append( (nj, nextNode) )
            nextNode +=1

    net = list(set(net))  # eliminate multiple edges
    for i in xrange(len(net)):
        if np.random.random() < randDir: # randomize some proportion of edges
            net[i] = (net[i][0]+1, net[i][1]+1) # if the edge dirs are not
        else:                                   # randomized, the ant solutions
            net[i] = (net[i][1]+1, net[i][0]+1) # are singular... and fail.
    return(net)


# writeNetwork(Bollobas(0.4, 0.2, 0.4, 0.1, 0.0, 0.2, 100), "net6.txt", "F")
#> l <- as.matrix(read.table("net6.txt")); g <- graph.edgelist(l[,1:2]+1)
#> plot(g, vertex.size=0.5, vertex.label.cex=0.5, edge.arrow.size=0.75)

    
def BAModel_Ver2(m0, nNodes, zeroAppeal):
    sample=np.random.choice
    deg = np.zeros(nNodes)
    net   = [(ni, sample(nodes)) for ni in nodes]
    for ni in net:
        deg[ni[0]] += 1
        deg[ni[1]] += 1
    for ni in range(m0,nNodes):
        for ei in xrange(m0):
            degsum = deg.sum()
            ps = np.array([d/degsum+zeroAppeal for d in deg])
            ps = ps/ps.sum()
            nj = randSelect(ps[0:ni]) # where to connect to
            deg[nj] += 1
            deg[ni] += 1
            net.append( (ni,nj) ) # connect ni to nj
    for i in xrange(len(net)):
        if np.random.random() < 0.5:
            net[i] = (net[i][0]+1, net[i][1]+1) # if the edge dirs are not
        else:                                   # randomized, the ant solutions
            net[i] = (net[i][1]+1, net[i][0]+1) # are singular... and fail.
    return(net)


        
def BAModel_Ver1(nodes):
    network = []
    degs = np.zeros(nodes)     # the degree of each node
    # start with a connected set of nodes
    network.append( (0,1) )
    degs[0] = 1
    degs[1] = 1
    for nodei in range(2,nodes):
        # we add a node starting at node 2
        for nodej in range(0,nodei):
            pj = float(degs[nodej])/float(sum(degs))
            if np.random.random() < pj:
                network.append( (nodei, nodej) )
                degs[nodei] += 1
                degs[nodej] += 1
    return(network)
                

    
def writeNetwork(net, fileout, wts):
    fout = open(fileout,'w')
    for ni in net:
        if wts == "T":
            wi = np.random.random()
            fout.write(str(ni[0]) + "\t" + str(ni[1]) + "\t" + str(wi) + "\n")
        else:
            fout.write(str(ni[0]) + "\t" + str(ni[1]) + "\n")
    fout.close()
