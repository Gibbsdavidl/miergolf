from numpy import array
from scipy.sparse import csc_matrix, lil_matrix, coo_matrix
from sys import exit
from pprint import pprint


# *** IDEA! *** #
# randomize the order of edges, so that the local search
# searches a random set of edges ... 

def buildGraph(s):
    # read in the network #
    graphText = open(s["graphfile"],'r').read().strip().split("\n")
    edgeText  = map(lambda x: x.split(), graphText)

    # for the ants # 
    edgeDict  = dict()
    nodeDict  = dict()
    matrixLookup = dict()
    nodeVect  = set()
    i = 0
    for node1, node2, wt in edgeText:
        edgeDict[i] = (node1,node2,float(wt),0.5,0.0)
        if node1 in nodeDict:
            nodeDict[node1] += [float(wt)]  # list of the outgoing edges
        else:
            nodeDict[node1] = [float(wt)]
        nodeVect.add(node1.strip()); nodeVect.add(node2.strip())
        i+=1

    nodeVect = list(nodeVect)
    i = 0
    for ni in nodeVect:
        matrixLookup[ni] = i
        i+=1

    data = []; row = []; col = [];

    if "line" in s["lineGraph"]:
        # then we search for edge-to-edge relationships
        for i in edgeDict.keys():
            for j in edgeDict.keys():
                nodeA = edgeDict[i]
                nodeB = edgeDict[j]
                if nodeA[1] == nodeB[0]: # nodeA leads to nodeB
                    data.append(nodeA[2]+nodeB[2])
                    row.append(i)
                    col.append(j)
        n = len(edgeDict)
    else:
        # then we use the regular graph structure #
        # the edgeDict becomes mean wts on nodes
        # and the matrix holds the regular edge weights
        i = 0
        newEdgeDict = dict()
        for nodei in nodeVect:  # the data structure used for sampling nodes
            if nodei in nodeDict:
                meanWt = sum(nodeDict[nodei]) / len(nodeDict[nodei])
                newEdgeDict[i] = (nodei, nodei, meanWt, 0.5, 0.0)
            else:
                newEdgeDict[i] = (nodei, nodei, 0.0, 0.5, 0.0)
            i += 1
        for edge in edgeDict.items():  # the matrix.
            i = matrixLookup[edge[1][0]]
            j = matrixLookup[edge[1][1]]
            data.append(edge[1][2])
            row.append(i)
            col.append(j)
        n = len(nodeVect)
        edgeDict = newEdgeDict

      # NEED TO BUILD A NODE CENTRIC DICT FOR PROBS #
        # maybe sum up the edge weights around a node,
        # and take the average, put that on the node.

    # build the sparse matrix
    sparseMat = csc_matrix( (data,(row,col)), shape=(n,n) )
    sparseLil = lil_matrix(sparseMat)

    # Then to normalize the rows ... and to apply dampening to the 
    # probabilities.
    for r in xrange(n):
        rowsum = (sparseLil[r,:]).sum()
        if rowsum != 0.0:
            for c in xrange(n):
                sparseLil[r,c] = (sparseLil[r,c]/rowsum) * s["damp"]

    return((edgeDict,csc_matrix(sparseLil)))


def weightsum(nodes,tup):
    totwt = 0.0
    for ti in tup:
        totwt += nodes[ti][2]
    return(totwt)
