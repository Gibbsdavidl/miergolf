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
    nodeVect  = []
    i = 0
    for node1, node2, wt in edgeText:
        edgeDict[i] = (node1,node2,float(wt),0.5,0.0)
        nodeVect.append(node1); nodeVect.append(node2)
        i+=1
    nodeVect = array(nodeVect)

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
        for i in edgeDict.keys():
            edgeA = edgeDict[i]
            a = (nodeVect == edgeA[0]).argmax()
            b = (nodeVect == edgeA[1]).argmax()
            data.append(edgeA[2])
            row.append(a)
            col.append(b)
        n = len(nodeVect)
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



