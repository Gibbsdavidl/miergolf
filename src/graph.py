from numpy import array
from scipy.sparse import csc_matrix, lil_matrix, coo_matrix
from sys import exit
from pprint import pprint

def buildGraph(s):
    # read in the network #
    graphText = open(s["graphfile"],'r').read().strip().split("\n")
    edgeText  = map(lambda x: x.split(), graphText)

    # for the ants # 
    nodeDict  = dict()
    nodeVect  = []
    i = 0
    for node1, node2, wt in edgeText:
        nodeDict[i] = (node1,node2,float(wt),0.5,0.0)
        nodeVect.append(node1); nodeVect.append(node2)
        i+=1
    nodeVect = array(nodeVect)

    data = []; row = []; col = [];
    if "line" in s["lineGraph"]:
        # then we search for edge-to-edge relationships
        for i in nodeDict.keys():
            for j in nodeDict.keys():
                nodeA = nodeDict[i]
                nodeB = nodeDict[j]
                if nodeA[1] == nodeB[0]: # nodeA leads to nodeB
                    data.append(nodeA[2]*nodeB[2])
                    row.append(i)
                    col.append(j)
        n = len(nodeDict)
    else:
        # then we use the regular graph structure #
        for i in nodeDict.keys():
            edgeA = nodeDict[i]
            a = (nodeVect == edgeA[0]).argmax()
            b = (nodeVect == edgeA[1]).argmax()
            data.append(edgeA[2])
            row.append(a)
            col.append(b)
        n = len(nodeVect)

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

    return((nodeDict,csc_matrix(sparseLil)))



