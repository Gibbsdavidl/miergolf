
from scipy.sparse import csc_matrix


def buildGraph(s):
    # read in the network #
    graphText = open(s["graphfile"],'r').read().strip().split("\n")
    edgeText  = map(lambda x: x.split(), graphText)

    # for the ants # 
    nodeDict  = dict()
    i = 0
    for node1, node2, wt in edgeText:
        nodeDict[i] = (node1,node2,float(wt),0.5,0.0)
        i+=1

    # then we search for edge-to-edge relationships
    data = []; row = []; col = [];
    for i in nodeDict.keys():
        for j in nodeDict.keys():
            nodeA = nodeDict[i]
            nodeB = nodeDict[j]
            if nodeA[1] == nodeB[0]: # nodeA leads to nodeB
                data.append(nodeA[2]*nodeB[2])
                row.append(i)
                col.append(j)

    # build the sparse matrix
    n = len(nodeDict)
    sparseMat   = csc_matrix( (data,(row,col)), shape=(n,n) )
    return((nodeDict,sparseMat))


