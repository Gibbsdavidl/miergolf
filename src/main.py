#-------------------------------------------------------------------------------------------------------------------------
#-- Influence Maximization Problem using hypercube-minmax ant optimization ----------------------------------------------
#-- for problems on biological pathway graphs 
#-------------------------------------------------------------------------------------------------------------------------

#-- Input:
#-- config file and directed pathways with edge weights --

#-- Output:
#-- The program state
#-- approximate best subset of size K, with the greatest influence on graph, by diffusion models.
#-------------------------------------------------------------------------------------------------------------------------


"""
mipdao

Maximization of Influence Problem by Diffusion models and Ant Optimization.

run it by: python main.py config_file graph_file
"""

import sys
import getopt
import multiprocessing as mp

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
    s = initProgramState (args)
    cpus = s["cpus"]
    pool = mp.Pool(cpus)
    (nodes, sparseMat) = buildGraph (s)
    s2 = optimize (pool, s, nodes, sparseMat)
    pool.close()
    printResults(s2,nodes)


def printResults(s,nodes):
    print "Graph Type: " + str(s["lineGraph"])
    print "Iterations: " + str(s["iters"])
    print "Solution  : " + str(s["bestEver"])
    print "Selected Edges: "
    for (i,k) in nodes.items():
        if i in s["bestEver"][1]:
            print str(i) +"   " + str(k)


def printGraph(flag, sparseMat):
    if flag == 1:
        pprint(sparseMat)
        sys.exit(0)
    else:
        return()


if __name__ == "__main__":
    from programState import *
    from graph import *
    from antopt import *
    main()
