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

Maximum Influence Problem by Diffusion models and Ant Optimization.

run it by: python main.py config_file graph_file
"""

import sys
import getopt

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
    (nodes, sparseMat) = buildGraph (s)
    s2 = optimize (s, nodes, sparseMat)
    printResults(s2,nodes)


def printResults(s,nodes):
    print s["bestEver"]
    for (i,k) in nodes.items():
        if i in s["bestEver"][1]:
            print str(i) +"   " + str(k)


if __name__ == "__main__":
    from programState import *
    from graph import *
    from antopt import *
    main()
