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

This is to derive the distribution of scores for random solutions.

run it by: python permutes.py config_file graph_file
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
    (s2,res0) = permute(pool, s, nodes, sparseMat)
    pool.close()
    printResults(s2,res0)


def printResults(s,res0):
    for p in res0:
        print str(p[0]) + "\t" + str(p[1][0]) + "\t" + str(p[1][1])
            

if __name__ == "__main__":
    from programState import *
    from graph import *
    from randomsolns import *
    main()


