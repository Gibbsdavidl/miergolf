#-------------------------------------------------------------------------------------------------------------------------
#-- Influence Maximization Problem                                          ----------------------------------------------
#-- for problems on biological pathway graphs 
#-------------------------------------------------------------------------------------------------------------------------

#-- Input:
#-- config file and directed network with edge weights --

#-- Output:
#-- size K subset of nodes with the greatest influence on graph.
#-------------------------------------------------------------------------------------------------------------------------


"""
miergolf

Influence Maximization Problem by Diffusion models and Ant Optimization.

run it by: python main.py config_file graph_file
"""

import sys
import getopt
import multiprocessing as mp
from programState import *
from graph import *
from optim import *

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
    if len(args) == 0:
            print __doc__
            sys.exit(0)
    s = initProgramState(args)
    cpus = s["cpus"]
    outf = s["outf"]
    pool = mp.Pool(cpus)
    (nodes, sparseMat) = buildGraph(s)
    s2 = optimize (pool, s, nodes, sparseMat)
    pool.close()
    printResults(s2,nodes,outf)


def printResults(state,nodes,outfile):
    fout = open(outfile,'w')
    wts = weightsum(nodes, state["bestEver"][2])
    #print ("Config\tGraph\tType\tOptimTo\tMode\tNodes\tAnts\tTx\tRx\tDamp\tLocal\tScore\tTouch\tSoln\tWt\tNodes")
    fout.write ( str(state["config"]) +"\t"+
            str(state["graphfile"]) +"\t"+
            str(state["lineGraph"]) +"\t"+
            str(state["opton"]) +"\t"+
            str(state["mode"]) +"\t"+
            str(len(nodes)) + "\t"+
            str(state["ants"]) + "\t"+
            str(state["tx"]) +"\t"+
            str(state["rx"]) +"\t"+
            str(state["damp"]) +"\t"+
            str(state["local"]) +"\t"+
            str(state["bestEver"][0]) +"\t"+
            str(state["bestEver"][1]) +"\t"+
            str(state["bestEver"][2]) +"\t"+
            str(wts) +"\t"+
            ":".join([nodes[i][0]+"-"+nodes[i][1] for i in state["bestEver"][2]] )+"\n"
     )


def printGraph(flag, sparseMat):
    if flag == 1:
        pprint(sparseMat)
        sys.exit(0)
    else:
        return()


if __name__ == "__main__":
    main()
