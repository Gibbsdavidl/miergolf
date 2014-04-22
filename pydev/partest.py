import sys
import numpy as np
from bisect import bisect
import random
import itertools as it
import getopt
import multiprocessing as mp


def generateSolutions(pool):
    antdat = it.izip(xrange(10), xrange(10))
    solns = pool.map(genSoln, antdat)
    return(solns)

    
def genSoln((x, y)):
    soln = []
    r = random.Random()
    r.jumpahead(int(1000*r.random()))
    for ki in xrange(10):
        soln.append(r.random())
    return(sum(soln))


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
    cpus = 4
    pool = mp.Pool(cpus)
    x = generateSolutions(pool)
    pool.close()
    printResults(x)


def printResults(x):
    for xi in x:
        print xi


if __name__ == "__main__":
    main()
