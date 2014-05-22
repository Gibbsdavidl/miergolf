
import sys
import scipy.sparse as sp
import scipy.sparse.linalg as lin
import numpy as np
from bisect import bisect
import random
import itertools as it

# Distribution of random solutions #
def permute (pool, s, nodes, sparseMat):

    # for the number of permutations, generate a solution
    solns = generateRandomSolutions(pool, int(s["k"]), int(s["permutes"]), len(nodes))
                            
    # score each solution #
    scores = computeScores(pool, s, solns, sparseMat)

    # the list to print out.
    res0 = buildReturnList (solns, scores)
    
    return((s, res0))



def generateRandomSolutions(pool, k, p, n):
    antdat = it.izip(xrange(p), it.repeat(k,p), it.repeat(n,p))
    solns = pool.map(genSoln, antdat)
    return(solns)

    
def genSoln((i, k, n)):
    soln = []
    r = random.Random()
    r.jumpahead(int(1000*r.random()))
    for ki in xrange(k):
        soln.append(int(r.random()*n))
    return(soln)


def computeScores(pool, s, solns, sparseMat):
    scoreDat = it.izip(solns, it.repeat( (s,sparseMat), len(solns)))
    scores = pool.map(scoreSoln, scoreDat)
    return(scores)


def scoreSoln( (solns, (s,smat)) ):
    # ss    -- the solutions set S
    # smat  -- nxn sparse matrix
    # ts    -- the set T
    n = (smat.shape[0]-1)
    ts = [i for i in xrange(n) if i not in solns]
    idn = sp.eye(len(ts))
    pst = subMatrix(solns, ts, smat)
    pts = subMatrix(ts, solns, smat)
    ptt = subMatrix(ts, ts, smat)
    lap = sp.csc_matrix(idn-ptt)
    pst_t = sp.csc_matrix(pst.transpose())
    lap_t = sp.csc_matrix(lap.transpose())    
    if s["mode"] == "both":
        f = lin.spsolve(lap, pts) 
        h = lin.spsolve(lap_t, pst_t)
        if type(f) == type(np.array([])): # came back as an array
            ftouch = sum(f > s["tx"])
            htouch = sum(h > s["rx"])
        else: # came back as a sparse matrix
            ftouch = sum(f.toarray().flatten() > s["tx"])
            htouch = sum(h.toarray().flatten() > s["rx"])
            #         best score  ... best touch       #
            return((f.sum()+h.sum(), ftouch+htouch))
    elif s["mode"] == "tx":
        f = lin.spsolve(lap, pts) 
        if type(f) == type(np.array([])): # came back as an array
            ftouch = sum(f > s["tx"])
        else: # came back as a sparse matrix
            ftouch = sum(f.toarray().flatten() > s["tx"])
            #         best score  ... best touch       #
            return((f.sum(), ftouch))
    elif s["mode"] == "rx":
        h = lin.spsolve(lap_t, pst_t)
        if type(h) == type(np.array([])): # came back as an array
            htouch = sum(h > s["rx"])
        else: # came back as a sparse matrix
            htouch = sum(h.toarray().flatten() > s["rx"])
            #         best score  ... best touch       #
            return((h.sum(), htouch))
    else:
        print "Error! mode must be rx, tx, or both."
        sys.exit(1)
        


def subMatrix(rows, cols, A):
    return(A.tocsr()[rows,:].tocsc()[:,cols])


def buildReturnList (solns, scores):
    solnList = []
    for (a, b) in it.izip(solns, scores):
        solnList.append([a,b])
    return solnList
