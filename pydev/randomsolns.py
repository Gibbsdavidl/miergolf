
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


def weightsum(nodes,tup):
    totwt = 0.0
    for ti in tup:
        totwt += nodes[ti][2]
    return(totwt)


def scoreSoln(soln, s, smat, nodes):
    # soln -- the solutions set S
    # smat  -- nxn sparse matrix
    # s     -- the program state
    wt = weightsum(nodes,soln)
    n = (smat.shape[0])
    ts = [i for i in xrange(n) if i not in soln] # set T
    idn = sp.eye(len(ts),len(ts))
    pst = subMatrix(soln, ts, smat) # prob of S to T
    pts = subMatrix(ts, soln, smat) # prob of T to S
    ptt = subMatrix(ts, ts, smat)    # prob of T to T
    lap = sp.csc_matrix(idn-ptt)
    pst_t = sp.csc_matrix(pst.transpose())
    lap_t = sp.csc_matrix(lap.transpose())    
    if s["mode"] == "both":
        return(scoreBoth(s,lap,pts,lap_t,pst_t,wt))
    elif s["mode"] == "tx":
        return(scoreTX(s,lap,pts,wt))
    elif s["mode"] == "rx":
        return(scoreRX(s,lap_t,pst_t,wt))
    else:
        print "ScoreSoln Error! mode must be rx, tx, or both."
        sys.exit(1)


def scoreBoth(s,lap,pts,lap_t,pst_t,wt):
    f = lin.spsolve(lap, pts) 
    h = lin.spsolve(lap_t, pst_t)
    if type(f) == type(np.array([])): # came back as an array
        fh = f+h
        score = fh.sum()
        touch = (fh > s["tx"]).sum()
    else: # came back as a sparse matrix
        fsum = np.array(f.sum(1)).flatten()
        hsum = np.array(h.sum(1)).flatten()
        fh = fsum + hsum
        touch = sum(fh > s["tx"])
        #  best score; best touch #
    return((touch+wt, touch))



def scoreRX(s,lap,pts,wt):
    f = lin.spsolve(lap, pts) 
    if type(f) == type(np.array([])): # came back as an array
        ftouch = sum(f > s["rx"])
    else: # came back as a sparse matrix
        fsum = np.array(f.sum(1)).flatten()
        ftouch = sum(fsum > s["rx"])
    return((ftouch+wt, ftouch))


def scoreTX(s, lap_t, pst_t, wt):
    h = lin.spsolve(lap_t, pst_t)
    if type(h) == type(np.array([])): # came back as an array
        htouch = sum(h > s["tx"])
    else: # came back as a sparse matrix
        hsum = np.array(h.sum(1)).flatten()
        htouch = sum(hsum > s["tx"])
    return((htouch+wt, htouch))


def subMatrix(rows, cols, A):
    return(A.tocsr()[rows,:].tocsc()[:,cols])


def buildReturnList (solns, scores):
    solnList = []
    for (a, b) in it.izip(solns, scores):
        solnList.append([a,b])
    return solnList
