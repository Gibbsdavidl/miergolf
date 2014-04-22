
import sys
import scipy.sparse as sp
import scipy.sparse.linalg as lin
import numpy as np
from bisect import bisect
import random
import itertools as it

# Hypercube MinMax Ant Optimization #
def optimize (pool, s, nodes, sparseMat):

    # for each new convergence #
    for run in xrange(s["runs"]):
        print "RUN: " + str(run) 

        # prepare for this reset-run #
        nodes = resetNodes(nodes)
        s     = resetState(s)
        iters = 0; maxIters = 100000;
        
        while iters < maxIters and s["c"] > s["ct"]:

            # generate the probability distribution # 
            ps = genProbs(s,nodes)

            # for each new ant, generate a solution
            solns = generateSolutions(pool, s, np.copy(ps))
                            
            # score each solution, choose the best #
            scores = computeScores(pool, s, solns, sparseMat)
            idx    = maxscore(scores)
            
            # perform local optimization ... if set in s["local"] #
            (bestSoln,bestScore) = parLocalOpt(pool, s, solns[idx], scores[idx], nodes, sparseMat) 
             
            # update the best/resetbest/iterbests #
            s = updateSolutionSet(s, bestSoln, bestScore, iters)

            # update the pheromones #
            nodes = updatePheromones(s, nodes)

            # check for convergence #
            s = checkConvergence(s, nodes)
            iters += 1
            s["iters"] = iters
        #(f,h) = scoreMats(s["bestEver"][2], s, sparseMat)
        #print "F"
        #print f
        #print "H"
        #print h
    return(s)


def resetNodes(nodes):
    for k in nodes.keys():
        (a,b,c,d,e) = nodes[k]
        nodes[k] = (a,b,c,0.5,e)
    return(nodes)


def resetState(s):
    s["bestRest"] = (0.0,0.0,[])
    s["bestIter"] = (0.0,0.0,[])
    s["c"] = 1.0
    return(s)

    
def genProbs(s,nodes):
    ps = np.zeros(len(nodes))
    for k in nodes.keys():
        ps[k] = (pow(nodes[k][2], s["alph"]) * pow(nodes[k][3], s["beta"]))
    return(ps)


def generateSolutions(pool, s, ps):
    antdat = it.izip(xrange(s["ants"]), it.repeat(s, s["ants"]), it.repeat(ps, s["ants"]))
    solns = pool.map(genSoln, antdat)
    return(solns)

    
def genSoln((i, s, ps)):
    soln = []
    r = random.Random()
    r.jumpahead(int(1000*r.random()))
    for ki in xrange(int(s["k"])):
        ps = ps/(sum(ps)) # after removing one ... renorm the probs
        cs = ps.cumsum()
        solni = bisect(cs,r.random())
        soln.append(solni)
        ps[solni] = 0
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
        

def scoreMats(solns, s, smat):
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
    f = lin.spsolve(lap, pts) 
    h = lin.spsolve(lap_t, pst_t)
    return((f,h))


def subMatrix(rows, cols, A):
    return(A.tocsr()[rows,:].tocsc()[:,cols])


def maxscore(scores):
    best = 0
    idx = 0
    m = -1
    for (a,b) in scores:
        if a > m:
            m = a
            best = idx
        idx += 1
    return(best)

    
def parLocalOpt(pool, s, bestSoln,
                (bestScore,bestTouch),
                nodes, sparseMat):
    if s["local"] == -1:   # none
        # go back! go back!
        return (bestSoln,(bestScore,bestTouch))
    elif s["local"] == 0:  # full optimization, slow!
         # build a list of possible solns
        newSolnList = []
        trySoln = []
        for i in xrange(len(bestSoln)):
            for ni in xrange(len(nodes)):
                if ni not in bestSoln:
                    trySoln = list(bestSoln)
                    trySoln[i] = ni
                    newSolnList.append(list(trySoln))
    else:
        # compute a smaller local area around the result.
        newSolnList = []
        trySolnUp = []; trySolnDn = [];
        for i in xrange(len(bestSoln)):
            for ni in range(1,s["local"]):  # the local search range
                trySolnUp = trySolnDn = list(bestSoln)
                trySolnUp[i] = max((bestSoln[i]+ni),len(nodes))
                trySolnDn[i] = min(0, (bestSoln[i]-ni))
                newSolnList.append(list(trySolnUp))
                newSolnList.append(list(trySolnDn))
                    
    # then score each potential solution
    scores = computeScores(pool, s, newSolnList, sparseMat)
    idx    = maxscore(scores)

    # anything better?
    if scores[idx][0] > bestScore:
        return( (newSolnList[idx], scores[idx]) )
    else:
        return (bestSoln,(bestScore,bestTouch))
        

def updateSolutionSet(s, bestSoln, (bestScore,bestTouch), iters):
    if s["opton"] == "score":
        if bestScore > s["bestEver"][0]:
            s["bestEver"] = (bestScore, bestTouch, bestSoln)
        if bestScore > s["bestRest"][0]:
            s["bestRest"] = (bestScore, bestTouch, bestSoln)
        if bestScore > s["bestIter"][0]:
            s["bestIter"] = (bestScore, bestTouch, bestSoln)
    elif s["opton"] == "touch": 
        if bestTouch > s["bestEver"][1]:
            s["bestEver"] = (bestScore, bestTouch, bestSoln)
        if bestTouch > s["bestRest"][1]:
            s["bestRest"] = (bestScore, bestTouch, bestSoln)
        if bestTouch > s["bestIter"][1]:
            s["bestIter"] = (bestScore, bestTouch, bestSoln)
    elif s["opton"] == "combo":
        if bestTouch*bestScore > s["bestEver"][1]*s["bestEver"][0]:
            s["bestEver"] = (bestScore, bestTouch, bestSoln)
        if bestTouch*bestScore > s["bestRest"][1]*s["bestRest"][0]:
            s["bestRest"] = (bestScore, bestTouch, bestSoln)
        if bestTouch*bestScore > s["bestIter"][1]*s["bestIter"][0]:
            s["bestIter"] = (bestScore, bestTouch, bestSoln)
    else:
        print "Error! config option 'opton' must be score, touch, or combo"
        sys.exit(1)
    return(s)


def updatePheromones(s, nodes):
    (iterp, restp, bestp) = pheroProportions(s)
    restartSoln = s["bestRest"][2]
    iterateSoln = s["bestIter"][2]
    bestSoln    = s["bestEver"][2]
    for k in nodes.keys():
        (a,b,w,p,ch) = nodes[k]
        inRest = int(k in restartSoln)
        inIter = int(k in iterateSoln)
        inBest = int(k in bestSoln)
        deposit = inIter * iterp + inRest * restp + inBest * bestp
        p2 = bounded(p + s["evap"]*(deposit - p))
        nodes[k] = (a,b,w,p2,(ch+(inIter+inRest+inBest)/3.0))
    return(nodes)


def pheroProportions(s):
    # return proportions of solutions to use
    # (iteration, restart, best)
    sc = s["c"]
    if sc > 0.8:
        x= (1.0, 0.0, 0.0) #- just started out, use iteration best
    elif sc >= 0.6 and sc < 0.8:
        x= (0.6669, 0.3331, 0.0)
    elif sc >= 0.4 and sc < 0.6:
        x= (0.3331, 0.6669, 0.0) # nearing the end - move to restart 
    elif sc >= 0.2 and sc > 0.4:
        x= (0.0, 0.6669, 0.3331)
    elif sc >= 0.1 and sc < 0.2:
        x= (0.0, 0.3331, 0.6669) # nearing the end - move to best ever
    else:
        x = (0.0,0.0,1.0)
    return(x)


def bounded(x):
    if x < 0.001:
        return(0.001)
    elif x > 0.999:
        return(0.999)
    else:
        return(x)


def checkConvergence(s, nodes):
    normfactor = (0.999-0.001) * len(nodes)
    ps = [p for (i,(a,b,c,p,e)) in nodes.items()]
    convergence = 1.0 - 2.0 * (( sum(map(convNum, ps)) / normfactor ) - 0.5 )
    pToKRatio = (sum (ps)) / ((0.999-0.001) * s["k"])
    s["c"] = convergence
    s["pTokRatio"] = pToKRatio
    return(s)

def divn(x, n):
    round(x/n)

def convNum(x):
    return(max(0.999-x, x-0.001))


