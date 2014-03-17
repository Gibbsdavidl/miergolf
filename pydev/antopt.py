
import sys
import scipy.sparse as sp
import scipy.sparse.linalg as lin
import numpy as np
from bisect import bisect
from random import random
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
            solns = [genSoln(i,s,np.copy(ps)) for i in xrange(s["ants"])]
                
            # score each solution, choose the best #
            scores = computeScores(pool, s, solns, sparseMat)
            idx    = scores.index(max(scores))

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

    return(s)
            

def resetNodes(nodes):
    for k in nodes.keys():
        (a,b,c,d,e) = nodes[k]
        nodes[k] = (a,b,c,0.5,e)
    return(nodes)


def resetState(s):
    s["bestRest"] = (0.0,[])
    s["bestIter"] = (0.0,[])
    s["c"] = 1.0
    return(s)

    
def genProbs(s,nodes):
    ps = np.zeros(len(nodes))
    for k in nodes.keys():
        ps[k] = (pow(nodes[k][2], s["alph"]) * pow(nodes[k][3], s["beta"]))
    return(ps)


def genSoln(i, s, ps):
    soln = []
    for ki in xrange(int(s["k"])):
        ps = ps/(sum(ps)) # after removing one ... renorm the probs
        cs = ps.cumsum()
        solni = bisect(cs,random())
        soln.append(solni)
        ps[solni] = 0
    return(soln)


def computeScores(pool, s, solns, sparseMat):
    scoreDat = it.izip(solns, it.repeat(sparseMat, len(solns)))
    scores = pool.map(scoreSoln, scoreDat)
    return(scores)


def scoreSoln( (ss,smat) ):
    # ss    -- the solutions set S
    # smat  -- nxn sparse matrix
    # ts    -- the set T
    n = (smat.shape[0]-1)
    ts = [i for i in xrange(n) if i not in ss]
    idn = sp.eye(len(ts))
    pst = subMatrix(ss, ts, smat)
    pts = subMatrix(ts, ss, smat)
    ptt = subMatrix(ts, ts, smat)
    lap = sp.csc_matrix(idn-ptt)
    pst_t = sp.csc_matrix(pst.transpose())
    lap_t = sp.csc_matrix(lap.transpose())    
    f = lin.spsolve(lap, pts) 
    h = lin.spsolve(lap_t, pst_t)
    return(f.sum() + h.sum())


def subMatrix(rows, cols, A):
    return(A.tocsr()[rows,:].tocsc()[:,cols])


def parLocalOpt(pool, s, bestSoln, bestScore, nodes, sparseMat):
    if s["local"] != 1:
        # go back! go back!
        return (bestSoln, bestScore)
    else:
         # build a list of possible solns
        newSolnList = []
        trySoln = []
        for i in xrange(len(bestSoln)):
            for ni in xrange(len(nodes)):
                if ni not in bestSoln:
                    trySoln = list(bestSoln)
                    trySoln[i] = ni
                    newSolnList.append(list(trySoln))

        # then score each potential solution
        scores = computeScores(pool, s, newSolnList, sparseMat)
        idx    = scores.index(max(scores))

        # anything better?
        if scores[idx] > bestScore:
            #print ("   local improvement found.. " 
            #       + str(bestSoln) +"   "+ str(trySoln) 
            #       + "  old score: " + str(bestScore) 
            #       + "  new score: " + str(scores[idx])) 
            return( (newSolnList[idx], scores[idx]) )
        else:
            return (bestSoln, bestScore)


def updateSolutionSet(s, bestSoln, bestScore, iters):
    if bestScore > s["bestEver"][0]:
        s["bestEver"] = (bestScore, bestSoln)
    if bestScore > s["bestRest"][0]:
        s["bestRest"] = (bestScore, bestSoln)
    if bestScore > s["bestIter"][0]:
        s["bestIter"] = (bestScore, bestSoln)
    return(s)


def updatePheromones(s, nodes):
    (iterp, restp) = pheroProportions(s)
    restartSoln = s["bestRest"][1]
    iterateSoln = s["bestIter"][1]
    for k in nodes.keys():
        (a,b,w,p,ch) = nodes[k]
        inRest = int(k in restartSoln)
        inIter = int(k in iterateSoln)
        deposit = inIter * iterp + inRest * restp
        p2 = bounded(p + s["evap"]*(deposit - p))
        nodes[k] = (a,b,w,p2,(ch+1))
    return(nodes)


def pheroProportions(s):
    # return proportions of solutions to use
    # (iteration, restart)
    sc = s["c"]
    if sc > 0.8:
        x= (1.0, 0.0) #- just started out, use iteration best
    elif sc >= 0.6 and sc < 0.8:
        x= (0.6669, 0.3331)
    elif sc >= 0.4 and sc < 0.6:
        x= (0.3331, 0.6669) 
    elif sc < 0.4:
        x= (0.0, 1.0) # nearing the end - just use restart best
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
    convergence = 1.0 - ( 2.0 * ( sum(map(convNum, ps)) / normfactor ) - 0.5)
    pToKRatio = (sum (ps)) / ((0.999-0.001) * s["k"])
    s["c"] = convergence
    s["pTokRatio"] = pToKRatio
    return(s)


def convNum(x):
    return(max(0.999-x, x-0.001))


