
import sys
import scipy.sparse as sp
import scipy.sparse.linalg as lin
import numpy as np
from bisect import bisect
import itertools as it
from diffusion import *

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

            # for each new ant, generate a solution, local search, and score
            (bestSoln,bestScore) = antWork(pool, s, np.copy(ps), sparseMat, nodes)

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
    s["bestRest"] = (0.0,0.0,[])
    s["bestIter"] = (0.0,0.0,[])
    s["c"] = 1.0
    return(s)


def genProbs(s,nodes):
    # generate the probablities for selecting each node
    ps = np.zeros(len(nodes))
    for k in nodes.keys():
        ps[k] = (pow(nodes[k][2], s["alph"]) * pow(nodes[k][3], s["beta"]))
    return(ps)


def antWork(pool, s, ps, sparseMat, nodes):
    # for the solns
    nants = s["ants"]
    seeds = [np.random.randint(1000000) for i in xrange(nants)] # seed each thread
    # generate the list of data for each ant
    antdat = (it.izip(xrange(nants), it.repeat(s, nants), it.repeat(ps, nants), 
                     it.repeat(sparseMat, nants), seeds, it.repeat(nodes, nants)))
    # send it to a pooled fun party
    solns = pool.map(poolParty, antdat)
    # return the best
    return(scoreMax(solns, s))


def poolParty( (i, s, ps, sparseMat, seedi, nodes) ):
    # start with a possible solution
    np.random.seed(seedi)
    soln = genSoln(i,s,np.copy(ps))
    # then do a local search
    (soln2, score) = localSearch(s, np.copy(ps), soln, sparseMat, nodes)
    # return the best
    return(soln2, score)

    
def genSoln(i, s, ps):
    # generate a solution, probabilistically for each ant.    
    soln = []
    for ki in xrange(int(s["k"])):
        ps = ps/(sum(ps)) # after removing one ... renorm the probs
        cs = ps.cumsum()
        solni = bisect(cs,np.random.random()) # should skip over the 0'ed ones...
        soln.append(solni)
        ps[solni] = 0
    return(soln)


def localSearch(s, ps, bestSoln, sparseMat, nodes):
    (bestScore,bestTouch) = scoreSoln(bestSoln, s, sparseMat, nodes)
    # soln now doesn't include the edge weights
    if s["local"] == -1:   # none
        # go back! go back!
        return (bestSoln,(bestScore,bestTouch))
    else:
        # hill climbing for a certain number of steps
        newSoln = list(bestSoln)
        newScore = bestScore
        newTouch = bestTouch
        ps = ps/(sum(ps))
        cs = ps.cumsum()
        n = s["local"] # the number of tries to make
        testsoln = list(newSoln)
        for i in xrange(n):
            remr  = testsoln[np.random.randint(0,len(testsoln),1)]        # the one to remove
            solnr = [xi for xi in testsoln if xi != remr]  # fragment list
            solni = testsoln[0];                           # pick a new one, not in the list already
            while solni in testsoln:
                solni = bisect(cs,np.random.random())      # the one to add, based on ps
            testsoln = list( (solnr + [solni]) )           # the new soln list
            score    = scoreSoln(testsoln, s, sparseMat, nodes)   # score it
            if s["opton"] == "touch":
                if score[1] > newTouch:
                    newScore = score[0]          # if better: keep it
                    newTouch = score[1]
                    newSoln = list(testsoln)
                else:
                    testsoln = list(newSoln)  # else: return to previous soln                    
            elif s["opton"] == "score":
                if score[0] > newScore:
                    newScore = score[0]          # if better: keep it
                    newTouch = score[1]
                    newSoln = list(testsoln)
                else:
                    testsoln = list(newSoln)  # else: return to previous soln                    

            else:        
                print "This opton mode is not implemented yet!, use score or touch."
                sys.exit(1)
        return (newSoln, (newScore, newTouch))


def combo(a,b):
    return(a*b)


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
        if combo(bestTouch,bestScore) > combo(s["bestEver"][1],s["bestEver"][0]):
            s["bestEver"] = (bestScore, bestTouch, bestSoln)
        if combo(bestTouch,bestScore) > combo(s["bestRest"][1],s["bestRest"][0]):
            s["bestRest"] = (bestScore, bestTouch, bestSoln)
        if combo(bestTouch*bestScore) > combo(s["bestIter"][1],s["bestIter"][0]):
            s["bestIter"] = (bestScore, bestTouch, bestSoln)
    else:
        print "Update Solution Error! config option 'opton' must be score, touch, or combo"
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


