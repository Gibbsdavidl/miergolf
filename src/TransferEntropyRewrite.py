
import numpy as np
from copy import copy, deepcopy
import multiprocessing as mp
from random import uniform, shuffle
from math import log
import itertools as it

# In this work, I am computing information theoretic values
# by, first, discretizing expression values into a given
# number of bins. Using those bins, the probability of a given
# interval is computed, and the joint probability over time 
# can also be computed (given two time series).

# Want  P(X_t+1, Y_k1, X_k2) * log (P(X_t+1,Y_k1,X_k2)*P(X_t+1)) / (P(X_t+1, X_k2)*P(X_k2,Y_K1))

# just get the joint, then get the others by marginalization

# parameters:
# ym: the markov order for Y = let it be 1
# xm: the markov order for x = let it be 1
# yk: the time delay for y
# xk: the time delay for x
# b : the number of bins

def safelog2(x):
    if (x > 0):
        return(log(x,2))
    else:
        return(0.0)


def genInd(x,xorg,k):
    for b in xrange(k):
        y = []
        for xi in x:  # [0]
            for xorgi in xorg:  # 1
                y.append([xorgi] + xi)
        x = copy(y)
    return(x)

        
def generateIndices(k, bins):
    # return all tuples of length k, 
    # with entries of 0 to (bins-1)
    x = [[a] for a in xrange(bins)]
    xorg = range(0,bins)
    l = genInd(x, xorg, k)
    tl = map(lambda z: tuple(z), l)
    return(tl)
    

def movingAverageSmoothing(x, wsize):
    # x: a numeric vector
    # wsize: the window size
    # repsize: 41 -- there's 41 time measurements 
    xint = []
    for i in range(0, (len(x)-3)):
        xint.append(sum( [x[j] for j in range(i,i+wsize)] )/wsize)
    return(xint)


def jointProb (x, y, yk, bins):  # just to lag Y #
    xint = movingAverageSmoothing(x,3)
    xidx = np.histogram(xint,bins)[1][range(0,bins)] # just need the left edges of bins
    xi   = np.digitize(xint, xidx)-1  # for each x_i, what bin? zero indexed
    yint = movingAverageSmoothing(y,3)
    yidx = np.histogram(yint,bins)[1][range(0,bins)] # just need the left edges of bins
    yi   = np.digitize(yint, yidx)-1  # for each y_i, what bin? zero indeyed
    p    = np.zeros(np.repeat(bins, 3)) # going to be a cube, 3D

    # generate the joint counts over time points (orders xk, yk) #
    for i in range(yk, len(xi)):
        kidx = range(i, i-2, -1) # i and the last K entries in time, for X
        lidx = i-yk # i and the last K entries in time, for Y
        tx = [xi[z] for z in kidx]
        ty = [yi[lidx]]
        p[tuple(tx+ty)] += 1 # count that state sequence
    return(p/(p.sum()))


def transferEntropy(a,b,yk,bins):
    joint = jointProb(a,b,yk,bins)
    marginalX = joint.sum(0).sum(1) # marginalize out the past events
    marginalXX = joint.sum(2) # marginalize out the Y
    marginalXY = joint.sum(0) # marginalzie out the current time step
    bigSum = 0    
    # then for each "column" of states (say t-1, t-2) sum it
    for tup in generateIndices(2, bins):
        bigSum += joint[tup] * safelog2((joint[tup]*marginalX[tup[0]]) / (marginalXX[tup[0:2]]*marginalXY[tup[1:3]]))
    return(bigSum)


def te((a,b,k,bins)):
    return(transferEntropy(a,b,k,bins))


def transferEntropyPermutationA(a,b,k,bins,p,cpus):
    pool = mp.Pool(cpus)
    observedTE = transferEntropy(a,b,k,bins)
    permutedList = it.repeat( (a,b,k,bins), p)
    permutedTE = pool.map(shufflete, permutedList)
    pvalueCount = sum( [1.0 for pte in permutedTE if pte > observedTE] )
    temean = np.mean(permutedTE)
    testd = np.std(permutedTE)
    zscore = (observedTE - temean)/testd
    pool.close()
    return([observedTE, temean, testd, zscore, pvalueCount, (pvalueCount/p)])


def tePermutationAScanner(a,b,k,bins,p,win,end,cpus):
    pool = mp.Pool(cpus)
    scan = np.array([0.0 for i in xrange(win)])
    for i in xrange(win):
        l = (end-win)+i
        aTob = transferEntropy(b[i:l],a[i:l],k,bins)
        bToa = transferEntropy(a[i:l],b[i:l],k,bins)
        scan[i] = (abs(aTob - bToa))
    idx = scan.argmax() # this is the window offset with the maximum distance
    ldx = (end-win)+idx
    observedTE = transferEntropy(a[idx:ldx],b[idx:ldx],k,bins)
    permutedList = it.repeat( (a[idx:ldx],b[idx:ldx],k,bins), p)
    permutedTE = pool.map(shufflete, permutedList)
    pvalueCount = sum( [1.0 for pte in permutedTE if pte > observedTE] )
    temean = np.mean(permutedTE)
    testd = np.std(permutedTE)
    zscore = (observedTE - temean)/testd
    pool.close()
    return([idx, ldx, observedTE, temean, testd, zscore, pvalueCount, (pvalueCount/p)])


def tescan(filename, fileout, endtime, k, bins, reps, win, cpus):
    dat = open(filename,'r').read().strip().split("\n")
    fout = open(fileout,'w')
    dats = map(lambda x: x.split("\t"), dat)
    marker = 0
    for i in xrange(len(dat)):
        for j in xrange(len(dat)):
            id1 = dats[i][0]
            id2 = dats[j][0]
            a = map(float,dats[i][1:endtime])
            b = map(float,dats[j][1:endtime])
            res0 = tePermutationAScanner(a,b,k,bins,reps,win,endtime,cpus)
            fout.write( str(i)+"\t"+ str(j) +"\t"+ id1 +"\t"+ id2 + "\t" + "\t".join(map(str,res0))+"\n")
            if marker % 997 == 0:
                print str(marker) +"\t"+ str(i)+"\t"+ str(j) +"\t"+ id1 +"\t"+ id2 + "\t" + "\t".join(map(str,res0))
            marker +=1
    fout.close()

