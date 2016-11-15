import sys
import numpy as np
from copy import copy, deepcopy
import multiprocessing as mp
from numpy.random import shuffle, random, normal
from math import log, sqrt, exp, pi
import itertools as it
from scipy.stats import gaussian_kde
from scipy.stats import ttest_1samp
from itertools import product
import gzip

# In this work, I am computing transfer entropies
# by, first, discretizing expression values into a given
# number of bins. Using those bins, the probability of a given
# interval is computed, and the joint probability over time
# can also be computed (given two time series).

# Want P(X_t+1, X_k2, Y_k1) * log (P(X_t+1,Y_k1,X_k2)*P(X_t+1)) / (P(X_t+1, X_k2)*P(X_k2,Y_K1))

# just get the joint, then get the others by marginalization

# parameters:
# yk: the markov order for Y = let it be 1
# xk: the markov order for x = let it be 1
# yl: the time delay for y
# xl: the time delay for x
# b : the number of bins

# autoTE is
# FOR TE (Y -> X)

def safelog(x):
    if (x > 0.0):
        return(log(x,2))
    else:
        return(0.0)


def safevec(x):
    a = []; b = []; c = []
    for xi in x:
        a.append(xi[0])
        b.append(xi[1])
        c.append(xi[2])
    return([a,b,c])


def makegrid2(d,n):
    # a d-dimensional grid
    c = 0; didx = dict()
    seq = np.linspace(-1,1,n)
    grid = []
    for idx in product(seq, repeat = d):
        grid.append(idx)
        didx[idx]=c
        c+=1
    return((grid, didx))


def makegrid(n):
    # a 3d grid
    idx = 0; didx = dict()
    seq = np.linspace(-1,1,n)
    a = []; b = []; c = []
    for i in xrange(n):
        for j in xrange(n):
            for k in xrange(n):
                a.append(seq[i]);
                b.append(seq[j]);
                c.append(seq[k])
                didx[(i,j,k)]=idx
                idx+=1
    return(([a,b,c], didx))


def marginalize(n, idx, probs):
    px  = np.zeros(n)
    pxx = np.zeros([n,n])
    pxy = np.zeros([n,n])
    pxyz = np.zeros([n,n,n])
    for i in xrange(n):
        for j in xrange(n):
            for k in xrange(n):
                l = idx[(i,j,k)]
                px[j] += probs[l]
                pxx[i][j] += probs[l]
                pxy[j][k] += probs[l]
                pxyz[i][j][k] += probs[l]
    return((px,pxx,pxy,pxyz))


def computeTE(n,(px,pxx,pxy,pxyz)):
    te = 0.0
    for i in xrange(n):
        for j in xrange(n):
            for k in xrange(n):
                num = pxyz[i][j][k]*px[j]
                den = pxx[i][j]*pxy[j][k]
                if den > 0.0:
                    te += pxyz[i][j][k]*safelog(num/den)
    return(te)


def kernelTE (y, x, yl, h, n):
    # normalize x, y
    x = np.array(x); x = (x-x.mean())/max(1,(x-x.mean()).max())
    y = np.array(y); y = (y-y.mean())/max(1,(y-y.mean()).max())
    l = [] # data list  (list of all)
    for i in range(yl, len(x)):
        kidx = [i,(i-1)]  # index for x_t and x_t-1
        lidx = i-yl       # index for y_t-yl
        l.append([x[z] for z in kidx]+[y[lidx]])
    lpdf = gaussian_kde(safevec(l), h)
    (grid,idx) = makegrid(n)         # 3D grid of coordinates
    lprobs  = lpdf(grid)             # these are density estimates
    lprobs  = lprobs/sum(lprobs)     # normalize to probabiliies
    marprobs = marginalize(n, idx, lprobs)  # marginalized prob
    te = computeTE(n, marprobs)
    return(te)


def shufflete((y,x,yl,h,n)):
    permutedA  = deepcopy(y)
    shuffle(permutedA)
    return(kernelTE(permutedA,x,yl,h,n))


def autoTE (y, x, yl, n):
    l = [] # data list  (list of all)
    for i in range(yl, len(x)):
        kidx = [i,(i-1)]  # index for x_t and x_t-1
        lidx = i-yl       # index for y_t-yl
        l.append([x[z] for z in kidx]+[y[lidx]])
    lpdf = gaussian_kde(safevec(l))
    (grid,idx) = makegrid(n)         # 3D grid of coordinates
    lprobs  = lpdf(grid)             # these are density estimates
    lprobs  = lprobs/sum(lprobs)     # normalize to probabiliies
    marprobs = marginalize(n, idx, lprobs)  # marginalized prob
    te = computeTE(n, marprobs)             # computed TE.
    return(te)


def sumLagTE(y,x,yl,n, fromProcess):
    res0 = 0.0
    for yli in range(1,(yl+1)):
        res0 += autoTE(y,x,yli,n)
    return(res0)


def autoshuff((y,x,yl,n)):
    permutedY = deepcopy(y)
    shuffle(permutedY)
    return(sumLagTE(permutedY,x,yl,n,"perm"))


def teShuff((y,x,yl,n)):
    permutedY = deepcopy(y)
    shuffle(permutedY)
    return(autoTE(permutedY,x,yl,n))


def mad(arr):
    """ Median Absolute Deviation: a "Robust" version of standard deviation.
    Indices variabililty of the sample.
    https://en.wikipedia.org/wiki/Median_absolute_deviation
    """
    arr = np.ma.array(arr).compressed() # should be faster to not use masked arrays.
    med = np.median(arr)
    return 1.483*np.median(np.abs(arr - med))


def autoPerm(y,x,yl,n,p,cpus):
    # autoTE is
    # FOR TE (Y -> X)
    # The NULL is that the observed TE is the same as permuted TEs
    #pool = mp.Pool(cpus)
    observedTE = sumLagTE(y,x,yl,n,"obs")
    #permutedList = it.repeat( (y,x,yl,n), p)
    #permutedTE = pool.map(autoshuff, permutedList)
    #robustDist = (observedTE - np.median(permutedTE)) / mad(permutedTE)
    #pool.close()
    #return([robustDist, observedTE] + permutedTE)
    #return([observedTE, robustDist])
    return([observedTE])


def tePerm(y,x,yl,n,p,cpus):
    # autoTE is NO SUMMATIONS
    # FOR TE (Y -> X)
    # The NULL is that the observed TE is the same as permuted TEs
    #pool = mp.Pool(cpus)
    observedTE = autoTE(y,x,yl,n)
    #permutedList = it.repeat( (y,x,yl,n), p)
    #permutedTE = pool.map(teShuff, permutedList)
    #robustDist = (observedTE - np.median(permutedTE)) / mad(permutedTE)
    #pool.close()
    #return([observedTE, robustDist])
    return([observedTE])


def randomWalk(n):
    # do a random walk
    x = [0]
    for i in range(n):
        u = np.random.normal()
        if u > 0:
            xt = x[i]+1
        else:
            xt = x[i]-1
        x.append(xt+np.random.normal(loc=0,scale=0.25))
    return(x)


def autoReg(n):
    # auto regression
    x = [np.random.normal()]
    for i in range(n):
        x.append(0.5*x[i] + np.random.normal())
    return(x)


def model_1_step_back(y):
    # auto regression with coupling
    # sample c -- coupling coefficient
    x = [0 for i in range(0,len(y))]
    for i in range(1,len(y)):
        x[i] = 0.6*x[i-1] + 0.4*y[i-1] + np.random.normal(loc=0, scale=0.25)
    return(x)


def model_3_steps_back(y):
    # auto regression with coupling
    # sample c -- coupling coefficient
    x = [0 for i in range(0,len(y))]
    for i in range(5,len(y)):
        x[i] = 0.6*x[i-1] + 0.4*y[i-3] + np.random.normal(loc=0, scale=0.25)
    return(x)


def model_1_and_3_steps_back(y):
    # auto regression with coupling
    # sample c -- coupling coefficient
    x = [0 for i in range(0,len(y))]
    for i in range(5,len(y)):
        x[i] = 0.4*x[i-1] + 0.4*y[i-1] + 0.2*y[i-3] + np.random.normal(loc=0, scale=0.25)
    return(x)


def main(argv):
    print("type\tforward\treverse")
    n  = 10
    p  = 30

    runs = 1000

    # The random unrelated case #
    yl = 5
    for i in xrange(runs):
        y = np.random.normal(loc=10, scale=1, size=100)
        x = np.random.normal(loc=10, scale=1, size=100)
        res0 = autoPerm(y,x,yl,n,p,6)
        res1 = autoPerm(x,y,yl,n,p,6)
        print("\t".join(["Normal_Random", str(res0[0]), str(res1[0])]))

    # random walk where info comes from 1 and 3 steps back
    yl = 5
    for i in xrange(runs):
        y = randomWalk(100)
        x = model_1_and_3_steps_back(y)
        res0 = autoPerm(y,x,yl,n,p,6)
        res1 = autoPerm(x,y,yl,n,p,6)
        print("\t".join(["RW_Model1_3", str(res0[0]), str(res1[0])]))

    # random walk where info comes from 1 step back
    yl = 5
    for i in xrange(runs):
        y = randomWalk(100)
        x = model_1_step_back(y)
        res0 = autoPerm(y,x,yl,n,p,6)
        res1 = autoPerm(x,y,yl,n,p,6)
        print("\t".join(["RW_Model1", str(res0[0]), str(res1[0])]))

    # random walk where info comes from 3 step back
    yl = 5
    for i in xrange(runs):
        y = randomWalk(100)
        x = model_3_steps_back(y)
        res0 = autoPerm(y,x,yl,n,p,6)
        res1 = autoPerm(x,y,yl,n,p,6)
        print("\t".join(["RW_Model3", str(res0[0]), str(res1[0])]))

    # random walk where info comes from 3 and 1 step back
    # but we DONT use the sum-lag.  just use three step back... 
    yl = 1
    for i in xrange(runs):
        y = randomWalk(100)
        x = model_1_and_3_steps_back(y)
        res0 = tePerm(y,x,yl,n,p,6)
        res1 = tePerm(x,y,yl,n,p,6)
        print("\t".join(["RW_Model_1_3_y1", str(res0[0]), str(res1[0])]))

    # random walk where info comes from 3 and 1 step back
    # but we DONT use the sum-lag.  just use three step back... 
    yl = 3
    for i in xrange(runs):
        y = randomWalk(100)
        x = model_1_and_3_steps_back(y)
        res0 = tePerm(y,x,yl,n,p,6)
        res1 = tePerm(x,y,yl,n,p,6)
        print("\t".join(["RW_Model_1_3_y3", str(res0[0]), str(res1[0])]))

    # random walk where info comes from 3 step back
    # but we DONT use the sum-lag.  just use one step back... 
    yl = 1
    for i in xrange(runs):
        y = randomWalk(100)
        x = model_3_steps_back(y)
        res0 = tePerm(y,x,yl,n,p,6)
        res1 = tePerm(x,y,yl,n,p,6)
        print("\t".join(["RW_Model3_y1", str(res0[0]), str(res1[0])]))


if __name__ == "__main__":
   main(sys.argv)
