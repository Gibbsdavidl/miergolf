
import numpy as np
from copy import copy, deepcopy
import multiprocessing as mp
from random import uniform, shuffle
from math import log

# In this work, I am computing information theoretic values
# by, first, discretizing expression values into a given
# number of bins. Using those bins, the probability of a given
# interval is computed, and the joint probability over time 
# can also be computed (given two time series).

    
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
    

def markovConditionalProb(x, k, bins):
    # markov transition matrix for order k
    # p(x(t) | x(t-1), x(t-2), x(t-3), ..., x(t-k)) #
    # the probability of the next state, given 
    # the previous k number of states.
    # the transition matrix is indexed in order of time,
    # from t to t-k in that order
 
    # for a markov matrix of order 2
    # p[x[t], x[t-1], x[t-2]] is how one would index.
    xidx = np.histogram(x,bins)[1][range(0,bins)] # just need the left edges of bins
    xi   = np.digitize(x, xidx)-1   # for each x_i, what bin? zero indexed
    p    = np.zeros(np.repeat(bins, k+1)) # K dimensional array of B number of bins
    # generate the joint counts over time points (order k) #
    for i in range(k,len(xi)):
        kidx = range(i, i-k-1, -1) # i and the last K entries in time
        p[tuple([xi[z] for z in kidx])] += 1 # count that state sequence

    # to compute the conditional probability
    # counts of t `div` sum over states (t-1 t-2 .. t-k)
    # first to compute the "column-sums"
    psum = sum(map(lambda i: p[i], xrange(bins)))
        
    # then for each "column" of states (say t-1, t-2) sum it
    for tup in generateIndices(k, bins):
        pcol = psum[tup[1:len(tup)]] # the tail of tup
        if pcol > 0:
            p[tup] = p[tup] / pcol
        else:
            p[tup] = 0.0
    return(p)


def twoSeriesJointProb (x, y, k, bins):
    # since x and y are time series, 
    # what is prob. of x==i and y==j and x(t-1) == m ..etc..
    # let k be the markov order
    xidx = np.histogram(x,bins)[1][range(0,bins)] # just need the left edges of bins
    xi   = np.digitize(x, xidx)-1  # for each x_i, what bin? zero indexed
    yidx = np.histogram(y,bins)[1][range(0,bins)] # just need the left edges of bins
    yi   = np.digitize(y, yidx)-1  # for each y_i, what bin? zero indeyed
    p    = np.zeros(np.repeat(bins, 2*k+1))
    # generate the joint counts over time points (order k) #
    for i in range(k,len(xi)):
        kidx = range(i,   i-k-1, -1) # i and the last K entries in time\        
        lidx = range(i-1, i-k-1, -1) # i and the last K entries in time
        tx = [xi[z] for z in kidx]
        ty = [yi[z] for z in lidx]
        p[tuple(tx+ty)] += 1 # count that state sequence
    return(p/(p.sum()))

    
    
def twoSeriesMarkovConditional (x, y, k, bins):
    # markov transition matrix for order k
    # p(x(t) | x(t-1), x(t-2), x(t-3), ..., x(t-k), y(t-1), y(t-2) .. y(t-k)) #

    # the probability of the next state, given 
    # the previous k number of states.
    # the transition matrix is indexed in order of time,
    # from t to t-k in that order
 
   # for a markov matrix of order 2
    xidx = np.histogram(x,bins)[1][range(0,bins)] # just need the left edges of bins
    xi   = np.digitize(x, xidx)-1  # for each x_i, what bin? zero indexed
    yidx = np.histogram(y,bins)[1][range(0,bins)] # just need the left edges of bins
    yi   = np.digitize(y, yidx)-1  # for each y_i, what bin? zero indeyed
    p    = np.zeros(np.repeat(bins, 2*k+1))

    # need t, t-1, t-2, t-k for X
    # need    t-1, t-2, t-k for Y
    
    # generate the joint counts over time points (order k) #
    for i in range(k,len(xi)):
        kidx = range(i,   i-k-1, -1) # i and the last K entries in time
        jidx = range(i-1, i-k-1, -1)
        tx = [xi[z] for z in kidx]
        ty = [yi[z] for z in jidx] 
        p[tuple(tx+ty)] += 1 # count that state sequence

    # to compute the conditional probability
    # counts of t `div` sum over states (t-1 t-2 .. t-k)
    # first to compute the "column-sums"
    psum = sum(map(lambda i: p[i], xrange(bins)))
        
    # then for each "column" of states (say t-1, t-2) sum it
    for tup in generateIndices(2*k, bins):
        pcol = psum[tup[1:len(tup)]] # the tail of tup
        if pcol > 0:
            p[tup] = p[tup] / pcol
        else:
            p[tup] = 0.0
    return(p)

   

def transferEntropy(a,b,k,bins):
    jointAB = twoSeriesJointProb(a,b,k,bins)
    conditionalAB = twoSeriesMarkovConditional(a,b,k,bins)
    conditionalA = markovConditionalProb(a,k,bins)
    bigSum = 0    
    # then for each "column" of states (say t-1, t-2) sum it
    for tup in generateIndices(2*k, bins):
        tupi = tuple([tup[i] for i in range(0,k+1)])
        if conditionalAB[tup] > 0 and conditionalA[tupi] > 0:
            bigSum += jointAB[tup] * log(conditionalAB[tup]/conditionalA[tupi],10)
    return(bigSum)


def te((a,b,k,bins)):
    return(transferEntropy(a,b,k,bins))

def transferEntropyPermutationA(a,b,k,bins,p,cpus):
    pool = mp.Pool(cpus)
    observedTE = transferEntropy(a,b,k,bins)
    permutedA  = deepcopy(a)
    permutedList = []
    for i in xrange(p):
        shuffle(permutedA)
        aprime = deepcopy(permutedA)
        permutedList.append((aprime,b,k,bins))
    permutedTE = pool.map(te, permutedList)
    pvalueCount = sum( [1.0 for pte in permutedTE if pte > observedTE] )
    temean = np.mean(permutedTE)
    testd = np.std(permutedTE)
    zscore = (observedTE - temean)/testd
    pool.close()
    return([observedTE, temean, testd, zscore, pvalueCount, (pvalueCount/p)])


def transferEntropyUniformBackground(a,b,k,bins,p,cpus):
    pool = mp.Pool(cpus)
    observedTE = transferEntropy(a,b,k,bins)
    permutedList = []
    for i in xrange(p):
        permutedN  = [uniform(min(a), max(a)) for i in xrange(len(a))]
        nprime = deepcopy(permutedN)
        permutedList.append((a,nprime,k,bins))
    permutedTE = pool.map(te, permutedList)
    pvalueCount = sum( [1.0 for pte in permutedTE if pte < observedTE] )
    temean = np.mean(permutedTE)
    testd = np.std(permutedTE)
    zscore = (observedTE - temean)/testd
    pool.close()
    return([observedTE, temean, testd, zscore, pvalueCount, (pvalueCount/p)])


def telist(filename, fileout):
    flag = 10
    dat = open(filename,'r').read().strip().split("\n")
    fout = open(fileout,'w')
    dats = map(lambda x: x.split("\t"), dat)
    for i in xrange(len(dat)):
        for j in xrange(len(dat)):
            id1 = dats[i][0]
            id2 = dats[j][0]
            a = map(float,dats[i][1:])
            b = map(float,dats[j][1:])
            res0 = transferEntropyPermutationA(a,b,2,3,200,4)
            print str(i)+"\t"+ str(j) +"\t"+ id1 +"\t"+ id2 + "\t" + "\t".join(map(str,res0))
            fout.write( str(i)+"\t"+ str(j) +"\t"+ id1 +"\t"+ id2 + "\t" + "\t".join(map(str,res0))+"\n")
    fout.close()


def teregionlist(filename, fileout, start, end):
    flag = 10
    dat = open(filename,'r').read().strip().split("\n")
    fout = open(fileout,'w')
    dats = map(lambda x: x.split("\t"), dat)
    for i in range(start,end):
        for j in range(start,end):
            id1 = dats[i][0]
            id2 = dats[j][0]
            a = map(float,dats[i][1:25])
            b = map(float,dats[j][1:25])
            res0 = transferEntropyPermutationA(a,b,2,3,200,4)
            print str(i)+"\t"+ str(j) +"\t"+ id1 +"\t"+ id2 + "\t" + "\t".join(map(str,res0))
            fout.write( str(i)+"\t"+ str(j) +"\t"+ id1 +"\t"+ id2 + "\t" + "\t".join(map(str,res0))+"\n")
    fout.close()


    #In [57]: transferEntropy(x18,x2,1,3)
    #Out[57]: 0.084551430762819499

    #In [58]: transferEntropy(x18,x7,1,3)
    #Out[58]: 0.11275223170972035

    #In [59]: transferEntropy(x18,r,1,3)
    #Out[59]: 0.12013085312581222
