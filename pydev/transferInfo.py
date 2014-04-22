
from math import log
import numpy as np
import copy

# In this work, I am computing information theoretic values
# by, first, discretizing expression values into a given
# number of bins. Using those bins, the probability of a given
# interval is computed, and the joint probability over time 
# can also be computed (given two time series).

def pr(x, bins):
    # break x into intervals, and return the prob. of 
    # being in one of those intervals.
    xi = (np.histogram(x, bins))[0][range(0,bins)]
    tablexi = map(float,xi)/sum(xi)
    return(tablexi)


def jointPr (x, y, bins):
    # since x and y are time series, 
    # what is prob. of x==i and y==j
    # at time t.
    if len(x) != len(y):
        print "x and y are not equal lengths"
        return()
    xidx = np.histogram(x,bins-1)[1][range(0,bins)]
    xi   = np.digitize(x, xidx)-1
    yidx = np.histogram(y,bins-1)[1][range(0,bins)]
    yi   = np.digitize(y, yidx)-1    
    p = np.zeros(shape=(bins,bins))
    for k in xrange(len(xi)):
        p[(xi[k]), (yi[k])] += 1
    return(p/p.sum())

    
def safelog2(x):
    if (x > 0):
        return(log(x,2))
    else:
        return(0.0)

def h(x, bins):
    # entropy using the above definitions
    px = pr(x,bins)
    logpx = map(safelog2, px)
    hx = -sum(px * logpx)
    return(hx)

def mi(x,y,bins):
    px = pr(x,bins)
    py = pr(y,bins)
    pxy = jointPr(x,y,bins)
    n = pxy.shape[0]
    mixy = 0
    for i in xrange(n):
        for j in xrange(n):
            if (pxy[i][j] > 0 and 
                px[i] > 0 and
                py[j] > 0):
                mixy += pxy[i][j]*log(pxy[i][j]/(px[i]*py[j]), 2)
    return(mixy)


def lagJointPr (x, y, bins, lag):
    # y is being lagged by lag
    # +lag means that y is after x
    # -lag means that y is before x 
    xidx = np.histogram(x,bins-1)[1][range(0,bins)]
    xi   = np.digitize(x, xidx)-1
    yidx = np.histogram(y,bins-1)[1][range(0,bins)]
    yi   = np.tile(np.digitize(y, yidx)-1, 2) # replicate
    p    = np.zeros(shape=(bins,bins))
    for k in xrange(len(xi)):
        p[(xi[k]), (yi[k+lag])] += 1
    return(p/sum(sum(p)))


def lagmi(x,y,bins,lag):
    # the lagged MI #
    x = np.array(x)
    y = np.array(y)
    px = pr(x,bins)
    py = pr(y,bins)
    pxy = lagJointPr(x,y,bins,lag)
    n = pxy.shape[0]
    mixy = 0
    for i in xrange(n):
        for j in xrange(n):
            if (pxy[i][j] > 0 and 
                px[i] > 0 and
                py[j] > 0):
                mixy += pxy[i][j]*log(pxy[i][j]/(px[i]*py[j]), 2)
    return(mixy)



def genInd(x,xorg,k):
    for b in xrange(k):
        y = []
        for xi in x:  # [0]
            for xorgi in xorg:  # 1
                y.append([xorgi] + xi)
        x = copy.copy(y)
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


def twoSeriesDiscreteJointProb (xi, yi, k, bins):
    # since x and y are time series, 
    # what is prob. of x==i and y==j and x(t-1) == m ..etc..
    # let k be the markov order
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

    
def twoSeriesMarkovConditionalDiscrete (xi, yi, k, bins):
    # markov transition matrix for order k
    # p(x(t) | x(t-1), x(t-2), x(t-3), ..., x(t-k), y(t-1), y(t-2) .. y(t-k)) #

    # the probability of the next state, given 
    # the previous k number of states.
    # the transition matrix is indexed in order of time,
    # from t to t-k in that order
 
   # for a markov matrix of order 2
    p = np.zeros(np.repeat(bins, 2*k+1))

    # need t, t-1, t-2, t-k for X
    # need    t-1, t-2, t-k for Y
    
    # generate the joint counts over time points (order k) #
    for i in range(k,len(xi)):
        kidx = range(i, i-k-1, -1) # i and the last K entries in time
        jidx = range(i-1, i-k-1,   -1)
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

a = [1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1]
b = [1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1]


from random import uniform
py = [[0.25,0.75],[0.1,0.9]]
px = [ [[0.5,0.5],[0.6,0.4]], [[0.8,0.2],[0.2,0.8]] ]

y = [0]
for i in xrange(100):
    if uniform(0,1) < py[y[i]][0]:
        y.append(0)
    else:
        y.append(1)


x = [0]
for i in xrange(100):
    if uniform(0,1) < px[y[i]][x[i]][0]:
        x.append(0)
    else:
        x.append(1)

