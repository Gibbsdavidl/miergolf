# David L Gibbs
# Institute for Systems Biology
# April 6 2017

# david.gibbs@systemsbiology.org

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

# In this work, I am computing transfer entropies with variable time lags.

#  P(X_t+1, X_k2, Y_k1) * log (P(X_t+1,Y_k1,X_k2)*P(X_t+1)) / (P(X_t+1, X_k2)*P(X_k2,Y_K1))

# just get the joint, then get the others by marginalization

# parameters:
# yk: the markov order for Y = let it be 1
# xk: the markov order for x = let it be 1
# yl: the time delay for y
# b : the number of bins

# autoTE is
# FOR TE (Y -> X)

#to run it:

# python3 TEDelay.py expression_file.tsv gene_file.tsv output.tsv 10 5 100 4 ADF1 QER2

# where:
# expression_file.tsv has one gene per row, and values for samples in columns, should represent a time series
# gene_file.tsv is a gene name per row, mapping to the rows of the expression_file
# output.tsv is the path to where the TE will be written
# gridsize of 10 seems to work well, can be increased, but does not change things too much in my experience
# number of lags to try... will take the lag that maximizes the robust distance from permutations
# reps1 is the number of permutations,
# g1 gene name found in gene_file
# g2 gene name found in gene_file.


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
    for i in range(n):
        for j in range(n):
            for k in range(n):
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
    for i in range(n):
        for j in range(n):
            for k in range(n):
                l = idx[(i,j,k)]
                px[j] += probs[l]
                pxx[i][j] += probs[l]
                pxy[j][k] += probs[l]
                pxyz[i][j][k] += probs[l]
    return((px,pxx,pxy,pxyz))


def computeTE(n,px,pxx,pxy,pxyz):
    te = 0.0
    for i in range(n):
        for j in range(n):
            for k in range(n):
                num = pxyz[i][j][k]*px[j]
                den = pxx[i][j]*pxy[j][k]
                if den > 0.0:
                    te += pxyz[i][j][k]*safelog(num/den)
    return(te)


def autoTE (y, x, yl, n):
    # x,y are vectors
    # yl is the lag time#\
    # n is the grid size
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
    te = computeTE(n, marprobs[0], marprobs[1], marprobs[2], marprobs[3]) # computed TE.
    return(te)


def autoTEList ( plist):
    l = [] # data list  (list of all)
    y=plist[0]; x=plist[1]; yl=plist[2]; n=plist[3]
    for i in range(yl, len(x)):
        kidx = [i,(i-1)]  # index for x_t and x_t-1
        lidx = i-yl       # index for y_t-yl
        l.append([x[z] for z in kidx]+[y[lidx]])
    lpdf = gaussian_kde(safevec(l))
    (grid,idx) = makegrid(n)         # 3D grid of coordinates
    lprobs  = lpdf(grid)             # these are density estimates
    lprobs  = lprobs/sum(lprobs)     # normalize to probabiliies
    marprobs = marginalize(n, idx, lprobs)  # marginalized prob
    te = computeTE(n, marprobs[0], marprobs[1], marprobs[2], marprobs[3]) # computed TE.
    return(te)


def chunkShuffle( dat, chunkSize ):
    chunks = [dat[x:x+chunkSize] for x in range(0, len(dat), chunkSize)]
    shuffle(chunks)
    return( [item for sublist in chunks for item in sublist] )


def autoshuff( plist ):
    y = plist[0] # the target sequence
    s = plist[4] # the chunk size
    permutedY = chunkShuffle( deepcopy(y), s )
    plist[0] = permutedY
    return(autoTEList(plist))


def autoPerm(y,x,yl,n,p,s,cpus):
    # x,y are vectors
    # yl is the lag time#\
    # n is the grid size
    # p is the number of permutations
    # s is the chunk size for shuffling
    pool = mp.Pool(cpus)
    observedTE = autoTE(y,x,yl,n)
    permutedList = it.repeat( [y,x,yl,n,s], p)
    permutedTE = pool.map(autoshuff, permutedList)
    robustDist = (observedTE - np.median(permutedTE)) / mad(permutedTE)
    pool.close()
    return([robustDist, observedTE] + permutedTE)


def pvals (lagres):
    p = (sum(lagres[2:] > lagres[1])+1) / (len(lagres)-2)
    return([lagres[1], lagres[0], p])


def scanLags (y,x,yl,n,p,s,cpus):
    lagResults = [autoPerm(y,x,i,n,p,s,cpus) for i in range(0,(yl+1))]
    return([pvals(lr) for lr in lagResults])


def maxLag (y,x,yl,n,p,s,cpus):
    scanRes = scanLags(y,x,yl,n,p,s,cpus)
    dists = np.array([x[1] for x in scanRes])
    idx = dists.argmax()
    return(scanRes[idx])


def mad(arr):
    """ Median Absolute Deviation: a "Robust" version of standard deviation.
    Indices variabililty of the sample.
    https://en.wikipedia.org/wiki/Median_absolute_deviation
    """
    arr = np.ma.array(arr).compressed() # should be faster to not use masked arrays.
    med = np.median(arr)
    return 1.483*np.median(np.abs(arr - med))


def geneindex(gene, genes):
    for i in range(0,len(genes)):
        if gene in genes[i]:
            return(i)
    return(-1)


def prepGeneDataGG(dats, genes, g1, g2):
    i = geneindex(g1, genes) # from
    j = geneindex(g2, genes) # to
    if (i > -1 and j > -1):
        x = dats[i].split("\t")
        y = dats[j].split("\t")
        x = [float(a) for a in x]
        y = [float(a) for a in y]
        x = np.array(x)
        y = np.array(y)
        x = (x-x.mean())/max(1,(x-x.mean()).max())
        y = (y-y.mean())/max(1,(y-y.mean()).max())
        return((x,y))
    else:
        return( ([],[]) )


def TEMaxLagOneEdge(exprfile, genefile, fileout, gridsize, ylmax, reps1, cpus, g1, g2):
    genes = open(genefile,'r').read().strip().split("\n")
    dat   = open(exprfile,'r').read().strip().split("\n")
    fout  = open(fileout,'w')
    try:
        (fromy,tox) = prepGeneDataGG(dat, genes, g1, g2)
        res0 = maxLag(fromy,tox,ylmax,gridsize,reps1,3,cpus)
        fout.write('\t'.join([g1, g2] + [str(x) for x in res0]) + "\n")
    except:
        fout.write('error\n')
        e = sys.exc_info()
        sys.stderr.write(str(e)+"\n")


def TEMaxEdgeList(exprfile, genefile, edgefile, fileout, gridsize, ylmax, reps1, cpus):
    fout  = open(fileout,'w')
    edges = open(edgefile,'r').read().strip().split("\n")
    for ei in edges:
        gs = ei.split('\t')
        print(gs)
        TEMaxLagOneEdge(exprfile, genefile, fout, gridsize, ylmax, reps1, gs[0], gs[1], cpus)
    fout.close()


def main(argv):
    #for i in range(1,len(argv)):
    #    print(str(i) +"  "+ argv[i])
    exprfile = argv[1]
    genefile = argv[2]
    edgefile = argv[3]
    fileout  = argv[4]
    gridsize = int(argv[5])
    ylmax    = int(argv[6])
    reps1    = int(argv[7])
    cpus =     int(argv[8])
    g1 = argv[9]
    g2 = argv[10]
    TEMaxLagOneEdge(exprfile, genefile, fileout, gridsize, ylmax, reps1, cpus, g1, g2)


if __name__ == "__main__":
   main(sys.argv)


#to run it:

# python3 TEDelay.py expression_file.tsv gene_file.tsv output.tsv 10 5 100 4 ADF1 QER2

# where:
# expression_file.tsv has one gene per row, and values for samples in columns, should represent a time series
# gene_file.tsv is a gene name per row, mapping to the rows of the expression_file
# output.tsv is the path to where the TE will be written
# gridsize of 10 seems to work well, can be increased, but does not change things too much in my experience
# number of lags to try... will take the lag that maximizes the robust distance from permutations
# reps1 is the number of permutations,
# g1 gene name found in gene_file
# g2 gene name found in gene_file.
