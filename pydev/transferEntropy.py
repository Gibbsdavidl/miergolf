import sys
import numpy as np
from copy import copy, deepcopy
import multiprocessing as mp
from numpy.random import shuffle, random, normal
from math import log, sqrt, exp, pi
import itertools as it
from scipy.stats import gaussian_kde

# In this work, I am computing information theoretic values
# by, first, discretizing expression values into a given
# number of bins. Using those bins, the probability of a given
# interval is computed, and the joint probability over time 
# can also be computed (given two time series).

# Want  P(X_t+1, Y_k1, X_k2) * log (P(X_t+1,Y_k1,X_k2)*P(X_t+1)) / (P(X_t+1, X_k2)*P(X_k2,Y_K1))

# just get the joint, then get the others by marginalization

# parameters:
# yk: the markov order for Y = let it be 1
# xk: the markov order for x = let it be 1
# yl: the time delay for y
# xl: the time delay for x
# b : the number of bins

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
    lprobs  = lpdf(grid)          # these are density estimates
    lprobs  = lprobs/sum(lprobs) # normalize to probabiliies
    marprobs = marginalize(n, idx, lprobs)  # marginalized prob
    te = computeTE(n, marprobs)
    return(te)


def shufflete((y,x,yl,h,n)):
    permutedA  = deepcopy(y)
    shuffle(permutedA)
    return(kernelTE(permutedA,x,yl,h,n))


def permTE(y,x,yl,h,n,p,cpus):
    pool = mp.Pool(cpus)
    observedTE = kernelTE(y,x,yl,h,n)
    permutedList = it.repeat( (y,x,yl,h,n), p)
    permutedTE = pool.map(shufflete, permutedList)
    pvalueCount = sum( [1.0 for pte in permutedTE if pte > observedTE] )
    temean = np.mean(permutedTE)
    testd = np.std(permutedTE)
    zscore = (observedTE - temean)/testd
    pool.close()
    return([observedTE, temean, testd, zscore, pvalueCount, (pvalueCount/p)])


def autoTE (y, x, yl):
    # normalize x, y
    n = 10 # grid size
    x = np.array(x); x = (x-x.mean())/max(1,(x-x.mean()).max())
    y = np.array(y); y = (y-y.mean())/max(1,(y-y.mean()).max())
    l = [] # data list  (list of all)
    for i in range(yl, len(x)):
        kidx = [i,(i-1)]  # index for x_t and x_t-1
        lidx = i-yl       # index for y_t-yl
        l.append([x[z] for z in kidx]+[y[lidx]])
    lpdf = gaussian_kde(safevec(l))
    (grid,idx) = makegrid(n)         # 3D grid of coordinates
    lprobs  = lpdf(grid)          # these are density estimates
    lprobs  = lprobs/sum(lprobs) # normalize to probabiliies
    marprobs = marginalize(n, idx, lprobs)  # marginalized prob
    te = computeTE(n, marprobs)
    return(te)


def autoshuff((y,x,yl)):
    permutedA = deepcopy(x)
    shuffle(permutedA)
    return(autoTE(y,permutedA,yl))


def autoPermTE(y,x,yl,p,cpus):
    n = 10
    pool = mp.Pool(cpus)
    observedTE = autoTE(y,x,yl)
    permutedList = it.repeat( (y,x,yl), p)
    permutedTE = pool.map(autoshuff, permutedList)
    pvalueCount = sum( [1.0 for pte in permutedTE if pte > observedTE] )
    temean = np.mean(permutedTE)
    testd = np.std(permutedTE)
    zscore = (observedTE - temean)/testd
    pool.close()
    return([observedTE, temean, testd, zscore, pvalueCount, (pvalueCount/p)])


def runte(filename, fileout, yl, h, n, reps, cpus):
    dat = open(filename,'r').read().strip().split("\n")
    fout = open(fileout,'w')
    dats = map(lambda x: x.split("\t"), dat)
    marker = 0
    for i in xrange(len(dat)):
        for j in xrange(len(dat)):
            id1 = dats[i][0]
            id2 = dats[j][0]
            a = map(float,dats[i][1:41])
            b = map(float,dats[j][1:41])
            res0 = permTE(a,b,yl,h,n,reps,cpus)
            fout.write( str(i)+"\t"+ str(j) +"\t"+ id1 +"\t"+ id2 + "\t" + "\t".join(map(str,res0))+"\n")
            if marker % 117 == 0:
                print str(marker) +"\t"+ str(i)+"\t"+ str(j) +"\t"+ id1 +"\t"+ id2 + "\t" + "\t".join(map(str,res0))
            marker +=1
    fout.close()


def edgelistTE(exprfile, genefile, edgefile, fileout, yl, start, end, reps, cpus):
    genes  = open(genefile,'r').read().strip().split("\n")
    dat    = open(exprfile,'r').read().strip().split("\n")
    dats   = map(lambda x: x.split("\t"), dat)
    edge   = open(edgefile,'r').read().strip().split("\n")
    edges  = map(lambda x: x.split("\t"), edge)
    edgel1 = [x[1] for x in edges]
    edgel2 = [x[2] for x in edges]
    fout   = open(fileout,'w')
    marker = 0
    flag = 1
    for e in xrange(len(edges)):
        try:
            i = genes.index(edgel1[e]) # from 
            j = genes.index(edgel2[e]) # to
            x = map(float,dats[i][start:end])
            y = map(float,dats[j][start:end])
            res0 = autoPermTE(x,y,yl,reps,cpus)
            res1 = autoPermTE(y,x,yl,reps,cpus)
            fout.write("forward"+"\t"+str(i)+"\t"+ str(j) +"\t"+ edgel1[e] +"\t"+ edgel2[e] + "\t" + "\t".join(map(str,res0))+"\n")
            fout.write("reverse"+"\t"+str(j)+"\t"+ str(i) +"\t"+ edgel2[e] +"\t"+ edgel1[e] + "\t" + "\t".join(map(str,res1))+"\n")
            if marker % 17 == 0:
                print str(marker) +"\t"+ str(i)+"\t"+ str(j) +"\t"+ edgel1[e] +"\t"+ edgel2[e] + "\t" + "\t".join(map(str,res0))
                print str(marker) +"\t"+ str(j)+"\t"+ str(i) +"\t"+ edgel2[e] +"\t"+ edgel1[e] + "\t" + "\t".join(map(str,res1))
        except:
            sys.stderr.write("error at: " + str(marker) + "\t" + str(e) + "\t" + edgel1[e] + "\t" + edgel2[e] + "\n")
        marker +=1
    fout.close()

#edgelistTE("../Max_Influence_Problem/Data/GRN/grn_expr_table.txt","../Max_Influence_Problem/Data/GRN/grn_gene_names.txt","../Max_Influence_Problem/Data/GRN/grn_yeastrac_edges.txt", "grn_weights.txt", 1, 1000, 4)

#edgelistTE("../Max_Influence_Problem/Data/GRN/grn_expr_table.txt","../Max_Influence_Problem/Data/GRN/grn_gene_names.txt","../Max_Influence_Problem/Data/GRN/grn_yeastrac_edges.txt", "grn_weights_lag2_cutat28", 2, 28, 1000, 3)

#edgelistTE("grn_expr_table.txt","grn_gene_names.txt","grn_yeastrac_edges.txt", "grn_weights_lag2_full", 2, 41, 1000, 8)
