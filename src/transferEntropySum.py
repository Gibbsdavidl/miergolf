import sys
import numpy as np
from copy import copy, deepcopy
import multiprocessing as mp
from numpy.random import shuffle, random, normal
from math import log, sqrt, exp, pi
import itertools as it
from scipy.stats import gaussian_kde
from itertools import product


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
    lprobs  = lpdf(grid)          # these are density estimates
    lprobs  = lprobs/sum(lprobs) # normalize to probabiliies
    marprobs = marginalize(n, idx, lprobs)  # marginalized prob
    te = computeTE(n, marprobs)
    return(te)


def sumLagTE(y,x,yl,n):
    res0 = 0.0
    for yl in range(1,(yl+1)):
        res0 += autoTE(y,x,yl,n)
    return(res0)


def autoshuff((y,x,yl,n)):
    permutedA = deepcopy(x)
    shuffle(permutedA)
    return(sumLagTE(y,permutedA,yl,n))


def autoPerm(y,x,yl,n,p,cpus):
    # autoTE is
    # FOR TE (Y -> X)
    pool = mp.Pool(cpus)
    observedTE = sumLagTE(y,x,yl,n)
    permutedList = it.repeat( (y,x,yl,n), p)
    permutedTE = pool.map(autoshuff, permutedList)
    pool.close()
    return([observedTE]+permutedTE)


def getEdges(edgefile):
    edge   = open(edgefile,'r').read().strip().split("\n")
    edges  = map(lambda x: x.split("\t"), edge)
    edge1 = [x[0] for x in edges]
    edge2 = [x[1] for x in edges]
    return( (edge1, edge2) )


def prepGeneData(dats, genes, edgeFrom, edgeTo, e):
    i = genes.index(edgeFrom[e]) # from
    j = genes.index(edgeTo[e]) # to
    x = map(float,dats[i]) #from
    y = map(float,dats[j]) # to
    x = np.array(x); x = (x-x.mean())/max(1,(x-x.mean()).max())
    y = np.array(y); y = (y-y.mean())/max(1,(y-y.mean()).max())
    return((x,y))



def summedTELag(exprfile, genefile, edgefile1, fileout, gridsize, ylmax, reps1, reps2, thresh1, cpus):
    genes  = open(genefile,'r').read().strip().split("\n")
    dat    = open(exprfile,'r').read().strip().split("\n")
    dats   = map(lambda x: x.split("\t"), dat)
    fout   = open(fileout,'w')
    (edgeFrom, edgeTo) = getEdges(edgefile1)

    #for e in xrange(len(edgeFrom)):
    for e in xrange(10):
        print(e)
        (fromy,tox) = prepGeneData(dats, genes, edgeFrom, edgeTo, e)
        # build list of TEs
        try:
            res0 = autoPerm(fromy,tox,ylmax,gridsize,reps1,cpus)
            # get number of perms greater than observed
            t = sum([a > res0[0] for a in res0[1:]]) / float(reps1)
            if (t < thresh1):
                # if there's fewer than our thresh CONTINUE with more perms
                res1 = autoPerm(fromy,tox,ylmax,gridsize,reps2,cpus)
                t = sum([a > res1[0] for a in res1[1:]]) / float(reps2)
                if (t < thresh1):
                    fout.write(str(edgeFrom[e]) +"\t"+ str(edgeTo[e]) +"\t"+ "\t".join(map(str,res1)) +"\n")
        except:
            sys.stderr.write("error at " + "\t".join([str(e), edgeFrom[e], edgeTo[e]]) + "\n")
    fout.close()


pref = "/Users/davidgibbs/Dropbox/Research/Projects/Influence_Maximization_Problem/Data/"
ef1 = pref+"Edges_w_Extra_TFs_Dec4_2014.txt"
ef2 = pref+"Edges_w_Extra_TFs_Dec4_2014.txt"
genes = pref+"yeast_array_genesymbols.csv"
gexpr = pref+"Eser_Averaged_Expression.txt"
tout = "/Users/davidgibbs/Desktop/test.txt"
summedTELag(gexpr, genes, ef1, tout, 10, 8, 3, 20, 0.05, 2)



#  step 1: small number of perms
#       2: if number of perms is higher
#       3:    do more terms











#forward	1025	1535	ACE2	BMH1	8	0.00509469519606	0.0410603541254	0.022660898014	-1.58712416901	18.0	1.0
#forward	1025	1910	ACE2	BUD9	4	0.400446640752	0.0166262164103	0.00628134427359	61.1048220929	0	0
#forward	1025	5183	ACE2	CDC6	2	0.301548888897	0.0421773756873	0.0153413053307	16.9067434366	0	0

#edgelistTElagRange("../Max_Influence_Problem/Data/GRN/grn_expr_table.txt","../Max_Influence_Problem/Data/GRN/grn_gene_names.txt","../Max_Influence_Problem/Data/GRN/grn_yeastrac_edges.txt", "grn_weights.txt", 1, 1000, 4)

#edgelistTE("../Max_Influence_Problem/Data/GRN/grn_expr_table.txt","../Max_Influence_Problem/Data/GRN/grn_gene_names.txt","../Max_Influence_Problem/Data/GRN/grn_yeastrac_edges.txt", "grn_weights_lag2_cutat28", 2, 28, 1000, 3)

#edgelistTE("grn_expr_table.txt","grn_gene_names.txt","grn_yeastrac_edges.txt", "grn_weights_lag2_full", 2, 41, 1000, 8)


#edgelistTE("Data/exprmat.csv","Data/genesymbols.csv","Data/Eser_Paper_SI/edgesToTest.tsv", "EserWeights", 2, 2, 41, 2000, 4)

x=[ 10.783598,  9.327544,  8.879752 , 9.127505,  9.159929 , 9.214980 , 9.190469 , 9.438197,  9.229540 , 9.975117, 10.302772, 11.134553, 11.125216, 11.139961,
    10.927145, 10.975110, 10.705192, 10.398936, 10.341666, 10.513413, 10.471136, 10.934537, 11.109571, 11.372016, 11.376848, 11.265910, 11.082069, 11.123436,
    10.576115, 10.717450, 10.342657, 10.196068, 10.167991, 10.627770, 10.810851, 10.707537, 10.623126, 10.678473, 10.480142, 10.374172, 10.264492, 10.872016]

y=[8.177953,  8.337438,  8.460188,  8.207932,  8.832201,  8.363152,  8.398193,  8.247786,  8.818389,  9.271919,  9.769506,  9.912917, 10.042233,  9.674573,
   9.777879,  9.119085,  8.963061,  8.841799,  8.984359,  8.626832,  9.065240,  9.018956,  9.532496,  9.652034,  9.515809,  9.819796,  9.647349,  9.075306,
   9.302837,  9.514659,  9.122192,  9.247045,  8.989483,  9.497320,  9.589326,  9.791079, 10.291129, 10.004832,  9.857048, 10.287642, 10.416445,  9.225077]


#paramSweep(x,y)

#laglistTE("/Volumes/YosemiteToast/Users/davidgibbs/Dropbox/Research/Projects/Influence_Maximization_Problem/Results/Yeast/new/Eser_Averaged_Expression.txt",
#          "/Volumes/YosemiteToast/Users/davidgibbs/Dropbox/Research/Projects/Influence_Maximization_Problem/Results/Yeast/new/genesymbols.csv",
#          "/Volumes/YosemiteToast/Users/davidgibbs/Dropbox/Research/Projects/Influence_Maximization_Problem/Results/Yeast/new/small_edge_list.txt",
#          "/Volumes/YosemiteToast/Users/davidgibbs/Dropbox/Research/Projects/Influence_Maximization_Problem/Results/Yeast/new/small_results.txt",
#          0, 42, 100, 1)



# edgelistTElagRange("Eser_Averaged_Expression.txt","genesymbols.csv","Edges_w_Extra_TFs_Dec4_2014.txt", "New_Lag_Results.txt", 5, 0, 42, 2000, 10)
