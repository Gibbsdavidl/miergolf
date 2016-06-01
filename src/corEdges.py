import sys
import numpy as np
from copy import copy, deepcopy
import multiprocessing as mp
from numpy.random import shuffle, random, normal
from math import log, sqrt, exp, pi
import itertools as it
from scipy.stats import gaussian_kde, pearsonr
from scipy.stats import ttest_1samp
from itertools import product

try:
    from Crypto.pct_warnings import PowmInsecureWarning
    import warnings
    warnings.simplefilter("ignore", PowmInsecureWarning)
except:
    pass


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


def autoshuff((x,y)):
    permutedY = deepcopy(y)
    shuffle(permutedY)
    return(pearsonr(x, permutedY)[0])


def autoCorr(x,y,reps1, cpus):
    pool = mp.Pool(cpus)
    observed = pearsonr(x,y)[0]
    permutedList = it.repeat( (x,y), reps1)
    permutedCor = pool.map(autoshuff, permutedList)
    pool.close()
    return([observed] + permutedCor)


def geneindex(gene, genes):
    for i in range(0,len(genes)):
        if gene in genes[i]:
            return(i)
    return(-1)


def prepGeneDataGG(dats, genes, g1, g2):
    i = geneindex(g1, genes) # from
    j = geneindex(g2, genes) # to
    if (i > -1 and j > -1):
        x = map(float,dats[i]) #from
        y = map(float,dats[j]) # to
        x = np.array(x); x = (x-x.mean())/max(1,(x-x.mean()).max())
        y = np.array(y); y = (y-y.mean())/max(1,(y-y.mean()).max())
        return((x,y))
    else:
        return( ([],[]) )


def corEdges(exprfile, genefile, fileout, reps, cpus, g1, g2):
    genes  = open(genefile,'r').read().strip().split("\n")
    dat    = open(exprfile,'r').read().strip().split("\n")
    dats   = map(lambda x: x.split("\t"), dat)
    fout   = open(fileout,'w')
    (fromx,toy) = prepGeneDataGG(dats, genes, g1, g2)
    res0 = autoCorr(fromx,toy,reps, cpus)
    fout.write(g1 +"\t"+ g2 +"\t"+ "\t".join(map(str,res0)) +"\n")
    fout.close()


def maxLagCorEdges(exprfile, genefile, fileout, reps, cpus, ylmax, g1, g2):
    genes  = open(genefile,'r').read().strip().split("\n")
    dat    = open(exprfile,'r').read().strip().split("\n")
    dats   = map(lambda x: x.split("\t"), dat)
    fout   = open(fileout,'w')
    (fromx,toy) = prepGeneDataGG(dats, genes, g1, g2)
    maxCorr = 0.0
    maxLag = 0.0
    for yl in range(0,(ylmax+1)):
        try:
            res0 = autoCorr(fromx,toy,reps, cpus)
            if (res0[0] > maxCorr):
                maxTE  = res0
                maxLag = yl
        except:
            e = sys.exc_info()
            sys.stderr.write(str(e)+"\n")
    fout.write(g1 +"\t"+ g2 +"\t"+ str(maxLag) +"\t"+ str(maxCorr) +"\t"+ "\t".join(map(str,res0)) +"\n")
    fout.close()


def main(argv):
    #for i in range(1,len(argv)):
    #    print(str(i) +"  "+ argv[i])
    exprfile = argv[1]
    genefile = argv[2]
    fileout  =  argv[3]
    reps     = int(argv[4])
    cpus =     int(argv[5])
    g1 =       argv[6]
    g2 =       argv[7]
    maxLagCorEdges(exprfile, genefile, fileout, reps, cpus, 6, g1, g2)


if __name__ == "__main__":
   main(sys.argv)

#pref="/Users/davidlgibbs/Dropbox/Research/Projects/Influence_Maximization_Problem/EserData/"
#pref  = "/users/dgibbs/EserData/"
#genes = pref +"yeast_array_genesymbols.csv"
#gexpr = pref +"Eser_Averaged_Expression.txt"
#tout  = "/Users/davidlgibbs/Desktop/x.txt"
#corEdges(gexpr, genes, tout, 20, 2, "YOX1", "MBP1")
