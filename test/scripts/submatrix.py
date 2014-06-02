
import sys
import scipy.sparse as sp
import scipy.sparse.linalg as lin
import numpy as np
from bisect import bisect
from numpy.random import rand, randint
import timeit

def randomMatrix():
    A = sp.lil_matrix((10000, 10000))
    A[0, :1000] = rand(1000)
    A[1, 1000:2000] = A[0, :1000]
    A.setdiag(rand(10000))
    return(A)

    
def subMatrix(rows, cols, A):
    return(A.tocsr()[rows,:].tocsc()[:,cols])


