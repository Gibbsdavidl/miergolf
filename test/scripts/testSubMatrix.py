import timeit

setup = '''
import scipy.sparse as sp
import numpy as np
from bisect import bisect
from numpy.random import rand, randint
import submatrix as s

r = [10,20,30]
A = s.randomMatrix()
'''
    
t = timeit.Timer("s.subMatrix(r,r,A)", setup).repeat(3, 10)
print t
#print t.timeit()
