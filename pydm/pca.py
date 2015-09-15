#!/usr/bin/env python

"""
Created on Sat Jul 21 2015

@author: Chao-Jun Feng

A class  for Principal Component Analysis (PCA)

Usage:
    p = PCA(A ,fraction = 0.90)

In:
    A : an array of e.g. 1000observations x 20 variables, 1000rows x 20 columns
    fraction: use principal components (PCs) that account for e.g. 90% of the total variance

Out:
    p.U, p.S, P.Vt : from scipy.linalg.svd, A =  U.S.Vt
    p.invS: 1/p.S, or 0
    p.eigen: the eigenvalues of A*A, in decreasing order (p.S**2)
    p.numPs:

"""

from __future__ import division, print_function

import numpy as np
#from scipy.linalg import svd
from numpy.linalg import svd

class PCA(object):
    """
    A: an array of e.g.
    metho:
        - SVD  (Sigularity Value Decomposition)
        - EIGD (Eigenvalues Decomposition)
    """
    def __init__(self, A, method = 'SVD'):
        self.A = A
        self.method = method
        if self.method == 'SVD':
            self.U, self.S, self.Vt = svd(A,full_matrices=False)
        elif self.method == 'EIGD':
            pass
        
    def pcs(self,fraction = 0.95 ):
        #assert 0 <= fraction <= 1 ('Fraction should be in range 0.~1.0')
        return self.S**2
