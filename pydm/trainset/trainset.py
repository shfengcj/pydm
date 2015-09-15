"""
Created on Sat Jul 19 12:49:34 2015

@author: Chao-Jun Feng

define trainset class

"""

import numpy as np

class trainset(object):
    """
    trainset class, to generate train set

    :param  z: redshift array, m = len(z)
    :param  n: number of models

    """
    def __init__(self, z,  n=20):
        self.z  = z
        self.n  = n
        self.m  = len(z)

        self.TS = np.zeros((self.m,self.n))

    def _genTS(self, cosModel, prange):
        # return  TS = n(row) x m(column) matrix
        #:param  cosModel: generator ( from some cosmological model)
        self.cosModel = cosModel
        p = np.zeros(len(prange))

        for i in range(self.n):
            for j, (l, u) in enumerate(prange):
                 p[j] = np.random.uniform(l, u)

            self.cosModel.parasUpdate(p)

            for j in range(self.m):
                self.TS[j,i] = self.cosModel.covDis(self.z[j])
        return self.TS
        
    def _genTSmix(self, cosModel, prange):
        # return  TS = n(row) x m(column) matrix
        #:param  cosModel: generator ( from some cosmological model)
        self.cosModel = cosModel
        
        assert len(self.cosModel) == self.n, ("Check the number of models")        
        
        

        for i in range(self.n):
            p = np.zeros(len(prange[i]))
            
            for j, (l, u) in enumerate(prange[i]):
                 p[j] = np.random.uniform(l, u)

            self.cosModel[i].parasUpdate(p)

            for j in range(self.m):
                self.TS[j,i] = self.cosModel[i].covDis(self.z[j])
        return self.TS
