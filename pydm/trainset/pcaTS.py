"""
Created on Sat Jul 19 12:51:34 2015

@author: Chao-Jun Feng

define PCA trainset class

"""
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as pl

from .trainset import trainset


from scipy.linalg import eigh


class pcaTS(trainset):
    """
    PCA analysis of TS

    :param  z: redshift array, m = len(z)
    :param  n: number of models

    """
    def __init__(self, z, n =20):
        super(pcaTS, self).__init__(z,n)

    def Anlysis(self, cosModel, prange, notMix = True, Info = 2):
        
        # Do the PCA Anlysis
        try:
            if notMix:
                self.TS = self._genTS(cosModel, prange)
            else:
                self.TS = self._genTSmix(cosModel, prange)
        except:
            print("Please tell me the generator!")

        for i in range(self.n):
            self.tsMean = self.TS.mean(axis=1)
            self.TS[:,i] -= self.tsMean

        SM = np.dot(self.TS, self.TS.T)/self.n

        eigens, eigvec = eigh(SM)
        
        
        sortIndex = np.argsort(-eigens)
        
        self.eigens = eigens[sortIndex]
        self.eigvec = eigvec[:,sortIndex]       
        
        self.eigtot = self.eigens.sum()

        if self.eigvec[1,0] < 0:
            self.eigvec = -self.eigvec
            

        self.pcsPercent = self.eigens/self.eigtot
        
        if Info > 2:

            print('\nPCA: Total eigenvalues: {0:.2f} '.format(self.eigtot) )
            print('='.join(['='*50]) )
    
            for i in range(5):
                print('The {0} eig: {1:.2f} in {2:.2f} %'
                    .format(i+1, self.eigens[i], 100*self.eigens[i]/self.eigtot))
    
            print('='.join(['='*50]) )
    
            for i in range(5):
                print("The first {0}s' cumulative eigs: {1:.2f} in {2:.2f} %"
                    .format(i+1, self.eigens[0:i+1].sum(), 100*self.eigens[0:i+1].sum()/self.eigtot))
            
        

    def pcEigs(self):
        return self.pcsPercent
    
    def pcs(self, fraction = 0.95):
        sumvar = np.cumsum(self.eigens)
        sumvar /= sumvar[-1]
        self.npcs = np.searchsorted(sumvar, fraction) + 1

        return self.eigvec[:,0:self.npcs]

    def pcsPlot(self, pcs):
        pl.figure()
        pl.plot(self.z, self.eigvec[:,pcs])
        
    def pc12Plot(self):
        pl.figure()
        
        x = np.dot(self.eigvec[:,0], self.TS)    
        y = np.dot(self.eigvec[:,1], self.TS)
        
        #xmid = x.max()-x.min()
        
        xmid = x.mean()
        
        ymin = -xmid/2.
        ymax =  xmid/2.
               
        
        pl.plot(x, y, 'bo')
        pl.ylim(ymin, ymax)
        
        pl.minorticks_on()
        pl.xlabel('$w_1$', fontsize='x-large')
        pl.ylabel('$w_2$', fontsize='x-large')
        
