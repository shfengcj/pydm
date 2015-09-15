"""
Created on Sat Jul 19 12:49:34 2015

@author: Chao-Jun Feng

define MDS trainset class

"""
from __future__ import print_function

from .trainset import trainset
import numpy as np

from scipy.linalg import eigh, inv, solve
import matplotlib.pyplot as pl

class mdsTS(trainset):
    """
    MDS analysis of TS

    :param  z: redshift array, m = len(z)
    :param  n: number of models

    """
    def __init__(self, z, n =20):
        super(mdsTS, self).__init__(z,n)

    def Anlysis(self, cosModel, prange, notMix = True, Info = 2):
        # Do the MDS Anlysis
        try:
             if notMix:
                 self.TS = self._genTS(cosModel, prange)
             else:
                 self.TS = self._genTSmix(cosModel, prange)
        except:
            print("Please tell me the generator!")
            
    
        self.disQ = np.zeros((self.n, self.n))
        
        for i in range(self.n):
            for j in range(i, self.n):
                self.disQ[i,j] = ( (self.TS[:,i] - self.TS[:,j])**2 ).sum()
                self.disQ[j,i] = self.disQ[i,j]
                
                
        matrixZ = np.identity(self.n) -  np.ones((self.n, self.n))/self.n
        matrixG = -0.5* np.dot( np.dot(matrixZ, self.disQ), matrixZ)
        
        
        eigens, eigvec = eigh(matrixG)
        
        
        sortIndex = np.argsort(-eigens)
        
        self.eigens = eigens.real[sortIndex]
        self.eigvec = eigvec.real[:,sortIndex]  
        
        
        self.eig1tot = abs(self.eigens).sum()
        self.eig2tot = (self.eigens**2).sum()
        

        self.mdsPercent1 = self.eigens/self.eig1tot
        self.mdsPercent2 = self.eigens**2/self.eig2tot
        
            
        if Info > 2:
        
        
            print('\nMDS: Total absolute eigenvalues: {0:.2f} '.format(self.eig1tot) )
            print('='.join(['='*50]) )
    
            for i in range(5):
                print('The {0} eig: {1:.2f} in {2:.2f} %'
                    .format(i+1, self.eigens[i], 100*self.eigens[i]/self.eig1tot))
    
            print('='.join(['='*50]) )
    
            for i in range(5):
                print("The first {0}s' cumulative eigs: {1:.2f} in {2:.2f} %"
                    .format(i+1, self.eigens[0:i+1].sum(), 100*self.eigens[0:i+1].sum()/self.eig1tot))
            
    
    
            print('\nMDS: Total square eigenvalues: {0:.2f} '.format(self.eig2tot) )
            print('='.join(['='*50]) )
    
            for i in range(5):
                print('The {0} square eig: {1:.2f} in {2:.2f} %'
                    .format(i+1, self.eigens[i]**2, 100*self.eigens[i]**2/self.eig2tot))
    
            print('='.join(['='*50]) )
    
            for i in range(5):
                print("The first {0}s' cumulative eigs: {1:.2f} in {2:.2f} %"
                    .format(i+1, (self.eigens[0:i+1]**2).sum(), 100*(self.eigens[0:i+1]**2).sum()/self.eig2tot))
                    

    def mdsEigs(self):
        return [self.mdsPercent1, self.mdsPercent2]
        
    def mds(self, fraction = 0.95):
        sumvar = np.cumsum(abs(self.eigens))
        sumvar /= sumvar[-1]
        n1 = np.searchsorted(sumvar, fraction) + 1
        
        sumvar = np.cumsum((self.eigens)**2)
        sumvar /= sumvar[-1]
        n2 = np.searchsorted(sumvar, fraction) + 1
                  
        self.nmds = max(n1,n2)
        
        SM = np.dot(self.TS, self.TS.T)        
        #tmp = np.dot(inv(SM), self.TS)
        
        tmp = np.dot(SM, self.TS)

        
        
        
        self.evec = np.zeros((self.m, self.nmds))
        
        
        
        
        
        for i in range(self.nmds):
            rhs = self.eigvec[:,i]*np.sqrt(self.eigens[i])
            self.evec[:,i] = np.dot(tmp,rhs)
            self.evec[:,i] /= self.evec[-1,i]
            
        
        
        
        return self.evec[:,0:self.nmds]
        
        
    
    



    def mdsPlot(self, mds):
        pl.figure()
        for i in mds:
            pl.plot(np.arange(self.m), self.evec[:,i])
        

     
    def mds12Plot(self):
        pl.figure()
        
        x = self.eigvec[:,0]*self.eigens[0]   
        y = self.eigvec[:,1]*self.eigens[1] 
        
        xmid = x.max()-x.min()
        
        #xmid = x.mean()
        
        ymin = -xmid/20.
        ymax =  xmid/20.
               
        
        pl.plot(x, y, 'bo')
        pl.ylim(ymin, ymax)
        
        pl.minorticks_on()
        pl.xlabel('$w_1$', fontsize='x-large')
        pl.ylabel('$w_2$', fontsize='x-large')
        
    
    

























