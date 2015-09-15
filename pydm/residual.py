# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 06:52:02 2015

@author: Chao-Jun Feng

To calculate the residuals from data.
One can added new residual (or data) by inheriting the resiCal class

"""

# 1pc = 3.0857e18 cm
PC2CM = 3.0857e18

import numpy as np

class ResiCal:

    def _resGen(self, cov, res):
        # To do the cholesky decomposition of the Covariance matrix
        lowtri = np.linalg.cholesky(cov)
        invlow = np.linalg.inv(lowtri)
        residual = np.dot(invlow, res.T)

        return residual

    def residual(p, fjac=None):
        res = 0
        return res

# residual for JLA supernova
class JLAresiCal(ResiCal):

    def __init__(self, cosModel, DATA_DIR_JLA):
        """
        cosModel: a object of class CosModel
        DATA_DIR_JLA: the directory of JLA data, which should look like this:

        DATA_DIR_JLA
            - JLA.paramnames
            - jla_lcparams.txt
            - jla.dataset
            - jla_mub.txt
            - jla_mub_covmatrix.dat
            - jla_v0b_covmatrix.dat
            - jla_simple.dataset
            - jla_va_covmatrix.dat
            - jla_v0_covmatrix.dat
            - jla_vab_covmatrix.dat
            - jla_v0a_covmatrix.dat
            - jla_vb_covmatrix.dat

        """
        self.dataDir = DATA_DIR_JLA
        self.cosModel = cosModel


        dataFile = open(DATA_DIR_JLA+'/jla_lcparams.txt', 'r')

        self.name      = np.array([])
        self.zcmb      = np.array([])
        self.zhel      = np.array([])
        self.dz        = np.array([])
        self.mb        = np.array([])
        self.dmb       = np.array([])
        self.x1        = np.array([])
        self.dx1       = np.array([])
        self.color     = np.array([])
        self.dcolor    = np.array([])
        self.h3rdvar   = np.array([])
        self.dh3rdvar  = np.array([])
        self.tmax      = np.array([])
        self.dtmax     = np.array([])
        self.cov_m_s   = np.array([])
        self.cov_m_c   = np.array([])
        self.cov_s_c   = np.array([])
        self.set       = np.array([])
        self.ra        = np.array([])
        self.dec       = np.array([])
        self.biascor   = np.array([])

        for line in dataFile:
            line = line.strip()
            if line.startswith('#'):
                temp=line.strip('#').split()
                self.snhead = temp
            else:
                temp=line.strip('\r\n').split()

                self.name      = np.append(self.name,       temp[0])
                self.zcmb      = np.append(self.zcmb,       float(temp[1]) )
                self.zhel      = np.append(self.zhel,       float(temp[2]) )
                self.dz        = np.append(self.dz,         float(temp[3]) )
                self.mb        = np.append(self.mb,         float(temp[4]) )
                self.dmb       = np.append(self.dmb,        float(temp[5]) )
                self.x1        = np.append(self.x1,         float(temp[6]) )
                self.dx1       = np.append(self.dx1,        float(temp[7]) )
                self.color     = np.append(self.color,      float(temp[8]) )
                self.dcolor    = np.append(self.dcolor,     float(temp[9]) )
                self.h3rdvar   = np.append(self.h3rdvar,    float(temp[10]))
                self.dh3rdvar  = np.append(self.dh3rdvar,   float(temp[11]))
                self.tmax      = np.append(self.tmax,       float(temp[12]))
                self.dtmax     = np.append(self.dtmax,      float(temp[13]))
                self.cov_m_s   = np.append(self.cov_m_s,    float(temp[14]))
                self.cov_m_c   = np.append(self.cov_m_c,    float(temp[15]))
                self.cov_s_c   = np.append(self.cov_s_c,    float(temp[16]))
                self.set       = np.append(self.set,        float(temp[17]))
                self.ra        = np.append(self.ra,         float(temp[18]))
                self.dec       = np.append(self.dec,        float(temp[19]))
                self.biascor   = np.append(self.biascor,    float(temp[20]))

        dataFile.close()

        self.numData = len(self.name)

        #Read Covariance Matrix
        self.C00 = self.__covRead('/jla_v0_covmatrix.dat')
        self.C11 = self.__covRead('/jla_va_covmatrix.dat')
        self.C22 = self.__covRead('/jla_vb_covmatrix.dat')
        self.C01 = self.__covRead('/jla_v0a_covmatrix.dat')
        self.C02 = self.__covRead('/jla_v0b_covmatrix.dat')
        self.C12 = self.__covRead('/jla_vab_covmatrix.dat')


    def __covRead(self, file_name):

        tmp = np.fromfile(self.dataDir + file_name, sep=" ")
        n = tmp[0]
        cov = tmp[1:]
        cov = np.reshape(cov, (n,n))
        return cov

    def zmuerr(self,p):
        scriptmcut = 10.

        self.cosModel.parasUpdate(p)



        covMatrix = self.C00 + p[0]*p[0]*self.C11 +p[1]*p[1]*self.C22 \
                  + 2.0*p[0]*self.C01 - 2.0*p[1]*self.C02 - 2.0*p[0]*p[1]*self.C12

        modulus = np.zeros(self.numData)
        modtheory = np.zeros(self.numData)

        for i in range(self.numData):
            if self.h3rdvar[i] > scriptmcut:
                prob = 1.0
            else:
                prob = 0.0

            modulus[i] = self.mb[i] - (p[2]-p[0]*self.x1[i] + p[1]*self.color[i] + p[3]*prob)
            covMatrix[i,i] += self.dmb[i]*self.dmb[i] + (p[0]*self.dx1[i])*(p[0]*self.dx1[i]) + \
                          (p[1]*self.dcolor[i])*(p[1]*self.dcolor[i]) + 2.0*p[0]*self.cov_m_s[i] - \
                           2.0*p[1]*self.cov_m_c[i] - 2.0*p[0]*p[1]*self.cov_s_c[i]

            cdz = self.cosModel.covDis(self.zcmb[i])
            modtheory[i] = 5.0*np.log10( (1.0+self.zhel[i])* cdz *3.e5/p[4] )+25.0

        self.modulus   = modulus   # record the modulus
        self.covMatrix = covMatrix # record the covMatrix
        zmuerr = {}
        zmuerr['z'] = self.zcmb
        zmuerr['mu'] = self.modulus
        zmuerr['err'] = np.diagonal(self.covMatrix).copy()
        zmuerr['mut'] = modtheory

        return zmuerr


    def residual(self, p, fjac=None):

        """
        Parameter informations:
        # JLA nuiance parameters
        # ==========================
        # p[0]:alpha
        # p[1]:beta
        # p[2]:M
        # p[3]:DeltaM
        #
        # Cosmological parameters
        # ==========================
        # p[4]:H0
        # p[5]:Omega_bh2
        # p[6]:Omega_m
        # p[7]:Omega_rh2
        # p[8]:Omega_k
        # p[9]... other Cosmological parameters
        """

        scriptmcut = 10.

        covMatrix = self.C00 + p[0]*p[0]*self.C11 +p[1]*p[1]*self.C22 \
                  + 2.0*p[0]*self.C01 - 2.0*p[1]*self.C02 - 2.0*p[0]*p[1]*self.C12

        # remember to update the model parameter before calculation
        self.cosModel.parasUpdate(p)

        res = np.zeros(self.numData)


        for i in range(self.numData):
            cdz = self.cosModel.covDis(self.zcmb[i])

            if self.h3rdvar[i] > scriptmcut:
                prob = 1.0
            else:
                prob = 0.0

            res[i]  = self.mb[i] - (p[2]-p[0]*self.x1[i] + p[1]*self.color[i] + p[3]*prob) \
                        - ( 5.0*np.log10( (1.0+self.zhel[i])* cdz *3.e5/p[4] )+25.0)

            covMatrix[i,i] += self.dmb[i]*self.dmb[i] + (p[0]*self.dx1[i])*(p[0]*self.dx1[i]) + \
                           (p[1]*self.dcolor[i])*(p[1]*self.dcolor[i]) + 2.0*p[0]*self.cov_m_s[i] - \
                            2.0*p[1]*self.cov_m_c[i] - 2.0*p[0]*p[1]*self.cov_s_c[i]



        residuals = self._resGen(covMatrix, res)

        return residuals

#
# residual for CMB
#
class CMBresiCal(ResiCal):

    def __init__(self,cosModel):

        self.cov = 1.0e-7 * np.array([[0.79039, -4.0042, 0.80608],
                                      [-4.0042,  66.950, -6.9243],
                                      [0.80608, -6.9243, 3.9712]])

        self.nu_cmb = np.array( [0.022065, 0.1199, 1.04131])
        self.numData = 3

        self.cosModel = cosModel

    def residual(self, p, fjac=None):

        """
        Parameter informations:
        # JLA nuiance parameters
        # ==========================
        # p[0]:alpha
        # p[1]:beta
        # p[2]:M
        # p[3]:DeltaM
        #
        # Cosmological parameters
        # ==========================
        # p[4]:H0
        # p[5]:Omega_bh2
        # p[6]:Omega_m
        # p[7]:Omega_rh2
        # p[8]:Omega_k
        # p[9]... other Cosmological parameters
        """

        h2 = ( p[4]/100.0 )**2
        omh2 = p[6]*h2        # dark matter only

        g1 = 0.0783*p[5]**(-0.238)/(1. + 39.5*p[5]**(0.763) )
        g2 = 0.560/(1.0 + 21.1 *p[5]**(1.81) )
        z_cmb = 1048.0*( 1.0 + 0.00124* p[5]**(-0.738) )* (1.0 + g1*omh2**(g2) )

        residuals = np.zeros(self.numData)

        residuals[0] =  p[5]      - self.nu_cmb[0]
        residuals[1] =  p[6]*h2   - self.nu_cmb[1]

        # remember to update the model parameter before calculation
        self.cosModel.parasUpdate(p)

        theta_cmb = self.cosModel.thetaFun(z_cmb)

        residuals[2] =  100.*theta_cmb - self.nu_cmb[2]

        res = self._resGen(self.cov, residuals)

        return res


class BAOresiCal(ResiCal):
    def __init__(self,cosModel):

        self.invCov = np.array([[4444.,  0.0,      0.0],
                                [0.0,    215156., 0.0],
                                [0.0,    0.0,      721487.]] )

        self.cov =  np.linalg.inv(self.invCov)

        self.zbao = np.array( [0.106, 0.35,   0.57]     )
        self.dbao = np.array( [0.336, 0.1126, 0.07315]  )
        self.numData = 3

        self.cosModel = cosModel


    def residual(self, p, fjac=None):

        """
        Parameter informations:
        # JLA nuiance parameters
        # ==========================
        # p[0]:alpha
        # p[1]:beta
        # p[2]:M
        # p[3]:DeltaM
        #
        # Cosmological parameters
        # ==========================
        # p[4]:H0
        # p[5]:Omega_bh2
        # p[6]:Omega_m
        # p[7]:Omega_rh2
        # p[8]:Omega_k
        # p[9]... other Cosmological parameters
        """


        h2 = ( p[4]/100.0 )**2
        omh2 = p[6]*h2        # dark matter only

        b1 = 0.313*omh2**(-0.419)*(1.0 + 0.607*omh2**(0.674) )
        b2 = 0.238*omh2**(0.223)
        z_drag = 1291.0*omh2**(0.251)*(1.0 + b1*p[5]**(b2) )/(1.0 + 0.659*omh2**(0.828) )

        residuals = np.zeros(self.numData)

        # remember to update the model parameter before calculation
        self.cosModel.parasUpdate(p)

        sound_h = self.cosModel.soundHor(z_drag)

        for i in range(self.numData):
            dv = self.cosModel.dzFun(self.zbao[i])
            residuals[i] = sound_h/dv - self.dbao[i]


        res = self._resGen(self.cov, residuals)

        return res


import random
# residual for gamma ray burst
class GamRayresiCal(ResiCal):

    def __init__(self, cosModel, DATA_DIR_GamRay):
        """
        cosModel: a object of class CosModel
        DATA_DIR_GamRay: the directory of Gamma ray data, which should look like this:

        """
    
                
        self.dataDir = DATA_DIR_GamRay
        self.cosModel = cosModel
        
        
        dataFile = open(DATA_DIR_GamRay+'/gammaray.txt', 'r')

        self.name      = np.array([])
        self.z         = np.array([])
        self.sb        = np.array([])
        self.sberr     = np.array([])
        self.ep        = np.array([])
        self.eperr     = np.array([])

        for line in dataFile:
            line = line.strip()
            if line.startswith('#'):
                temp=line.strip('#').split()
                self.snhead = temp
            else:
                temp=line.strip('\r\n').split()

                self.name      = np.append(self.name,       temp[0])
                self.z         = np.append(self.z,       float(temp[1]) )
                self.sb        = np.append(self.sb,      float(temp[2]) )
                self.sberr     = np.append(self.sberr,   float(temp[3]) )
                self.ep        = np.append(self.ep,      float(temp[4]) )
                self.eperr     = np.append(self.eperr,   float(temp[5]) )



        dataFile.close()

        self.numData = len(self.name)
        
        self.syserr = 0.7571
        
        # To get system residual when self.isAllRes = False
        self.isAllRes = True   
        self.ndata   = 100   
        self.rindex = np.array(random.sample(range(0,self.numData),self.ndata))
        



    def residualAll(self, p, fjac=None):

        """
        Parameter informations:
        # JLA nuiance parameters
        # ==========================
        # p[0]:alpha
        # p[1]:beta
        # p[2]:M
        # p[3]:DeltaM
        #
        # Cosmological parameters
        # ==========================
        # p[4]:H0
        # num = cosModel.nparams  (not include H0)
        # p[5]:Omega_bh2
        # p[6]:Omega_m
        # p[7]:Omega_rh2
        # p[8]:Omega_k
        # p[9]... other Cosmological parameters
        # ...
        # p[4+num]
        
        # Gammary buerst parameters
        # ==============================
        # p[4+num+1] : lambda   p[5+num]
        # p[4+num+2] : b        p[6+num] 
        """
        
        num = self.cosModel.nparams
        

        # remember to update the model parameter before calculation
        self.cosModel.parasUpdate(p)

        res = np.zeros(self.numData)
        err = np.zeros(self.numData)
        
        

        for i in range(self.numData):
            cdz = self.cosModel.covDis(self.z[i])
            
            powe  = (self.ep[i]/300.0)**p[6+num]
            sbolo = self.sb[i]*(PC2CM**2)*1.e-3 
            tmp = powe*(1.+self.z[i])/(4.*np.pi)/sbolo
            
            res[i] = 2.5*np.log10(tmp)+2.5*p[5+num]- ( 5.0*np.log10( (1.0+self.z[i])* cdz *3.e5/p[4] )+25.0)            
            
            sig1 =  p[6+num]*self.eperr[i]/self.ep[i] 
            sig2 =  self.sberr[i]/self.sb[i]
            
            staerr2 = (2.5/np.log(10))**2*(sig1**2+sig2**2)
            syserr2 = (2.5/np.log(10))**2*self.syserr**2
            
            err[i] = np.sqrt(staerr2 + syserr2)
        
            
            
        residuals = abs(res/err)
            

        return residuals
        
    def residual(self, p, fjac=None):
        
        if self.isAllRes:
            return self.residualAll(p)
        else:        
            return self.residualAll(p)[self.rindex]
        
        
        
    def __func(self, serr, args):
                
        
        
        self.syserr = serr
        
        
        res = self.residual(args[1])
        
        y = (res**2).sum()-args[0]        
        
        return y
        
        
    
    def syserrCal(self, freePara, p ):
        from scipy.optimize import root
        
        num     = len(p)-5
        dof     = self.ndata-num
        
        args = [ dof, p]
        #print args
        
        sol = root(self.__func, 0.3, args=args)
        
        self.syserr = sol.x
        
        return self.syserr
        
        
        
        
        
        
        
        
        

        
