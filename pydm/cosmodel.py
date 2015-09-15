# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 15:28:10 2015

@author: Chao-Jun Feng
"""


#===============================================================================
"""
                            cosmological model data struct
"""
#===============================================================================

"""
Class CosModel

attribute: name, paras

method:
    __invHubble : Calculate inverse of dimensionless Hubble parameter 1/E = H_0/H
    __soundHubble: Integrand of the sound horizon

    covDis: Hubble-independent comoving diameter distance D(z) = H0*(1+z)*D_A(z)
    soundHor: Sound horizon
    thetaFun, dzFun: two functions to calculate BAO distance

    parasUpdate: very important function to update new values of parameters


Users can build their own cosmological model by inherit this class

Example:

from numpy import sqrt, sin, sinh
from scipy.integrate import romberg


class lcdm(CosModel):

    def __init__(self, name = 'LCDM', isFlat = True, divmax = 15):
            self.name     = name
            self.numParas = 3
            self.isFlat   = isFlat

            self.ogh2     = 2.469e-5 # light
            self.Omega_bh2 = 0.022 # baryon * h**2

            self.Omega_r=(1.0+0.2271*3.046)*self.ogh2/( (0.6834)**2 )  # radiation
            self.Omega_M = 0.3  # dark matter + baryon

            if self.isFlat:
                self.Omega_k = 0.0 # spatial curvature
            else:
                self.Omega_k = -0.01

            # romberg intergral
            self.divmax = divmax


    def __invHubble(self,z):
        # Inverse of dimensionless Hubble 1/E = H_0/H


        omega_mgz = self.Omega_M*(1.+z)**3 + self.Omega_r*(1.+z)**4 + self.Omega_k*(1.+z)**2
        rhs = omega_mgz + (1. - self.Omega_M - self.Omega_r - self.Omega_k)
        invHubble = 1.0/sqrt(rhs)

        return invHubble

    def __soundHubble(self,a):

        omega_mgz = self.Omega_M*a**(-3) + self.Omega_r*a**(-4) + self.Omega_k*a**(-2)
        rhs = omega_mgz + (1. - self.Omega_M - self.Omega_r - self.Omega_k)

        invHubble = 1.0/sqrt(rhs)

        cs = 1.0/(3.0* (1.0 + 0.75*self.Omega_bh2*a/self.ogh2 ) )**(0.5)

        soundH  = cs*invHubble/a**2

        return soundH



    def covDis(self, z):
         #Hubble - independent Comoving diameter distance D(z) = H0*(1+z)*D_A(z)

        cdz = romberg(self.__invHubble, 0, z,divmax=self.divmax)

        if self.isFlat == False:

            if self.Omega_k > 0 :
                cdz =(1./sqrt(self.Omega_k))*sinh(sqrt(self.Omega_k)*cdz)
            if self.Omega_k < 0 :
                cdz = (1./sqrt(abs(self.Omega_k)))*sin(sqrt(abs(self.Omega_k))*cdz)

        return cdz

    def soundHor(self, z):
        # sound horizon

        a = 1.0/(1. + z)

        resh = romberg(self.__soundHubble, 1.0e-8,a,divmax=self.divmax)
        return resh

    def thetaFun(self,z):
        soundh  = self.soundHor(z)
        cdz     = self.covDis(z)
        theta   = soundh/cdz

        return theta

    def dzFun(self,z):
        cdz     = self.covDis(z)
        invh    = self.__invHubble(z)
        dv      = (cdz**2 * z *invh)**(1.0/3.0)

        return dv

    def parasUpdate(self,p):

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

        h2 = ( p[4]/100.0 )**2

        self.Omega_bh2  = p[5]
        self.Omega_M    = p[6] + p[5]/h2 # Omega_M = dark matter + baryon
        self.Omega_r    = p[7]/h2
        self.Omega_k    = p[8]


        return

"""
class CosModel(object):
    def __init__(self, name, nparams, paras=[]):
        """
        name: the name of cosmological model
        nparams: the number of params
        nparams: the list of parameters of the model

        one could also define parameters here: e.g.

        self.Omega_m = 0.2
        ...
        """
        self.name    = name
        self.nparams = nparams


    def __invHubble(self,z):
        """
        z: the redshift
        """
        res = 0
        return res


    def __soundHubble(self,a):
        """
        a: the scale factor
        """
        res = 0
        return res


    def covDis(self, z):
        res = 0
        return res

    def soundHor(self, z):
        res = 0
        return res

    def thetaFun(self,z):
        res = 0
        return res

    def dzFun(self,z):
        res = 0
        return res

    def parasUpdate(p):
        """
        p: new parameter list or array

        e.g.
        self.Omega_m = p[0]
        """

        pass
        return
