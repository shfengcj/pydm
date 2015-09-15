# -*- coding: utf-8 -*-
"""
Created on Sat Jul 11 08:21:34 2015

@author: Chao-Jun Feng

define mcmcEngine class

"""

from __future__ import print_function

import numpy as np
import matplotlib.pyplot as pl
import sys
import copy

from .mpfit import mpfit
from .emcee import EnsembleSampler
from .utility import progress_bar as probar
from .utility import triangle



class mcmcEngine(object):
    """
    class mcmcEngine

    methods:

    __init__   : initial the class
    __residual : residual function
    __lnlike   : log likelyhood function
    __lnprior  : log prior function
    __lnprob   : log post probability function

    readSample : read samples from existing chain files
    LMmin      : LM minization based on mpfit
    MCMCsim    : MCMC simulation based on emcee
    FIGplot    : plot figures
    """

    def __init__(self, startp, freePara, labelPara, resobj,
                 nwalkers = 200, iterations = 300 ,burnIn = 100):
        """
        __init__ function:

         :param startp: start position of parameters
            e.g. startp = np.array([0.141 , 3.099 , -19.10 ,-0.07 , \
                                    68.34 , 0.0221, 0.305,          \
                                    orh2,   0.001 ],dtype='float64')
         :param freePara: which parameters to be fitted, 1 for fit, 0 for fixed
            e.g. freePara    = [1,1,1,1,1,1,1,0,0]
         :param labelPara: labels of parameters
            e.g. labelPara = ['$\\alpha$','$\\beta$', '$M$', '$\Delta M$', \
                         '$H_0$', '$\Omega_{b}h^2$','$\Omega_{m}$', \
                         '$\Omega_{r}h^2$', '$\Omega_{k}$']
         :param resobj: list of residual objects

        """

        self.startp = startp
        self.n      = len(startp) # n: nuber of parameters

        self.freePara   = freePara
        self.labelPara  = labelPara

        self.freeIndex   = [i for i, j in enumerate(self.freePara) if j>0]
        self.freeLabel   = [i for i, j in zip(self.labelPara, self.freePara) if j > 0 ]


        # Start position for MCMC simulation
        self.p0  = self.startp[self.freeIndex]


        self.resobj = resobj

        # set the sampler
        self.nwalkers   = nwalkers
        self.iterations = iterations
        self.burnIn     = burnIn



    # Residual function
    def __residual(self, p,  fjac=None):

        res = np.array([])

        for obj in self.resobj:
            tmp = obj.residual(p)
            res = np.append(res, tmp)

        status = 0
        return [status, res]
        

    # Log likelyhood function
    def __lnlike(self,p):

        status, res = self.__residual(p)
        llike = -0.5*np.sum(res**2)

        return llike


    # Define the probability function as likelihood * prior.
    def __lnprior(self, p):

        if self.plimit is None:
            return 0.0

        for i, (j, k) in enumerate(self.plimit):
            if j < p[i] < k:
                pass
            else:
                return -np.inf

        return 0.0




    def __lnprob(self, para):
        """
        log post probability function
        para: free paramters that may be changed during MCMC simulation
        while, the following p is the total paramters of the model which includes the free ones and fixed ones.
        """

        p = copy.deepcopy(self.startp)

        p[self.freeIndex] = para

        lp =  self.__lnprior(p)
        
        if not np.isfinite(lp):
            return -np.inf
        return lp + self.__lnlike(p)


    def readSample(self, fileName):

        # n: number of parameters
        n = len(self.freeIndex)
        nline = self.iterations*self.nwalkers
        sample = np.zeros((nline, n))

        try:
            sampleFile = open(fileName, 'r')

            for i, line in enumerate(sampleFile):
                temp=line.strip('\r\n').split()
                oneSample = np.array([float(val) for val in temp])
                sample[i,:] = oneSample
            if i <> (nline-1) or len(temp) <> n:
                raise mcmcErr("{0} is not the correct file, please check whether it has {1} rows and {2} columns".format(fileName, nline, n))
        except:
            #print("can not read {0}, plese check !".format(fileName))
            return None
        finally:
            sampleFile.close()

        return sample[self.burnIn*self.nwalkers:, :]

    def readParaVal(self, fileName):
        paraVal = np.array([])
        try:
            f = open(fileName, 'r')

            for line in f:
                line = line.strip()
                if line.startswith('#'):
                    continue
                else:
                    temp=line.strip('\r\n').split()
                    if temp == []:
                        continue
                    paraVal = np.append(paraVal,  float(temp[2]) )

            if len(paraVal) <> self.n:
                 print("{0} is not the correct file".format(fileName))

        except:
            raise mcmcErr("can not read {0}!".format(fileName))
        finally:
            f.close()



        return paraVal



    def LMmin(self,  parinfo = None, LMminSave = True, LMminFile = None):

        # Set up the parameters
        if parinfo is None:
            self.parbase={'value':0., 'fixed':0, 'limited':[0,0], 'limits':[0.,0.]}
            self.parinfo=[]

            for i in range(self.n):
                self.parinfo.append(copy.deepcopy(self.parbase))
            for i in range(self.n):
                self.parinfo[i]['value']=self.startp[i]
                self.parinfo[i]['fixed']= 1- self.freePara[i]
        else:
            self.parinfo = parinfo
            
        # Do the LM process
        m = mpfit(self.__residual, self.startp, parinfo=self.parinfo )

        if (m.status <= 0):
            print ('error message = ', m.errmsg)

        chisq = (self.__residual(m.params)[1]**2).sum()
        llike = self.__lnlike(m.params)

        # To show the results:
        print ('Best fitting values: ', m.params)
        print ('1 sigma error : ', m.perror)
        print ('chi2 min = ', chisq)
        print ('lnlike = ', llike)

        if LMminSave:
            try:
                f = open(LMminFile, "w")
                f.write("# LM fitting result: \n")
                f.write('# chi2 min = ' + str(chisq) +'/' + str(m.dof)+" DOF\n")
                f.write('# lnlike = ' + str(llike)+"\n" )
                f.write('# ==========================\n')
                for i in range(len(m.params)):
                    f.write("\n"+self.labelPara[i].strip('$\\') + ' = ' + str(m.params[i]) \
                    +' +- '+str(m.perror[i]) )
                    f.write('\n')
            except:
                print("can not read {0}!".format(LMminFile))
            finally:
                f.close()


        self.p0 = m.params[self.freeIndex]

        self.bestPara = copy.deepcopy(self.p0)

    def MCMCsim(self, plimit = None, MCMCsimSave =True, MCMCsimFile = None, chainSave = True, chainFile = None):

        self.ndim     = len(self.p0)

        self.plimit   = plimit

        # The positions of all walkers
        if not 0 in self.p0:
            diff = 0.01*self.p0
        else:
            diff = 0.01


        self.pos0 = [self.p0 + diff*np.random.randn(self.ndim) for i in range(self.nwalkers)]

        print("Get sampler...")

        self.sampler = EnsembleSampler(self.nwalkers, self.ndim, self.__lnprob)

        print("Running MCMC...")
        sys.stdout.flush()

        # Set the progressing bar
        pbar = probar(self.iterations)
        self.pbar = pbar


        if chainSave:
            # Sample, outputting to a file
            try:
                f = open(chainFile, "w")
                for pos, prob, rstate, it in self.sampler.sample(self.pos0, iterations=self.iterations):
                    pbar.update(it)
                    
                    # Write the current position to a file, one line per walker
                    f.write("\n".join(["\t".join([str(q) for q in p]) for p in pos]))
                    f.write("\n")
            except:
                print("can not read {0}!".format(chainFile))
            finally:
                f.close()

        # Print out the mean acceptance fraction.
        print("Mean acceptance fraction:", np.mean(self.sampler.acceptance_fraction))

        # Estimate the integrated autocorrelation time for the time series in each parameter.
        print("Autocorrelation time:", self.sampler.get_autocorr_time())

        self.MCMCsamp = self.sampler.chain[:, self.burnIn:, :].reshape((-1, self.ndim ))

        # e.g. oneSig = [50.*(1-0.683), 50., 50.*(1+0.683) ]
        CL = [0.683, 0.954, 0.997]
        SigErr = zip((1.-np.array(CL))*50.,np.ones(len(CL))*50. , (1.+np.array(CL))*50.  )

        # v: best, +1s, -1s
        self.MCMCresult = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                         zip(*np.percentile(self.MCMCsamp,SigErr[0], axis=0)))


        if self.bestPara is None:
            self.bestPara = np.array([len(self.freeIndex)])

        print("MCMC result: ")
        print("===============\n")
        for i in range(len(self.freeIndex)):
            self.bestPara[i] = self.MCMCresult[i][0]

            print("""{0} = {1[0]}+{1[1]} -{1[2]}
            """.format(self.freeLabel[i].strip('$\\'), self.MCMCresult[i]))
            print('\n')

        if MCMCsimSave:
            lnprob = self.sampler.lnprobability.max()
            chisq  = -2.0*lnprob
            try:
                f = open(MCMCsimFile, "w")
                f.write('MCMC result with {0} walkers {1} steps : \n'.format(self.nwalkers,self.iterations))
                f.write('chi2 min = ' + str(chisq) +"\n")
                f.write('lnlike = ' + str(lnprob)+"\n" )
                f.write("===============\n")
                for i in range(len(self.freeIndex)):
                    f.write("""\n {0} = {1[0]}+{1[1]} -{1[2]}
                    """.format(self.freeLabel[i].strip('$\\'), self.MCMCresult[i]))
                    f.write('\n')
                f.write("===============\n")
                f.write("\n Mean acceptance fraction: "+str(np.mean(self.sampler.acceptance_fraction)) )
                f.write("\n Autocorrelation time:"+ str(self.sampler.get_autocorr_time()))
            except:
                print("can not read {0}!".format(MCMCsimFile))
            finally:
                f.close()


        print("MCMC simulation done.")


    def FIGplot(self, plotPara, sampleFile = None,
                plotMode = 0,   isSmooth = True, isSmooth1d = True,
                color    = 'b', bins     = 30,   levels     = [0.683, 0.954],
                isFillCont  = False, isNormHist  = False, saveDir = './output/',
                preFig = 'lcdm-', sufFig = '.png', dataKwargs = None, truths = None, truthFile=None):
        """
        Plot options
         - plotMode
            0: triangle,
            1: hist2d between two parameters,
            2: hist for one parameter
         - plotPara : Which paramters to be plotted, 1 for plot 0 for not plot
            e.g. plotPara  = [1,1,0,0,1,1,1,0,0]
         - sampleFile : data sample from MCMC simulation or exiting chain file provided
         - isSmooth : smooth the contours
         - isSmooth1d : smooth the 1d histogram
         - color: contour color
         - bins: histogram bins
         - levels: contour level, e.g. levels     = [0.683, 0.954, 0.997] for 1,2,3 sigma CL
         - isFillCont: True for filled contours
         - isNormHist: True for nomalized histogram

         - **dataKwargs: datapoint properties
                e.g. dataKwargs  = {"color":'r'}

        """
        self.plotModeDic = {0:'triangle', 1:'contours', 2:'histogram'}

        self.plotPara    = plotPara

        assert len(self.plotPara) == self.n, (
            "Parameters plotting vector must has the same dimension as  "
            "that of the start parameters!")
        for i in range(self.n):
            assert self.freePara[i] >= self.plotPara[i], (
            "Parameters to be plotted must be a free parameter!")

        self.plotIndex   = [i for i, j in zip(self.plotPara, self.freePara)  if j>0]
        self.plotIndex   = [i for i, j in enumerate(self.plotIndex)          if j>0]
        self.plotLabel   = [i for i, j in zip(self.labelPara, self.plotPara) if j>0]
        # prefix and suffix of figure names
        self.preFig = preFig
        self.sufFig = sufFig

        print ('Plot process with {0} mode :'.format(self.plotModeDic[plotMode]))

        if sampleFile is None:
            try:
                samp = self.MCMCsamp[:,self.plotIndex]
                truths = self.bestPara[self.plotIndex]
            except:
                raise mcmcErr('Please run LMmin or MCMCsim to get samples')
        else:
            # sample from saved file
            try:
                samp = self.readSample(fileName = sampleFile)[:,self.plotIndex]
            except:
                raise mcmcErr("Samples are not correctly read! Check the chain file!")


        if  truths is None:
            if truthFile <> None:

                truths = self.readParaVal(fileName = truthFile)[self.freeIndex][self.plotIndex]
            else:
                print("WARNING: Truth values are not given.")

        if plotMode == 0:

            figTri  = triangle.corner(samp, labels        = self.plotLabel,
                                            truths        = truths,
                                            bins          = bins,
                                            levels        = levels,
                                            color         = color,
                                            smooth        = isSmooth,
                                            smooth1d      = isSmooth1d,
                                            fill_contours = isFillCont,
                                            data_kwargs   = dataKwargs)
            figTri.savefig(saveDir+preFig+"triangle" + sufFig)

        elif plotMode == 1:

            n_plot = len(self.plotIndex)

            for i in range(n_plot-1):
                for j in range(i+1, n_plot):
                    #print( i ,j)

                    figCon = triangle.hist2d(samp[:,i], samp[:,j],
                                             levels = levels,
                                             color = color,
                                             smooth = isSmooth,
                                             fill_contours= isFillCont,
                                             data_kwargs = dataKwargs)

                    xl = pl.gca().set_xlabel(self.plotLabel[i])
                    yl = pl.gca().set_ylabel(self.plotLabel[j])

                    fileName = "{0}-{1}-contour".format(self.plotLabel[i].strip('$\\'),
                                                        self.plotLabel[j].strip('$\\'))

                    pl.savefig(saveDir+preFig+fileName+sufFig)

        elif plotMode == 2:
            n_plot = len(self.plotIndex)
            for i in range(n_plot):
                figHist =  pl.hist(samp[:,i], bins =bins, color=color, normed = isNormHist)
               # pl.xlabel = plotLabel[i]
                xl = pl.gca().set_xlabel(self.plotLabel[i])
                fileName = "{0}-hist".format(self.plotLabel[i].strip('$\\') )
                pl.savefig(saveDir+preFig+fileName+sufFig)
                pl.cla()

        else:
            print("Figures are not plotted, please choice a correct plot mode!")

class mcmcErr(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return repr(self.msg)
