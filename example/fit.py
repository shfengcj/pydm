# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 06:53:32 2015

@author: Chao-Jun Feng
"""


from __future__ import print_function

import numpy as np
import os


from lcdm   import lcdm
from pydm   import JLAresiCal, CMBresiCal, BAOresiCal
from pydm   import mcmcEngine

from genset import *

try:
    os.makedirs(OUT_DIR)
except os.error:
    pass

#========================================================================

labelPara = ['$\\alpha$','$\\beta$', '$M$', '$\Delta~M$', \
             '$H_0$', '$\Omega_{b}h^2$','$\Omega_{m}$', \
             '$\Omega_{r}h^2$', '$\Omega_{k}$']
resobj=[]

if LMmin or MCMCsim:

    print ('Reading data ...')

    # Cosmological model
    model = lcdm(divmax = divMax)

    if use_sn_data : resobj.append( JLAresiCal(cosModel = model, DATA_DIR_JLA = JLA_DIR) )
    if use_cmb_data: resobj.append( CMBresiCal(cosModel = model) )
    if use_bao_data: resobj.append( BAOresiCal(cosModel = model) )


# start point for LM minization
startp  = np.array([0.141 , 3.099 , -19.10 ,-0.07 , \
                    68.34 , 0.0221, 0.305,          \
                    orh2,   0.001 ],dtype='float64')

# prior set
pmin = np.array([0.01, 0.9, -25.0, -1.0, 20.0, 0.005, 0.1, 0.5*orh2, -0.5])
pmax = np.array([2.00, 4.6, -15.0,  0.1, 100., 0.100, 1.0, 1.5*orh2,  0.5])
plimit = zip(pmin, pmax)


# Which parameters to be fitted (plotted)
# 1 for fit(plot), 0 for fixed(not plot)
# Note: plotPara must in freePara
freePara  = [1,1,1,1,1,1,1,0,0]
plotPara  = [1,1,0,0,1,1,1,0,0]


#========================================================================

sim = mcmcEngine(startp = startp, freePara = freePara, labelPara = labelPara, resobj = resobj,
             nwalkers = nwalkers, iterations = iterations ,burnIn = burnIn)
if LMmin:
    sim.LMmin( LMminSave = LMminSave, LMminFile = LMminFile)
if MCMCsim:
    sim.MCMCsim(plimit = plimit, MCMCsimSave = MCMCsimSave, MCMCsimFile = MCMCsimFile,
            chainSave   = chainSave,   chainFile   = chainFile)

if FIGplot:
   sim.FIGplot(plotPara = plotPara, sampleFile = None,
           plotMode = plotMode,   isSmooth = isSmooth, isSmooth1d = isSmooth1d,
           color    = color, bins     = bins,   levels     = levels,
           isFillCont  = isFillCont, isNormHist  = isNormHist, saveDir = OUT_DIR,
           preFig = 'lcdm-', sufFig = '.png', dataKwargs = dataKwargs)
#
# if FIGplot:
#     sim.FIGplot(plotPara = plotPara, sampleFile = chainFile,
#             plotMode = plotMode,   isSmooth = isSmooth, isSmooth1d = isSmooth1d,
#             color    = color, bins     = bins,   levels     = levels,
#             isFillCont  = isFillCont, isNormHist  = isNormHist, saveDir = OUT_DIR,
#             preFig = 'lcdm-', sufFig = '.png', dataKwargs = dataKwargs, truthFile = LMminFile)

print("Jobs are all finished.")
