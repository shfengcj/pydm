# -*- coding: utf-8 -*-


# Genearl setting

"""
Some constants
"""
divMax = 15  # for romberg integral
ogh2     = 2.469e-5   # light
orh2     = ogh2*(1.+0.2271*3.046) # radiation = light + neutrino


"""
The directory of JLA data, which should look like this:

JLA_DIR
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
# Modify the floowing line to you data directory.
JLA_DIR = '/Users/chaojun/Documents/Research/2015/grb/pycode/data/jla'


"""
Data setting
"""
use_sn_data     = True
use_cmb_data    = True
use_bao_data    = True



"""
MCMC work mode set
"""

# Set True to use Levenberg-Marquardt least-squares minimization
LMmin = True

# Set True to use MCMC simulation
MCMCsim = True

# Set True to plot figgures
FIGplot = True

# Result save directory
OUT_DIR = './output/'
# Save the fitting result
LMminSave = True
LMminFile = OUT_DIR + "lcdm.LM.res"

MCMCsimSave = True
MCMCsimFile = OUT_DIR + "lcdm.MCMC.res"

# Set up the sampler.
nwalkers   = 200
iterations = 300
burnIn     = 100

# Save the chain
chainSave = True
chainFile = OUT_DIR + "lcdm.chain.res"

# contour properties
# Plot options ---->
# plot mode
#   0: triangle,
#   1: hist2d between two parameters,
#   2: hist for one parameter
plotMode = 0
#plotModeDic = {0:'triangle', 1:'contours', 2:'histogram'}

# smooth or not
isSmooth    = True
isSmooth1d  = True

# color, bin, level for contour and hist
color       = 'b'
bins        = 30
levels      = [0.683, 0.954]#, 0.997]
isFillCont  = False
isNormHist  = False

# datapoint color
dataKwargs  = {"color":'r'}
