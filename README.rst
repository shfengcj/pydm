
=======================================
pydm: Python Data Mining for Astronomy
=======================================

pydm will be a Python module for data mining in astronomy and other fields. At present version, it provides two methods for cosmological model parameter fitting: the Levenberg-Marquardt (LM) technique and the Markov chain Monte Carlo (MCMC) approach. If you want to constrain a cosmological model with observations like supernovae, BAO, CMB, it is very easy and  convenient to use "pydm" package, with which you can:

  - Constrain parameters, find the best fitting values, and their covariance.
  - Make some contours between each of the parameters after marginalized others.
  - Plot some distribution, histograms of the given parameters.
  - Easily add new models.
  - Easily add new data sets.
  - Easily generalize to other fields.

Introduction
============
:Version: 1.0.0a
:Authors: Chao-Jun Feng <shfengcj@yahoo.com>
:License: "pydm" is released under the GNU 2.0 License.


Installation
============

This package uses distutils, which is the default way of installing python
modules.  **Before installation, make sure your system meets the prerequisites
listed in Dependencies, listed below.**

To install the  ``pydm`` package in your home directory, use::

  pip install pydm

To install from source, use::

  python setup.py install

You can specify an arbitrary directory for installation using::

  python setup.py install --prefix='/some/path'

To install system-wide on Linux/Unix systems::

  python setup.py build
  sudo python setup.py install


Dependencies
============

The authors of "pydm" would like to thank all the contributors of the following depending packages and data-base. Some much relevant packages and data are already included in "pydm" with some modifications for convenience. If the authors does't want their packages or data to be included in "pydm", just tell us: shfengcj AT yahoo.com, and we will remove them ASAP.

Required packages
-----------------
- Python >= 2.7.x (2.7.10 is preffered)
- Numpy  >= 1.9.2
- Scipy  >= 0.15.1
- Matplotlib >= 1.4.3

Required Astronomical data


Optional
------------
- agpy  >= 0.1.4
  agpy include the mpfit module which can perform Levenberg-Marquardt (LM) least-squares minimization, based on MINPACK-1. For convenience, mpfit is already included in the "pydm".
- emcee >= 2.1.0
  emcee is a stable, well tested Python implementation of the affine-invariant ensemble sampler for Markov chain Monte Carlo (MCMC) approach. For convenience, emcee is also included in the "pydm". Also, some modification was made to show the progress bar during the MCMC procedure.


TODO
=======
- develop more useful data mining function in Astronomy.
- ...


Links
=======

- Python: http://www.python.org
- Numpy: http://www.numpy.org
- Scipy: http://www.scipy.org
- Scikit-learn: http://scikit-learn.org
- Matplotlib: http://matplotlib.org
- AstroPy: http://www.astropy.org/
- agpy: http://packages.python.org/agpy
- emcee: http://dan.iel.fm/emcee/
- PyMC: http://pymc-devs.github.com/pymc/
- HEALPy: https://github.com/healpy/healpy>
...

Copyright (C) 2015- Chao-Jun Feng and contributors.
