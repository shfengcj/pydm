#! /usr/bin/env python
# -*- coding: utf-8 -*-

try:
    from setuptools import setup
    setup
except ImportError:
    from distutils.core import setup


DESCRIPTION = "A Python Module for Data Mining in Astronomy and other Fields"
LONG_DESCRIPTION = open('README.rst').read()
NAME = "pydm"
AUTHOR = "Chao-Jun Feng"
AUTHOR_EMAIL = "shfengcj@yahoo.com"
MAINTAINER = "Chao-Jun Feng"
MAINTAINER_EMAIL = "shfengcj@yahoo.com"
URL = ''
DOWNLOAD_URL = ''
LICENSE = 'GNU 2.0'

import pydm
VERSION = pydm.__version__

setup(name=NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      maintainer=MAINTAINER,
      maintainer_email=MAINTAINER_EMAIL,
      url=URL,
      download_url=DOWNLOAD_URL,
      license=LICENSE,
      packages=['pydm',
                'pydm.emcee',
                'pydm.mpfit',
                'pydm.utility',
            ],
      classifiers=[
            'Development Status :: 1 - Alpha',
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU 2.0 License',
            'Natural Language :: English',
            'Programming Language :: Python :: 2.7.10',
            'Topic :: Scientific/Engineering :: Astronomy'],
    )
