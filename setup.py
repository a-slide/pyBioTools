#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup
import pyBioTools as pbt

# Collect info in a dictionnary for setup.py
setup(
    name =  pbt.__name__,
    version = pbt.__version__,
    description = pbt.__description__,
    url = "https://github.com/a-slide/pyBioTools",
    author = 'Adrien Leger',
    author_email = 'aleg@ebi.ac.uk',
    license = 'GPLv3',
    python_requires ='>=3.5',
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3'],
    install_requires = [
        'tqdm==4.32.1',
        'numpy==1.16.4',
        'pysam==0.15.2',
        'pandas==0.24.2'],
    packages = [pbt.__name__],
    entry_points = {'console_scripts': ['pyBioTools=pyBioTools.__main__:main']})
