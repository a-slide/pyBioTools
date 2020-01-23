#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup

# Define package info
name = "pyBioTools"
version = '0.2.0.5'
description = "pyBioTools is a collection of tools to manipulate biological sequences"
with open("README.md", "r") as fh:
    long_description = fh.read()

# Collect info in a dictionnary for setup.py
setup(
    name = name,
    description = description,
    version = version,
    long_description = long_description,
    long_description_content_type="text/markdown",
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
        'tqdm>=4.32.1',
        'numpy>=1.16.4',
        'pysam>=0.15.2',
        'pandas>=0.24.2'],
    packages = [name],
    entry_points = {'console_scripts': ['pyBioTools=pyBioTools.__main__:main']})
