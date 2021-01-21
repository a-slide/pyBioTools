#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup

# Long description from README file
with open("README.md", "r") as fh:
    long_description = fh.read()

# Collect info in a dictionary for setup.py
setup(
    name="pyBioTools",
    description="pyBioTools is a collection of python tools to manipulate biological sequences",
    version="0.2.7",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/a-slide/pyBioTools",
    author="Adrien Leger",
    author_email="aleg@ebi.ac.uk",
    license="GPLv3",
    python_requires=">=3.6",
    classifiers=["Development Status :: 3 - Alpha", "Intended Audience :: Science/Research", "Topic :: Scientific/Engineering :: Bio-Informatics", "License :: OSI Approved :: GNU General Public License v3 (GPLv3)", "Programming Language :: Python :: 3"],
    install_requires=["tqdm>=4.51.0", "numpy>=1.19.4", "pysam>=0.16.0", "pandas>=1.1.4",  "colorlog>=4.1.0"],
    packages=["pyBioTools"],
    entry_points={"console_scripts": ["pyBioTools=pyBioTools.__main__:main"]},
)
