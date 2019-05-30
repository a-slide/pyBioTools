# -*- coding: utf-8 -*-

# Define self package variable
__version__ = '0.0.0.1a'
__all__ = [""]
__description__="""
pyBioTools is a collection of tools to manipulate biological sequences
"""
# Collect info in a dictionnary for setup.py
setup_dict = {
    "name": __name__,
    "version": __version__,
    "description": __description__,
    "url": "https://github.com/a-slide/pyBioTools",
    "author": 'Adrien Leger',
    "author_email": 'aleg {at} ebi.ac.uk',
    "license": 'GPLv3',
    "python_requires":'>=3.4',
    "classifiers": [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3'],
    "install_requires": ['tqdm==4.32.1', 'numpy==1.16.4', 'pysam==0.15.2', 'pandas==0.24.2'],
    "packages": [__name__],
    "entry_points": {'console_scripts': ['pyBioTools=pyBioTools.__main__:main']}
}
