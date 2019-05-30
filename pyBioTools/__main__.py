#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import argparse
import sys
import os
import json
import datetime
import textwrap

# Local imports
from pyBioTools import Alignment
# from pyBioTools.Annotation import Annotation
# from pyBioTools.Fastq import Fastq
# from pyBioTools.Fasta import Fasta
from pyBioTools.common import *
from pyBioTools import __version__ as package_version
from pyBioTools import __name__ as package_name
from pyBioTools import __description__ as package_description

#~~~~~~~~~~~~~~TOP LEVEL ENTRY POINT~~~~~~~~~~~~~~#
def main ():
    """
    Main entry point for command line interface
    """

    # Parser and subparsers for command
    parser = argparse.ArgumentParser (description=package_description)
    parser.add_argument("--version", action="version", version="{} v{}".format(package_name, package_version))
    subparsers = parser.add_subparsers (description="%(prog)s implements the following subcommands", dest="subcommands")
    subparsers.required = True

    #~~~~~Alignment suparser~~~~~#
    sp_a = subparsers.add_parser("Alignment", description="Aligned reads related functions")
    sp_a_subparsers = sp_a.add_subparsers (description="Alignment implements the following subcommands", dest="alignment_subcommands")
    sp_a_subparsers.required = True

    # alignment_index_reads subparser
    f = Alignment.index_reads
    sp_air = sp_a_subparsers.add_parser("index_reads", description=doc_func(f))
    sp_air.set_defaults(func=f)
    sp_air.add_argument("-b", "--bam_fn", **arg_opt(f, "bam_fn"))
    sp_air.add_argument("-m", "--min_mapq", **arg_opt(f, "min_mapq"))
    sp_air.add_argument("-u", "--keep_unmapped", **arg_opt(f, "keep_unmapped"))
    sp_air.add_argument("-s", "--keep_secondary", **arg_opt(f, "keep_secondary"))
    sp_air.add_argument("-p", "--keep_supplementary", **arg_opt(f, "keep_supplementary"))

    # alignment_downsample_reads
    f = Alignment.sample_reads
    sp_asr = sp_a_subparsers.add_parser("sample_reads", description=doc_func(f))
    sp_asr.set_defaults(func=f)
    sp_asr.add_argument("-b", "--bam_fn", **arg_opt(f, "bam_fn"))
    sp_asr.add_argument("-f", "--out_folder", **arg_opt(f, "out_folder"))
    sp_asr.add_argument("-p", "--out_prefix", **arg_opt(f, "out_prefix"))
    sp_asr.add_argument("-r", "--n_reads", **arg_opt(f, "n_reads"))
    sp_asr.add_argument("-s", "--n_samples", **arg_opt(f, "n_samples"))
    sp_asr.add_argument("--rand_seed", **arg_opt(f, "rand_seed"))

    # Add common group parsers
    for sp in [sp_air, sp_asr]:
        sp_verbosity = sp.add_mutually_exclusive_group()
        sp_verbosity.add_argument("-v", "--verbose", action="store_true", default=False, help="Increase verbosity (default: %(default)s)")
        sp_verbosity.add_argument("-q", "--quiet", action="store_true", default=False, help="Reduce verbosity (default: %(default)s)")

    # Parse args and call subfunction
    args = parser.parse_args()
    args.func(**vars(args))
