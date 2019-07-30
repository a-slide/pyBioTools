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
from pyBioTools import Fastq
# from pyBioTools.Annotation import Annotation
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
    sp_a_ir = sp_a_subparsers.add_parser("index_reads", description=doc_func(f))
    sp_a_ir.set_defaults(func=f)
    sp_a_ir.add_argument("-b", "--bam_fn", **arg_opt(f, "bam_fn"))
    sp_a_ir.add_argument("-m", "--min_mapq", **arg_opt(f, "min_mapq"))
    sp_a_ir.add_argument("-u", "--keep_unmapped", **arg_opt(f, "keep_unmapped"))
    sp_a_ir.add_argument("-s", "--keep_secondary", **arg_opt(f, "keep_secondary"))
    sp_a_ir.add_argument("-p", "--keep_supplementary", **arg_opt(f, "keep_supplementary"))

    # alignment_downsample_reads
    f = Alignment.sample_reads
    sp_a_sr = sp_a_subparsers.add_parser("sample_reads", description=doc_func(f))
    sp_a_sr.set_defaults(func=f)
    sp_a_sr.add_argument("-b", "--bam_fn", **arg_opt(f, "bam_fn"))
    sp_a_sr.add_argument("-f", "--out_folder", **arg_opt(f, "out_folder"))
    sp_a_sr.add_argument("-p", "--out_prefix", **arg_opt(f, "out_prefix"))
    sp_a_sr.add_argument("-r", "--n_reads", **arg_opt(f, "n_reads"))
    sp_a_sr.add_argument("-s", "--n_samples", **arg_opt(f, "n_samples"))
    sp_a_sr.add_argument("--rand_seed", **arg_opt(f, "rand_seed"))

    #~~~~~Fastq suparser~~~~~#
    sp_fq = subparsers.add_parser("Fastq", description="Fastq format related functions")
    sp_fq_subparsers = sp_fq.add_subparsers (description="Fastq implements the following subcommands", dest="fastq_subcommands")
    sp_fq_subparsers.required = True

    # Fastq filter_reads subparser
    f = Fastq.filter_reads
    sp_fq_fr = sp_fq_subparsers.add_parser("filter_reads", description=doc_func(f))
    sp_fq_fr.set_defaults(func=f)
    sp_fq_fr.add_argument("-s", "--src_fn", **arg_opt(f, "src_fn"))
    sp_fq_fr.add_argument("-d", "--dest_fn", **arg_opt(f, "dest_fn"))
    sp_fq_fr.add_argument("-l", "--min_len", **arg_opt(f, "min_len"))
    sp_fq_fr.add_argument("-u", "--min_qual", **arg_opt(f, "min_qual"))
    sp_fq_fr.add_argument("-r", "--remove_duplicates", **arg_opt(f, "remove_duplicates"))
    sp_fq_fr.add_argument("-o", "--qual_offset", **arg_opt(f, "qual_offset"))

    # Add common group parsers
    for sp in [sp_a_ir, sp_a_sr, sp_fq_fr]:
        sp_verbosity = sp.add_mutually_exclusive_group()
        sp_verbosity.add_argument("-v", "--verbose", action="store_true", default=False, help="Increase verbosity (default: %(default)s)")
        sp_verbosity.add_argument("-q", "--quiet", action="store_true", default=False, help="Reduce verbosity (default: %(default)s)")

    # Parse args and call subfunction
    args = parser.parse_args()
    args.func(**vars(args))
