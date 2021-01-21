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
from pyBioTools import Seqsum
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

    # Alignment Reads_index subparser
    f = Alignment.Reads_index
    sp_a_ir = sp_a_subparsers.add_parser("Reads_index", description=doc_func(f))
    sp_a_ir.set_defaults(func=f)
    arg_from_docstr(sp_a_ir, f, "input_fn", "i")
    arg_from_docstr(sp_a_ir, f, "skip_unmapped", "u")
    arg_from_docstr(sp_a_ir, f, "skip_secondary", "s")
    arg_from_docstr(sp_a_ir, f, "skip_supplementary", "p")

    f = Alignment.Reads_sample
    sp_a_sr = sp_a_subparsers.add_parser("Reads_sample", description=doc_func(f))
    sp_a_sr.set_defaults(func=f)
    arg_from_docstr(sp_a_sr, f, "input_fn", "i")
    arg_from_docstr(sp_a_sr, f, "output_folder" , "o")
    arg_from_docstr(sp_a_sr, f, "output_prefix" , "p")
    arg_from_docstr(sp_a_sr, f, "n_reads" , "r")
    arg_from_docstr(sp_a_sr, f, "n_samples" , "s")
    arg_from_docstr(sp_a_sr, f, "rand_seed")

    f = Alignment.References_sample
    sp_a_rs = sp_a_subparsers.add_parser("References_sample", description=doc_func(f))
    sp_a_rs.set_defaults(func=f)
    arg_from_docstr(sp_a_rs, f, "input_fn", "i")
    arg_from_docstr(sp_a_rs, f, "output_fn", "o")
    arg_from_docstr(sp_a_rs, f, "selected_reads_fn", "s")
    arg_from_docstr(sp_a_rs, f, "frac_reads", "f")
    arg_from_docstr(sp_a_rs, f, "min_reads_ref", "r")
    arg_from_docstr(sp_a_rs, f, "sorting_threads", "t")
    arg_from_docstr(sp_a_rs, f, "rand_seed")

    f = Alignment.Filter
    sp_a_fr = sp_a_subparsers.add_parser("Filter", description=doc_func(f))
    sp_a_fr.set_defaults(func=f)
    arg_from_docstr(sp_a_fr, f, "input_fn", "i")
    arg_from_docstr(sp_a_fr, f, "output_fn", "o")
    arg_from_docstr(sp_a_fr, f, "selected_reads_fn", "l")
    arg_from_docstr(sp_a_fr, f, "skip_unmapped", "u")
    arg_from_docstr(sp_a_fr, f, "skip_secondary", "s")
    arg_from_docstr(sp_a_fr, f, "skip_supplementary", "p")
    arg_from_docstr(sp_a_fr, f, "index_reads", "x")
    arg_from_docstr(sp_a_fr, f, "orientation", "t")
    arg_from_docstr(sp_a_fr, f, "min_read_len", "r")
    arg_from_docstr(sp_a_fr, f, "min_align_len", "a")
    arg_from_docstr(sp_a_fr, f, "min_mapq", "m")
    arg_from_docstr(sp_a_fr, f, "min_freq_identity", "f")
    arg_from_docstr(sp_a_fr, f, "select_ref_fn")
    arg_from_docstr(sp_a_fr, f, "exclude_ref_fn")

    f = Alignment.To_fastq
    sp_a_tf = sp_a_subparsers.add_parser("To_fastq", description=doc_func(f))
    sp_a_tf.set_defaults(func=f)
    arg_from_docstr(sp_a_tf, f, "input_fn", "i")
    arg_from_docstr(sp_a_tf, f, "output_r1_fn" , "1")
    arg_from_docstr(sp_a_tf, f, "output_r2_fn" , "2")
    arg_from_docstr(sp_a_tf, f, "ignore_paired_end" , "s")

    f = Alignment.Split
    sp_a_st = sp_a_subparsers.add_parser("Split", description=doc_func(f))
    sp_a_st.set_defaults(func=f)
    arg_from_docstr(sp_a_st, f, "input_fn", "i")
    arg_from_docstr(sp_a_st, f, "output_dir" , "o")
    arg_from_docstr(sp_a_st, f, "n_files" , "n")
    arg_from_docstr(sp_a_st, f, "output_fn_list" , "l")
    arg_from_docstr(sp_a_st, f, "index", "x")

    #~~~~~Fastq suparser~~~~~#
    sp_fq = subparsers.add_parser("Fastq", description="Fastq format related functions")
    sp_fq_subparsers = sp_fq.add_subparsers (description="Fastq implements the following subcommands", dest="fastq_subcommands")
    sp_fq_subparsers.required = True

    # Fastq Filter subparser
    f = Fastq.Filter
    sp_fq_fr = sp_fq_subparsers.add_parser("Filter", description=doc_func(f))
    sp_fq_fr.set_defaults(func=f)
    arg_from_docstr(sp_fq_fr, f, "input_fn" , "i")
    arg_from_docstr(sp_fq_fr, f, "output_fn" , "o")
    arg_from_docstr(sp_fq_fr, f, "min_len" , "l")
    arg_from_docstr(sp_fq_fr, f, "min_qual" , "u")
    arg_from_docstr(sp_fq_fr, f, "remove_duplicates" , "r")
    arg_from_docstr(sp_fq_fr, f, "qual_offset" , "f")

    #~~~~~Seqsum suparser~~~~~#
    sp_ss = subparsers.add_parser("Seqsum", description="Sequencing summary format related functions")
    sp_ss_subparsers = sp_ss.add_subparsers (description="Seqsum implements the following subcommands", dest="fastq_subcommands")
    sp_ss_subparsers.required = True

    # Fastq Filter subparser
    f = Seqsum.Merge
    sp_ss_mg = sp_ss_subparsers.add_parser("Merge", description=doc_func(f))
    sp_ss_mg.set_defaults(func=f)
    arg_from_docstr(sp_ss_mg, f, "input_fn" , "i")
    arg_from_docstr(sp_ss_mg, f, "output_fn" , "o")
    arg_from_docstr(sp_ss_mg, f, "old_filename_synthax")

    # Add common group parsers
    for sp in [sp_a_ir, sp_a_sr, sp_a_rs, sp_a_fr, sp_a_tf, sp_fq_fr, sp_a_st, sp_ss_mg]:
        sp.add_argument("-v", "--verbose", action="store_true", default=False, help="Increase verbosity (default: %(default)s)")
        sp.add_argument("-q", "--quiet", action="store_true", default=False, help="Reduce verbosity (default: %(default)s)")
        sp.add_argument("--progress", action="store_true", default=False, help="Display a progress bar")

    # Parse args and call subfunction
    args = parser.parse_args()
    args.func(**vars(args))
