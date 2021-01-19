# -*- coding: utf-8 -*-

# Standard library imports
from collections import *
import gzip

# Third party library imports
import pysam
from tqdm import tqdm
import numpy as np

## Local imports
from pyBioTools.common import *

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~filter_reads~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def Filter (
    input_fn:[str],
    output_fn:str,
    min_len:int=None,
    min_qual:float=None,
    remove_duplicates:bool=False,
    qual_offset:int=33,
    verbose:bool=False,
    quiet:bool=False,
    progress:bool=False,
    **kwargs):
    """
    Filter fastq reads based on their length, mean quality and the presence of duplicates. Can also be used to concatenate reads from multiple files in a single one.
    * input_fn
        Fastq file path or directory containing fastq files or list of files, or regex or list of regex. It is quite flexible.
    * output_fn
        Destination fastq file. Automatically gzipped if the .gz extension is found
    * min_len
        Minimal reads length
    * min_qual
        Minimal mean read PHRED quality
    * remove_duplicates
        If true duplicated reads with the same read id are discarded
    * qual_offset
        Quality scoring system off set. Nowadays pretty much everyone uses +33
    """

    # Define logger
    logger = get_logger (name="Fastq_Filter", verbose=verbose, quiet=quiet)
    logger.warning("Running Fastq Filter")

    # Define output mode
    open_fun, open_mode = (gzip.open, "wt") if is_gziped(output_fn) else (open, "w")

    # Parse reads
    logger.info("Parsing reads")
    c = Counter()
    read_ids = set()

    try:
        with tqdm (desc="Reads processed ", unit=" reads", disable=not progress) as p:
            with open_fun(output_fn, open_mode) as fp_out:
                for fn in super_iglob (input_fn, regex_list=["*.fastq","*.fq","*.fastq.gz","*.fq.gz"]):
                    c["source files"]+=1
                    with pysam.FastxFile(fn) as fp_in:
                        logger.debug("Reading file {}".format(fn))
                        for read in fp_in:

                            if min_len and len(read.sequence) < min_len:
                                c["short_reads"]+=1

                            elif min_qual and np.mean(read.get_quality_array(qual_offset)) < min_qual:
                                c["low_qual_reads"]+=1

                            elif remove_duplicates and read.name in read_ids:
                                c["duplicate_reads"]+=1

                            else:
                                # Write valid read to dest file
                                fp_out.write("{}\n".format(read))
                                c["valid_reads"]+=1

                            # Update counters
                            read_ids.add(read.name)
                            c["total_reads"]+=1
                            p.update()

                        logger.debug("End of file {}".format(fn))

    except (StopIteration, KeyboardInterrupt):
        pass

    # Print read count summary
    logger.info("Read counts summary")
    log_dict(c, logger.info)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Helper writer class~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Writer ():
    def __init__(self, fn, verbose=False, quiet=False):
        """
        Very basic fastq writting class
        * fn
            Output fastq path. Automatically gzipped if.gz
        """
        # Init logger and counter
        self.log = get_logger (name="Fastq_Writer", verbose=verbose, quiet=quiet)
        self.n_seq = 0

        # Define output mode
        open_fun, open_mode = (gzip.open, "wt") if is_gziped(fn) else (open, "w")

        # Open input file in writing mode
        self.fn = fn
        self.log.debug("Opening file {} in writing mode".format(self.fn))
        self.fp = open_fun(fn, open_mode)

    #~~~~~~~~~~~~~~MAGIC METHODS~~~~~~~~~~~~~~#
    def __repr__ (self):
        return "Output target: {}\n\tSequences writen: {}".format(self.fn, self.n_seq)

    def __enter__ (self):
        return self

    def __exit__(self, exception_type, exception_val, trace):
        self.close()

    #~~~~~~~~~~~~~~PUBLIC METHODS~~~~~~~~~~~~~~#
    def close (self):
        try:
            self.log.debug ("Closing file:{}".format(self.fn))
            self.fp.close()
        except Exception as E:
            print (E)
        finally:
            self.log.debug("Sequences writen: {}".format(self.n_seq))

    def write (self, read_id, seq, qual):
        """Add new sequence to fastq record"""
        if not read_id.startswith("@"):
            read_id = "@"+read_id
        # Write sequence
        self.fp.write("{}\n{}\n+\n{}\n".format(read_id, seq, qual))
        self.n_seq+=1
