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
def filter_reads (
    src_fn:str,
    dest_fn:str,
    min_len:int=None,
    min_qual:float=None,
    remove_duplicates:bool=False,
    qual_offset:int=33,
    **kwargs):
    """
    Filter fastq reads based on their length, mean quality and the presence of duplicates. Can also be used to concatenate reads from multiple files in a single one.
    * src_fn
        Fastq file path or directory containing fastq files or list of files, or regex or list of regex. It is quite flexible.
    * dest
        Destination fastq file. Automatically gzipped if the .gz extension is found
    * min_len
        Minimal reads length
    * min_qual
        Minimal mean read PHRED quality
    * remove_duplicates
        If true duplicated reads with the same read id are discarded
    * qual_offset
        Quality scoring system off set. Nowadays pretty much everyone uses +33
    * kwargs
        Allow to pass extra options such as verbose and quiet
    """

    # Define logger
    logger = set_logger (verbose=kwargs.get("verbose", False), quiet=kwargs.get("quiet", False))

    # Define source if directory given
    if os.path.isdir(src_fn):
        src_fn = [os.path.join(src_fn, "*.fastq"), os.path.join(src_fn, "*.fq")]

    # Define output mode
    if is_gziped(dest_fn):
        open_fun, open_mode = gzip.open, "wt"
    else:
        open_fun, open_mode = open, "w"

    # Parse reads
    logger.info("Parsing reads")
    c = Counter()
    read_ids = set()

    try:
        with tqdm (desc="Reads processed ", unit=" reads") as p:
            with open_fun(dest_fn, open_mode) as fp_out:
                for fn in super_iglob (src_fn):
                    c["source files"]+=1
                    with pysam.FastxFile(fn) as fp_in:
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

    except (StopIteration, KeyboardInterrupt):
        pass

    # Print read count summary
    logger.debug("\nRead counts summary")
    logger.debug(dict_to_str(c, ntab=1))
