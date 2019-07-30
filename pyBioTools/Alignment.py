# -*- coding: utf-8 -*-

# Standard library imports
from collections import *
import gzip
import random

# Third party library imports
import pysam
from tqdm import tqdm

## Local imports
from pyBioTools.common import *

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~INDEX READS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def index_reads (
    bam_fn:str,
    min_mapq:int=0,
    keep_unmapped:bool=False,
    keep_secondary:bool=False,
    keep_supplementary:bool=False,
    **kwargs):
    """
    Index reads found in a coordinated sorted bam file by read_id.
    The created index file can be used to randon access the alignment file per read_id
    * bam_fn
        Path to the bam file to index
    * min_mapq
        Minimal mapq quality of a read to be included in the index
    * keep_unmapped
        Unmapped reads are included in the index
    * keep_secondary
        Secondary alignment are included in the index
    * keep_supplementary
        Supplementary alignment are included in the index
    * kwargs
        Allow to pass extra options such as verbose and quiet
    """
    # Define logger
    logger = set_logger (verbose=kwargs.get("verbose", False), quiet=kwargs.get("quiet", False))

    # Check bam
    logger.info("Checking Bam file")
    check_bam (bam_fn)

    # Parse reads
    logger.info("Parsing reads")
    c = defaultdict (Counter)
    idx_fn = bam_fn+".idx.gz"
    with pysam.AlignmentFile(bam_fn) as bam, gzip.open(idx_fn, "wt") as idx, tqdm(unit=" Reads", disable=logger.level>=30) as pbar:
        try:
            while True:
                # Save pointer and read_id
                p = bam.tell()
                read = next(bam)
                pbar.update()

                # Filter reads
                if read.is_unmapped:
                    if not keep_unmapped:
                        c["skipped"]["unmapped"]+=1
                        continue
                    else:
                        c["retained"]["unmapped"]+=1
                elif read.is_secondary:
                    if not keep_secondary:
                        c["skipped"]["secondary"]+=1
                        continue
                    else:
                        c["retained"]["secondary"]+=1
                elif read.is_supplementary:
                    if not keep_supplementary:
                        c["skipped"]["supplementary"]+=1
                        continue
                    else:
                        c["retained"]["supplementary"]+=1
                elif read.mapq < min_mapq:
                    c["skipped"]["low_mapq"] +=1
                    continue
                else:
                    c["retained"]["primary"] +=1

                # Write read if passed filter
                idx.write("{}\t{}\n".format(read.query_name, p))

        except (StopIteration, KeyboardInterrupt):
            pass

    # Print read count summary
    logger.debug("\nRead counts summary")
    logger.debug(dict_to_str(c, ntab=1))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~INDEX READS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def sample_reads (
    bam_fn:str,
    out_folder:str="./",
    out_prefix:str="out",
    n_reads:int=1000,
    n_samples:int=1,
    rand_seed:int=42,
    **kwargs):
    """
    Randomly sample `n_reads` reads from a bam file and write downsampled files in `n_samples` bam files.
    If the input bam file is not indexed by read_id `index_reads` is automatically called.
    * bam_fn
        Path to the indexed bam file
    * out_folder
        Path to a folder where to write sample files
    * out_prefix
        Path to a folder where to write sample files
    * n_reads
        Number of randomly selected reads in each sample
    * n_samples
        Number of samples to generate files for
    * rand_seed
        Seed to use for the pseudo randon generator. For non deterministic behaviour set to 0
    * kwargs
        Allow to pass extra options such as verbose and quiet
    """

    # Define logger
    logger = set_logger (verbose=kwargs.get("verbose", False), quiet=kwargs.get("quiet", False))

    # Checking bam and index reads if needed
    check_bam (bam_fn)
    idx_fn = bam_fn+".idx.gz"
    if not file_readable(idx_fn):
        logger.info("Read index not found. Creating one")
        index_reads (bam_fn, **kwargs)

    # Create the output folder if it does exist yet
    mkdir(out_folder, exist_ok=True)

    # Load index file
    logger.info("Load index")
    l = []

    with gzip.open(idx_fn, "rt") as idx:
        for line in tqdm (idx, desc="\tIndex", disable=logger.level>=30):
            l.append (int(line.rstrip().split("\t")[1]))

    # Get random reads
    logger.info("Write sample reads")

    with pysam.AlignmentFile(bam_fn) as bam_in:

        for sample_num in range(1, n_samples+1):
            sample_id = "{:0{}}".format(sample_num, len(str(n_samples)))

            # Generate output file name and open file for writing
            bam_out_fn = os.path.join(out_folder, "{}_{}.bam".format(out_prefix, sample_id))
            with pysam.AlignmentFile(bam_out_fn, "wb", header=bam_in.header) as bam_out:

                # Set random seed for sample and draw random reads
                if rand_seed:
                    random.seed(rand_seed+sample_num)
                p_list = sorted(random.sample(l, n_reads))

                for p in tqdm(p_list, desc="\tSample {}".format(sample_id), unit=" Reads", disable=logger.level>=30):
                    bam_in.seek(p)
                    read = next(bam_in)
                    bam_out.write(read)

            logger.info("\tIndex output file")
            pysam.index(bam_out_fn)

def check_bam (bam_fn):
    """Check bam file readability, index and sorting"""
    if not file_readable(bam_fn):
        raise pyBioToolsError ("Cannot read input bam file")

    with pysam.AlignmentFile("./data/sample_1.bam") as bam:
        try:
            if bam.header["HD"]["SO"] != "coordinate":
                raise pyBioToolsError ("Bam file not sorted by coordinates")
        except KeyError:
            raise pyBioToolsError ("Bam header missing '@HD SO' section")

        # Index bam file with pysam is needed
        if not bam.has_index():
             pysam.index(bam_fn)
