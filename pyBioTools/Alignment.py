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
from pyBioTools import Fastq

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Reads_index~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def Reads_index (
    input_fn:str,
    skip_unmapped:bool=False,
    skip_secondary:bool=False,
    skip_supplementary:bool=False,
    verbose=False,
    quiet=False,
    progress=False,
    **kwargs):
    """
    Index reads found in a coordinated sorted bam file by read_id.
    The created index file can be used to randon access the alignment file per read_id
    * input_fn
        Path to the bam file to index
    * skip_unmapped
        Filter out unmapped reads
    * skip_secondary
        Filter out secondary alignment
    * skip_supplementary
        Filter out supplementary alignment
    """
    # Define logger
    logger = set_logger (verbose=verbose, quiet=quiet)

    # Check bam
    logger.info("Checking Bam file")
    _check_bam (input_fn)

    # Parse reads
    logger.info("Parsing reads")
    c = defaultdict (Counter)
    idx_fn = input_fn+".idx.gz"
    with pysam.AlignmentFile(input_fn) as bam, gzip.open(idx_fn, "wt") as idx, tqdm(unit=" Reads", disable=not progress) as pbar:
        try:
            while True:
                # Save pointer and read_id
                p = bam.tell()
                read = next(bam)
                pbar.update()

                valid_read, read_status = _eval_read(
                    read=read,
                    skip_unmapped = skip_unmapped,
                    skip_secondary = skip_secondary,
                    skip_supplementary = skip_supplementary)

                if valid_read:
                    c["Reads retained"][read_status]+=1
                    c["Reads retained"]["total"]+=1
                    idx.write("{}\t{}\n".format(read.query_name, p))
                else:
                    c["Reads discarded"][read_status]+=1
                    c["Reads discarded"]["total"]+=1

        except (StopIteration, KeyboardInterrupt):
            pass

    # Print read count summary
    logger.debug("\nRead counts summary")
    logger.debug(dict_to_str(c, ntab=1))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Reads_sample~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def Reads_sample (
    input_fn:str,
    output_folder:str="./",
    output_prefix:str="out",
    n_reads:int=1000,
    n_samples:int=1,
    rand_seed:int=42,
    verbose=False,
    quiet=False,
    progress=False,
    **kwargs):
    """
    Randomly sample `n_reads` reads from a bam file and write downsampled files in `n_samples` bam files.
    If the input bam file is not indexed by read_id `index_reads` is automatically called.
    * input_fn
        Path to the indexed bam file
    * output_folder
        Path to a folder where to write sample files
    * output_prefix
        Path to a folder where to write sample files
    * n_reads
        Number of randomly selected reads in each sample
    * n_samples
        Number of samples to generate files for
    * rand_seed
        Seed to use for the pseudo randon generator. For non deterministic behaviour set to 0
    """

    # Define logger
    logger = set_logger (verbose=verbose, quiet=quiet)

    # Checking bam and index reads if needed
    logger.info("Checking Bam and index file")
    _check_bam (input_fn)
    idx_fn = input_fn+".idx.gz"
    if not file_readable(idx_fn):
        logger.info("Read index not found. Creating one")
        index_reads (input_fn, **kwargs)

    # Create the output folder if it does exist yet
    mkdir(output_folder, exist_ok=True)

    # Load index file
    logger.info("Load index")
    l = []

    with gzip.open(idx_fn, "rt") as idx:
        for line in tqdm (idx, desc="\tIndex", disable=not progress):
            l.append (int(line.rstrip().split("\t")[1]))

    # Get random reads
    logger.info("Write sample reads")

    with pysam.AlignmentFile(input_fn) as bam_in:

        for sample_num in range(1, n_samples+1):
            sample_id = "{:0{}}".format(sample_num, len(str(n_samples)))

            # Generate output file name and open file for writing
            output_fn = os.path.join(output_folder, "{}_{}.bam".format(output_prefix, sample_id))
            with pysam.AlignmentFile(output_fn, "wb", header=bam_in.header) as bam_out:

                # Set random seed for sample and draw random reads
                if rand_seed:
                    random.seed(rand_seed+sample_num)
                p_list = sorted(random.sample(l, n_reads))

                for p in tqdm(p_list, desc="\tSample {}".format(sample_id), unit=" Reads", disable=not progress):
                    bam_in.seek(p)
                    read = next(bam_in)
                    bam_out.write(read)

            logger.info("\tIndexing output bam file")
            pysam.index(output_fn)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Filter~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def Filter (
    input_fn:str,
    output_fn:str,
    skip_unmapped:bool=False,
    skip_secondary:bool=False,
    skip_supplementary:bool=False,
    orientation:str=".",
    min_read_len:int=0,
    min_align_len:int=0,
    min_mapq:int=0,
    min_freq_identity:float=0,
    select_ref:[str]=[],
    exclude_ref:[str]=[],
    verbose=False,
    quiet=False,
    progress=False,
    **kwargs):
    """
    * input_fn
        Path to the bam file to filter
    * output_fn
        Path to the write filtered bam file
    * skip_unmapped
        Filter out unmapped reads
    * skip_secondary
        Filter out secondary alignment
    * skip_supplementary
        Filter out supplementary alignment
    * orientation
        Orientation of alignment on reference genome {"+","-" ,"."}
    * min_read_len
        Minimal query read length (basecalled length)
    * min_align_len
        Minimal query alignment length on reference
    * min_mapq
        Minimal mapping quality score (mapq)
    * min_freq_identity
        Minimal frequency of alignment identity [0 to 1]
    * select_ref
        List of references on which the reads have to be mapped.
    * exclude_ref
        List of references on which the reads should not be mapped.
    """
    # Define logger
    logger = set_logger (verbose=verbose, quiet=quiet)

    # Check bam
    logger.info("Checking input bam file")
    _check_bam (input_fn)

    # Get random reads
    logger.info("Parsing reads")
    c = defaultdict (Counter)
    with pysam.AlignmentFile(input_fn) as bam_in:
        with pysam.AlignmentFile(output_fn, "wb", header=bam_in.header) as bam_out:
            for read in tqdm(bam_in, desc="\t", unit=" Reads", disable=not progress):
                valid_read, read_status = _eval_read(
                    read = read,
                    skip_unmapped = skip_unmapped,
                    skip_secondary = skip_secondary,
                    skip_supplementary = skip_supplementary,
                    orientation = orientation,
                    min_read_len = min_read_len,
                    min_align_len = min_align_len,
                    min_mapq = min_mapq,
                    min_freq_identity = min_freq_identity,
                    select_ref = select_ref,
                    exclude_ref = exclude_ref)

                if valid_read:
                    c["Reads retained"][read_status]+=1
                    c["Reads retained"]["total"]+=1
                    bam_out.write(read)
                else:
                    c["Reads discarded"][read_status]+=1
                    c["Reads discarded"]["total"]+=1

    # Check Bam readability and index
    logger.info("Indexing output bam file")
    pysam.index(output_fn)

    # Print read count summary
    logger.debug("\nRead counts summary")
    logger.debug(dict_to_str(c, ntab=1))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~To_fastq~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def To_fastq (
    input_fn:[str],
    output_r1_fn:str,
    output_r2_fn:str=None,
    ignore_paired_end:bool=False,
    verbose=False,
    quiet=False,
    progress=False,
    **kwargs):
    """
    Dump reads from an alignment file or set of alignment file(s) to a fastq or pair of fastq file(s).
    Only the primary alignment are kept and paired_end reads are assumed to be interleaved.
    Compatible with unmapped or unaligned alignment files as well as files without header.
    * input_fn
        Path (or list of paths) to input BAM/CRAM/SAM file(s)
    * output_r1_fn
        Path to an output fastq file (for Read1 in paired_end mode of output_r2_fn is provided). Automatically gzipped if the .gz extension is found
    * output_r2_fn
        Optional Path to an output fastq file. Automatically gzipped if the .gz extension is found
    * ignore_paired_end
        Ignore paired_end information and output everything in a single file.
    """
    # Define logger
    logger = set_logger (verbose=verbose, quiet=quiet)
    counter = Counter()

    try:
        # Open output fastq file writers
        fastq_r1_fp = Fastq.Writer(output_r1_fn, verbose=verbose)
        fastq_r2_fp = Fastq.Writer(output_r2_fn, verbose=verbose) if output_r2_fn else None

        logger.info("Parsing reads")
        input_fn_list = input_fn if isinstance(input_fn, (list, tuple)) else [input_fn]
        for input_fn in input_fn_list:
            logger.info("Reading input file {}".format(input_fn))
            try:
                with pysam.AlignmentFile(input_fn, check_sq=False) as al_fp, tqdm(desc="\tReading", unit=" Reads", disable=not progress) as pbar:
                    # Can not use a for loop in no header mode
                    while True:
                        r1 = _next_valid(al_fp)
                        # Paired-end mode
                        if ignore_paired_end or r1.is_paired:
                            r2 = _next_valid(al_fp)
                            if r1.query_name != r2.query_name or not r1.is_read1 or not r2.is_read2:
                                raise pyBioToolsError ("Reads are paired in dataset but reads do not seem to be interleaved")
                            # If r2 file provided
                            if output_r2_fn:
                                fastq_r1_fp.write(r1.query_name, r1.query_sequence, r1.qual)
                                fastq_r2_fp.write(r2.query_name, r2.query_sequence, r2.qual)
                                if progress: pbar.update()
                            # else write interleaved fastq file
                            else:
                                fastq_r1_fp.write(r1.query_name+"_1", r1.query_sequence, r1.qual)
                                fastq_r1_fp.write(r2.query_name+"_2", r2.query_sequence, r2.qual)
                                if progress: pbar.update()

                        # Single-end mode
                        else:
                            fastq_r1_fp.write(r1.query_name, r1.query_sequence, r1.qual)
                            if progress: pbar.update()

            except StopIteration:
                logger.debug("\tReached end of input file {}".format(input_fn))

    # Close all output files
    finally:
        for fp in (fastq_r1_fp, fastq_r2_fp):
            try:
                fp.close()
            except:
                pass

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PRIVATE FUNCTION~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def _next_valid (fp):
    """Auto skip secondary and supplementary alignments"""
    while True:
        r = next(fp)
        if not r.is_secondary and not r.is_supplementary:
            return r

def _check_bam (fn):
    """Check bam file readability, index and sorting"""
    if not file_readable(fn):
        raise pyBioToolsError ("Cannot read input bam file")

    with pysam.AlignmentFile(fn) as bam:
        try:
            if bam.header["HD"]["SO"] != "coordinate":
                raise pyBioToolsError ("Bam file not sorted by coordinates")
        except KeyError:
            raise pyBioToolsError ("Bam header missing '@HD SO' section")

        # Index bam file with pysam is needed
        if not bam.has_index():
             pysam.index(fn)

def _eval_read (
    read,
    skip_unmapped = False,
    skip_secondary = False,
    skip_supplementary = False,
    orientation = ".",
    min_read_len = 0,
    min_align_len = 0,
    min_mapq = 0,
    min_freq_identity = 0,
    select_ref = [],
    exclude_ref = []):
    """"""
    # Get unmapped out of the way first
    if read.is_unmapped:
        if skip_unmapped:
            return (False, "unmapped")
        else:
            return (True, "unmapped")

    # Define is read is primary, secondary or supplementary
    if read.is_secondary:
        read_status = "secondary"
        if skip_secondary:
            return (False, read_status)
    elif read.is_supplementary:
        read_status = "supplementary"
        if skip_supplementary:
            return (False, read_status)
    else:
        read_status = "primary"

    # Filter based on orientation
    if (orientation=="+" and read.is_reverse) or (orientation=="-" and not read.is_reverse):
        return (False, "wrong_orientation")
    if min_read_len and read.infer_query_length() < min_read_len :
        return (False, "short_read")
    if min_align_len and read.query_alignment_length < min_align_len :
        return (False, "short_alignment")
    if min_mapq and read.mapping_quality < min_mapq :
        return (False, "low_mapping_quality")
    if min_freq_identity and read.has_tag("NM"):
        qlen = read.query_alignment_length
        edit_dist = read.get_tag("NM")
        if (qlen-edit_dist)/qlen < min_freq_identity:
            return (False, "low_identity")
    if select_ref and read.reference_name not in select_ref :
        return (False, "invalid_reference")
    if exclude_ref and read.reference_name in exclude_ref :
        return (False, "invalid_reference")

    # Return valid code if read went through all the filters
    return (True, read_status)
