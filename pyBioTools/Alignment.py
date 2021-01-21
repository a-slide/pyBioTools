# -*- coding: utf-8 -*-

# Standard library imports
from collections import *
import gzip
import random
import math
import tempfile

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
    verbose:bool=False,
    quiet:bool=False,
    progress:bool=False,
    **kwargs):
    """
    Index reads found in a coordinated sorted bam file by read_id.
    The created index file can be used to randon access the alignment file per read_id
    * input_fn
        Path to the bam file to index
    * skip_unmapped
        Do not include unmapped reads in index
    * skip_secondary
        Do not include secondary alignment in index
    * skip_supplementary
        Do not include supplementary alignment in index
    """
    # Define logger
    logger = get_logger (name="Alignment_Reads_index", verbose=verbose, quiet=quiet)
    logger.warning("Running Alignment Reads_index")

    # Check bam
    logger.info("Checking Bam file")
    _check_bam (input_fn)

    # Parse reads
    logger.info("Parsing reads")
    c = defaultdict (Counter)
    idx_fn = input_fn+".idx.gz"
    try:
        with pysam.AlignmentFile(input_fn) as bam, gzip.open(idx_fn, "wt") as idx, tqdm(unit=" Reads", disable=not progress) as pbar:
            idx.write("read_id\tref_id\tpointer\tindex\n")
            try:
                i=0
                while True:

                    # Save pointer and read_id
                    p = bam.tell()
                    read = next(bam)
                    pbar.update()

                    valid_read, read_status = _eval_read(read=read, skip_unmapped = skip_unmapped, skip_secondary = skip_secondary, skip_supplementary = skip_supplementary)

                    if valid_read:
                        c["Reads retained"][read_status]+=1
                        c["Reads retained"]["total"]+=1
                        idx.write("{}\t{}\t{}\t{}\n".format(read.query_name, read.reference_name, p, i))
                    else:
                        c["Reads discarded"][read_status]+=1
                        c["Reads discarded"]["total"]+=1

                    i+=1

            except (StopIteration, KeyboardInterrupt):
                pass

    # Print read count summary
    finally:
        logger.info("Read counts summary")
        log_dict(c, logger.info)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Reads_sample~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def Reads_sample (
    input_fn:str,
    output_folder:str="./",
    output_prefix:str="out",
    n_reads:int=1000,
    n_samples:int=1,
    rand_seed:int=42,
    verbose:bool=False,
    quiet:bool=False,
    progress:bool=False,
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
    logger = get_logger (name="Alignment_Reads_sample", verbose=verbose, quiet=quiet)
    logger.warning("Running Alignment Reads_sample")

    # Checking bam and index reads if needed
    logger.info("Checking Bam and index file")
    _check_bam (input_fn)
    idx_fn = input_fn+".idx.gz"
    if not file_readable(idx_fn):
        logger.info("Read index not found. Creating one")
        Reads_index (input_fn, **kwargs)

    # Create the output folder if it does exist yet
    mkdir(output_folder, exist_ok=True)

    # Load index file
    logger.info("Load index")
    l = []

    with gzip.open(idx_fn, "rt") as idx_fp:
        _= next(idx_fp)
        for line in tqdm (idx_fp, desc="\tIndex", disable=not progress):
            l.append (int(line.rstrip().split("\t")[2]))

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

            logger.info("Indexing output bam file")
            pysam.index(output_fn)

def References_sample (
    input_fn:[str],
    output_fn:str="out.bam",
    selected_reads_fn:str="select_ref.txt",
    frac_reads:float=0.5,
    min_reads_ref:int=30,
    rand_seed:int=42,
    sorting_threads:int=4,
    verbose:bool=False,
    quiet:bool=False,
    progress:bool=False,
    **kwargs):
    """
    Randomly sample reads per references according to a fraction of the reads mapped to this reference for one or
    several files and write selected reads in a new bam file
    * input_fn
        Bam file path or directory containing bam files or list of files, or regex or list of regex.
        It is quite flexible. All files need to be sorted and aligned to the same reference file.
    * output_fn
        Path to the output bam file (sorted and indexed)
    * selected_reads_fn
        Path to the output text file containing all the read id selected
    * frac_reads
        Fraction of reads mapped to sample for each reference
    * min_reads_ref
        Minimal read coverage per file and reference before sampling
    * rand_seed
        Seed to use for the pseudo randon generator. For non deterministic behaviour set to None
    * sorting_threads
        Number of threads to use for bam file sorting
    """

    # Define logger
    logger = get_logger (name="Alignment_Ref_sample", verbose=verbose, quiet=quiet)
    logger.warning("Running Alignment Ref_sample")

    # Set random seed
    random.seed(rand_seed)

    # Create the output folders if they do not exist yet
    for fn in [output_fn, selected_reads_fn]:
        mkdir(os.path.dirname(fn), exist_ok=True)

    # Load index file
    logger.warning("Index files")

    header = None
    d = OrderedDict()
    c =  Counter()
    for fn in super_iglob(input_fn):
        logger.info("Indexing alignment file {}".format(fn))

        # Check Bam file compliance
        _check_bam (fn)

        try:
            with pysam.AlignmentFile(fn) as bam, tqdm(desc="\tReading ", unit=" Reads", disable=not progress) as pbar:
                # Collect header and verify that all files where aligned with the same reference and parameters
                if header is None:
                    header = bam.header
                elif bam.header.references != header.references:
                    raise pyBioToolsError("Bam Header from file {fn} is not consistent with previous files")

                # Parse reads and save index
                while True:
                    loc = bam.tell()
                    read = next(bam)
                    pbar.update()

                    if read.is_unmapped:
                        c["unmapped reads"]+=1
                    elif read.is_secondary:
                        c["secondary reads"]+=1
                    elif read.is_supplementary:
                          c["supplementary reads"]+=1
                    else:
                        c["primary reads"]+=1
                        ref_id = read.reference_name
                        if not ref_id in d:
                            d[ref_id] = OrderedDict()
                        if not fn in d[ref_id]:
                            d[ref_id][fn] = []
                        d[ref_id][fn].append(loc)

        except (StopIteration, KeyboardInterrupt):
            pass

    logger.info("Raw read counts summary")
    log_dict(c, logger.info)

    # Select reads per ref
    logger.warning("Randomly pick reads per references")
    select_loc_d = defaultdict(list)
    c = Counter()
    for ref_id, fn_d in d.items():

        valid_cov = True
        if min_reads_ref is not None:
            for fn, loc_list in fn_d.items():
                n_reads = len(loc_list)
                if n_reads < min_reads_ref:
                    valid_cov = False
                    break

        if valid_cov is True:
            c["valid references"]+=1
            for fn, loc_list in fn_d.items():
                n_reads = len(loc_list)
                c["valid reads"]+=n_reads
                n_samples = math.ceil(len(loc_list)*frac_reads)
                sample_loc_list = random.sample(loc_list, n_samples)
                select_loc_d[fn].extend(sample_loc_list)
                c["valid sampled reads"]+=n_samples

        else:
            c["low coverage references skipped"]+=1

    logger.warning("Sample reads and write to output file")
    with tempfile.NamedTemporaryFile() as temp_bam:
        with pysam.AlignmentFile(temp_bam.name, "wb", header=header) as bam_out, open(selected_reads_fn, "w") as selected_reads_fp:
            for fn, loc_list in select_loc_d.items():
                loc_list = sorted(loc_list)
                logger.info("Writing selected reads for bam file {}".format(fn))
                with pysam.AlignmentFile(fn) as bam_in:
                    for loc in tqdm(loc_list, desc="\tWriting ", unit=" Reads", disable=not progress):
                        bam_in.seek(loc)
                        read = next(bam_in)
                        bam_out.write(read)
                        selected_reads_fp.write("{}\n".format(read.query_name))

        logger.info("Sort BAM File")
        pysam.sort(temp_bam.name, "-o", output_fn, "-@", str(sorting_threads))

    logger.info("Index sorted BAM File")
    pysam.index(output_fn)

    logger.info("Selected read counts summary")
    log_dict(c, logger.info)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Filter~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def Filter (
    input_fn:str,
    output_fn:str,
    selected_reads_fn:str=None,
    skip_unmapped:bool=False,
    skip_secondary:bool=False,
    skip_supplementary:bool=False,
    index_reads:bool=False,
    orientation:str=".",
    min_read_len:int=0,
    min_align_len:int=0,
    min_mapq:int=0,
    min_freq_identity:float=0,
    select_ref_fn:str=None,
    exclude_ref_fn:str=None,
    verbose:bool=False,
    quiet:bool=False,
    progress:bool=False,
    **kwargs):
    """
    * input_fn
        Path to the bam file to filter
    * output_fn
        Path to the write filtered bam file
    * selected_reads_fn
        Optional file where to write ids of selected reads
    * skip_unmapped
        Filter out unmapped reads
    * skip_secondary
        Filter out secondary alignment
    * skip_supplementary
        Filter out supplementary alignment
    * index_reads
        Index bam file with both pysam and pybiotools reads_index
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
    * select_ref_fn
        File containing a list of references on which the reads have to be mapped.
    * exclude_ref_fn
        File containing a list of references on which the reads should not be mapped.
    """
    # Define logger
    logger = get_logger (name="Alignment_Filter", verbose=verbose, quiet=quiet)
    logger.warning("Running Alignment Filter")

    # Check bam
    logger.info("Checking input bam file")
    _check_bam (input_fn)

    # Import selected and excuded references list
    select_ref=[]
    if select_ref_fn:
        logger.info("Loading selected reference file ids")
        with open(select_ref_fn, "r") as fp:
            for l in fp:
                select_ref.append(l.strip())
            logger.debug("Found {} references to select".format(len(select_ref)))
    exclude_ref=[]
    if exclude_ref_fn:
        logger.info("Loading excluded reference file ids")
        with open(exclude_ref_fn, "r") as fp:
            for l in fp:
                exclude_ref.append(l.strip())
            logger.debug("Found {} references to exclude".format(len(exclude_ref)))

    # Make output dir if needed
    logger.debug("Make output directories")
    mkdir(os.path.dirname(output_fn), exist_ok=True)
    if selected_reads_fn:
        mkdir(os.path.dirname(selected_reads_fn), exist_ok=True)

    # Parse reads
    logger.info("Parsing reads")
    c = defaultdict (Counter)
    selected_read_ids = set()
    try:
        with pysam.AlignmentFile(input_fn) as bam_in:
            with pysam.AlignmentFile(output_fn, "wb", header=bam_in.header) as bam_out:
                for read in tqdm(bam_in, desc="\tReading ", unit=" Reads", disable=not progress):
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
                        selected_read_ids.add(read.query_name)
                    else:
                        c["Reads discarded"][read_status]+=1
                        c["Reads discarded"]["total"]+=1

        if index_reads:
            logger.info("Indexing output bam file with pysam")
            pysam.index(output_fn)

            logger.info("Indexing output bam file with pyBiotools")
            Reads_index(output_fn)

        if selected_reads_fn:
            logger.info("Write selected read_id list")
            with open (selected_reads_fn, "w") as fp:
                for read_id in sorted(selected_read_ids):
                    fp.write("{}\n".format(read_id))

    # Print read count summary
    finally:
        logger.info("Read counts summary")
        log_dict(c, logger.info)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~To_fastq~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def To_fastq (
    input_fn:[str],
    output_r1_fn:str,
    output_r2_fn:str=None,
    ignore_paired_end:bool=False,
    verbose:bool=False,
    quiet:bool=False,
    progress:bool=False,
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
    logger = get_logger (name="Alignment_To_fastq", verbose=verbose, quiet=quiet)
    logger.warning("Running Alignment To_fastq")

    # Open output fastq file writers
    fastq_r1_fp = Fastq.Writer(output_r1_fn, verbose=verbose, quiet=quiet)
    fastq_r2_fp = Fastq.Writer(output_r2_fn, verbose=verbose, quiet=quiet) if output_r2_fn else None

    try:
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
                logger.debug("Reached end of input file {}".format(input_fn))

    # Close all output files
    finally:
        for fp in (fastq_r1_fp, fastq_r2_fp):
            try:
                fp.close()
            except:
                pass

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Split~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def Split (
    input_fn:str,
    output_dir:str="",
    n_files:int=10,
    output_fn_list:[str]=[],
    index:bool=False,
    verbose:bool=False,
    quiet:bool=False,
    progress:bool=False,
    **kwargs):
    """
    Split reads in a bam file in N files. The input bam file has to be sorted by coordinates and indexed.
    The last file can contain a few extra reads.
    * input_fn
        Path to the bam file to filter
    * output_dir
        Path to the directory where to write split bam files.
        Files generated have the same basename as the source file and are suffixed with numbers starting from 0
    * n_files
        Number of file to split the original file into
    * output_fn_list
        As an alternative to output_dir and n_files one can instead give a list of output files.
        Reads will be automatically split between the files in the same order as given
    * index
        Index output BAM files
    """
    # Define logger
    logger = get_logger (name="Alignment_Split", verbose=verbose, quiet=quiet)
    logger.warning("Running Alignment Split")

    # Check bam
    logger.info("Checking input bam file")
    _check_bam (input_fn)

    # Generate list of output files if not given
    if not output_fn_list:
        fn_basename = os.path.basename(input_fn).rpartition(".")[0]
        for chunk in range(n_files):
            output_fn_list.append(os.path.join(output_dir, "{}_{}.bam".format(fn_basename, chunk)))
    else:
        n_files= len(output_fn_list)

    logger.debug("List of output files to generate:")
    log_list(output_fn_list, logger.debug)

    logger.info("Parsing reads")
    c = Counter()
    try:
        with pysam.AlignmentFile(input_fn) as bam_in:
            logger.debug("Counting reads")
            total = bam_in.mapped+bam_in.unmapped
            n_reads_per_chunk = total // n_files
            c["Reads from index"] = total
            c["Reads per file"] = n_reads_per_chunk

            with tqdm(total=total, desc="\tReading", unit=" Reads", disable=not progress) as pbar:
                for chunk, output_fn in enumerate(output_fn_list):
                    logger.debug("Open ouput file '{}'".format(output_fn))
                    mkbasedir (output_fn, exist_ok=True)
                    with pysam.AlignmentFile(output_fn, "wb", template=bam_in) as bam_out:

                        n_reads = 0
                        while True:
                            try:
                                read = next(bam_in)
                            except StopIteration:
                                logger.debug("Reached end of input file")
                                break

                            bam_out.write(read)
                            pbar.update()
                            c["Reads writen"]+=1
                            n_reads+=1

                            if n_reads == n_reads_per_chunk and chunk < n_files-1:
                                break

                    # Close file and index
                    logger.debug("Close output file '{}'".format(output_fn))
                    logger.debug("Reads written: {:,}".format(n_reads))
                    if index:
                        logger.debug("index output file '{}'".format(output_fn))
                        pysam.index (output_fn)

    # Print read count summary
    finally:
        logger.info("Read counts summary")
        log_dict(c, logger.info)

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
