# -*- coding: utf-8 -*-

# Standard library imports
from collections import *
import gzip
import csv
import os

# Third party library imports
from tqdm import tqdm

## Local imports
from pyBioTools.common import *

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~filter_reads~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

def Merge (
    input_fn:[str],
    output_fn:str,
    old_filename_synthax:bool=False,
    verbose:bool=False,
    quiet:bool=False,
    progress:bool=False,
    **kwargs):
    """
    * input_fn
        Sequencing summary file path or directory containing Sequencing summary file  or list of files,
        or regex or list of regex. It is quite flexible. Files can also be gzipped
    * output_fn
        Destination sequencing summary file. Automatically gzipped if the .gz extension is found
    * old_filename_synthax
        Replace the `filename_fast5` field by `filename` as in older versions.
        Useful for nanopolish index compatibility
    """

    # Define logger
    logger = get_logger (name="Seqsum_Merge", verbose=verbose, quiet=quiet)
    logger.warning("Running Seqsum Merge")

    # Parse reads
    logger.info("Parsing reads")
    c = Counter()
    main_header = None
    read_id_list = set()
    duplicated_ids=0
    try:
        # Define output mode
        open_fun_w, open_mode_w = (gzip.open, "wt") if is_gziped(output_fn) else (open, "w")
        with open_fun_w(output_fn, open_mode_w) as fp_out:
            for fn_in in super_iglob (input_fn, regex_list=["*sequencing_summary.txt.gz", "*sequencing_summary.txt"]):
                c["Files found"]+=1

                # Define input mode
                open_fun_r, open_mode_r = (gzip.open, "rt") if is_gziped(fn_in) else (open, "r")
                with open_fun_r(fn_in, open_mode_r) as fp_in:
                    logger.debug("Reading file {}".format(fn_in))

                    # First header line case
                    header = next(fp_in).rstrip().split("\t")
                    if old_filename_synthax and "filename_fast5" in header and "filename" not in header:
                        i = header.index("filename_fast5")
                        header[i]="filename"

                    if main_header is None:
                        read_id_index = header.index("read_id")
                        main_header = header
                        fp_out.write("\t".join(header)+"\n")

                    elif header != main_header:
                        logger.error("Header of file `{}` is not consistant".format(fn_in))
                        logger.debug("Skipping file {}".format(fn_in))
                        c["Invalid files"]+=1
                        continue

                    # Rest of the lines
                    for line in tqdm(fp_in, desc="Line processed ", unit=" lines", disable=not progress):
                        read_id = line.split("\t")[read_id_index]
                        if read_id in read_id_list:
                            duplicated_ids+=1
                        else:
                            read_id_list.add(read_id)
                            fp_out.write(line)
                            c["Valid lines"]+=1

                    c["Valid files"]+=1
                    logger.debug("End of file {}".format(fn_in))

    except (StopIteration, KeyboardInterrupt):
        pass

    except:
        logger.error("An error occured. Atempting to delete corrupted output file")
        try:
            os.remove(output_fn)
            logger.error("Output file removed")
        except:
            pass
            logger.error("Could not remove output file")
        raise


    # Print read count summary
    if duplicated_ids:
        logger.error("{} duplicated read ids in input files".format(duplicated_ids))

    logger.info("Read counts summary")
    log_dict(c, logger.info)
