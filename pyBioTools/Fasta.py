# -*- coding: utf-8 -*-

# Standard library imports
import gzip

## Local imports
from pyBioTools.common import *

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~API~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

def Reader (fn):
    """very simple generator API for fasta files, including gzipped file support"""
    open_fun, open_mode = (gzip.open, "rt") if is_gziped(fn) else (open, "r")
    with open_fun(fn, open_mode) as fp:
        name = ""
        seq = []
        for line in fp:
            if line.startswith(">"):
                if name and seq:
                    yield FastaRecord(name, "".join(seq))
                # start new sequence
                name = line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())

        # last sequence
        if name and seq:
            yield FastaRecord(name, "".join(seq))

class FastaRecord ():
    """Simple object containing data for single fasta record"""
    def __init__ (self,name, seq):
        self.short_name = name.partition(" ")[0]
        self.long_name = name
        self.seq = seq

    def __len__ (self):
        return len(self.seq)

    def __repr__ (self):
        if len(self) > 50:
            return "{}: {}...".format(self.short_name, self.seq[0:50])
        else:
            return "{}: {}".format(self.short_name, self.seq)

    def __str__ (self):
        return (">{}\n{}\n".format(self.short_name, self.seq))
