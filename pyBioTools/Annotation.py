# -*- coding: utf-8 -*-

# Strandard library imports
import gzip
from collections import OrderedDict, namedtuple, Counter
import sys

# Third party import
import pandas as pd
import numpy as np
from pycltools.pycltools import *

#~~~~~~~~~~~~~~CLASS~~~~~~~~~~~~~~#

class Annotation():
    """
    Parse data from a file containing genomic annotation in GFF3, GTF or BED format.
    Can return the list of annotations for a given interval
    """

    @ classmethod
    def concat (self, annotation_list, name="Merged Annotations", full=False, verbose=False):
        """
        Concatenate a list of Annotation objects together
        * annotation_list
            list of Annotation objects
        * name: str (default "Merged Annotations")
            Name for the Reference object to be created
        * full: BOOL (default False)
            If True all the columns of the annotation dataframe will be stacked (Could be a mess when mixing many
            different kind of annotations). Else only the "core" fields ("refid","start","end","ID","score","strand")
            will be kept
        """
        if verbose: print("Concatenate Annotations")

        df_list = []
        for a in annotation_list:
            if verbose: print ("\t{}".format(a.name))
            try:
                if full:
                    df_list.append (a.feature_df)
                else:
                    df_list.append (a.feature_df[["refid","start","end","ID","score","strand"]])
            except:
                if verbose: print("Skipping invalid annotation object")

        if not df_list:
            return

        # Concatenate all the dataframe and make a hard copy
        df = pd.concat(df_list).copy()

        # Reorder the columns to have the core ones at the begining
        if full:
            core_cols = ["refid","start","end","ID","score","strand"]
            other_cols = [col for col in df.columns if not col in core_cols]
            df = df [core_cols+other_cols]

        A = Annotation(fp=df, name=name, verbose=verbose)
        return A

    #~~~~~~~~~~~~~~FUNDAMENTAL METHODS~~~~~~~~~~~~~~#

    def __init__ (self,
        fp,
        name = "annot",
        min_len = None,
        max_len = None,
        extend_coord_left = None,
        extend_coord_right = None,
        verbose = False,
        invert_selection = None,
        **kwargs):
        """
        Create an Annotation object from a standard genomic annotation file or a previously pickled file.
        * name: str (default "annot")
            Name for the annotation object to be created
         * fp : str
            A path to a standard genomic file containing features annotations among the following format
            gff3: http://www.ensembl.org/info/website/upload/gff3.html
            gtf: http://www.ensembl.org/info/website/upload/gff.html
            bed:  http://www.ensembl.org/info/website/upload/bed.html
            Alternatively, one can use a python pickle file (.pkl) generated during a previous run.
            The file can eventually be compressed in ‘gzip’ format
        * min_len : int (default None)
            Minimal size (start to end) of a feature to be selected
        * max_len : int (default None)
            Maximal size (start to end) of a feature to be selected
        * extend_coord_left : int (default None)
            Number of nucleotides to extend coordinates to the left of each feature
        * extend_coord_right : int (default None)
            Number of nucleotides to extend coordinates to the right of each feature
        * invert_selection: str or list (default None)
            Invert the selection for any of the kwargs fields listed here. For example:
            Annotation(fp=XXX, invert_selection = "type", type = "exon") will select all non exon features
        * kwargs
            In addition, annotation features can be selected based on a list of values (or a single value) for any
            matching fields name. Field names are not predefined in the function arguments and need be one of the
            standard "main fields" or of one of the "attribute fields" defined in gtf and gff3 files. Non-matching
            field names will be ignored.
            * Main fields
                "refid","start","end","ID","score","strand","type"
            * Example of Attribute fields:
                "gene_id", "gene_type", "gene_name", "transcript_id", "transcript_type", "exon_id" ...
        """
        # Store the name if given
        self.name=name

        # For the class concat method or direct instantiation from a dataframe
        if type(fp) == pd.DataFrame:
            if verbose: print ("Create Annotation object '{}' from an existing features dataframe".format(self.name))
            self.feature_df=fp

        # For normal instantiation
        else:
            if verbose: print ("Create Annotation object '{}' from file '{}'".format(self.name, fp))
            # Verify that the file is readable
            is_readable_file(fp)

            # Find if gziped
            if has_extension (fp, pos=-1, ext=["gz","tgz"]):
                if verbose: print("\tFile is gziped")
                compression="gzip"
                ext_pos=-2
            else:
                if verbose: print("\tFile is not compressed")
                compression=None
                ext_pos=-1

            # Find extension type
            if has_extension (fp, pos=ext_pos, ext="gtf"):
                self.feature_df = self._gtf_parser(fp=fp, compression=compression, verbose=verbose)
            elif has_extension (fp, pos=ext_pos, ext="gff3"):
                self.feature_df = self._gff3_parser(fp=fp, compression=compression, verbose=verbose)
            elif has_extension (fp, pos=ext_pos, ext="bed"):
                self.feature_df = self._bed_parser(fp=fp, compression=compression, verbose=verbose)

            # Else try to import as a pickled file
            else:
                try:
                    self.feature_df = self._pickle_parser(fp=fp, verbose=verbose)
                except Exception as E:
                    raise ValueError("\tCannot open file or the file is not in a valid format")

        # Optional filterig steps
        if kwargs:
            self.select_features (**kwargs, invert_selection=invert_selection, sort=False, reindex=False, metaindex=False, verbose=verbose)
        if min_len or max_len:
            self.select_len (min_len=min_len, max_len=max_len, sort=False, reindex=False, metaindex=False, verbose=verbose)

        # Sort the dataframe and reset index
        if not self.feature_df.empty:
            if verbose: print("\tSorting by coordinates and reset index")
            self.feature_df.dropna(axis=1, how="all", inplace=True)
            self.feature_df.sort_values(by=["refid","start","end"], inplace=True)
            self.feature_df.reset_index(drop=True, inplace=True)

        # Optional extension of coordinates
        if extend_coord_left or extend_coord_right:
            self.extend_coord (left=extend_coord_left, right=extend_coord_right, verbose=verbose)

        # Meta index for fast feature access
        self.metaindex = self._make_metaindex()
        if verbose: print("\tTotal imported features: {:,}".format(self.feature_count))

    def __repr__ (self):
        return ("{} - {} - {} features".format(self.__class__.__name__, self.name, self.feature_count))

    #~~~~~~~~~~~~~~PROPERTY METHODS~~~~~~~~~~~~~~#

    @property
    def feature_count(self):
        """Number of features collected
        """
        return len(self.feature_df)

    @property
    def refid_list(self):
        """List of unique reference sequence ids found
        """
        return set(self.feature_df["refid"])

    #~~~~~~~~~~~~~~MAGIC AND GETTER METHODS~~~~~~~~~~~~~~#

    def __len__ (self):
        return self.feature_count

    def __iter__(self):
        """Iterate over the df
        """
        return self.feature_df.itertuples()

    def __getitem__(self, key):
        """Return a feature by position over the df
        """
        return self.feature_df.iloc[key]

    def head (self, n=5):
        """ return the head of the feature dataframe
        """
        return None if self.feature_df.empty else self.feature_df.head(n)

    def tail (self, n=5):
        """ return the head of the feature dataframe
        """
        return None if self.feature_df.empty else self.feature_df.tail(n)

    def sample (self, n=5):
        """ return the random sample of the feature dataframe sorted by index
        """
        if self.feature_df.empty:
            return None
        else:
            if len(self.feature_df)<n:
                n = len(self.feature_df)
            df = self.feature_df.sample(n)
            df.sort_index(inplace=True)
            return df

    def count_field(self, field_name):
        """Number of unique values found for a particular field
        """
        if field_name in self.feature_df:
            return self.feature_df[field_name].nunique()

    def list_field(self, field_name):
        """List of unique values found for a particular field
        """
        if field_name in self.feature_df:
            return set(self.feature_df[field_name])

    def count_uniq_field(self, field_name):
        """List of unique  values found for a particular field and count per values
        """
        if field_name in self.feature_df:
            return pd.DataFrame(self.feature_df.groupby(field_name).size().sort_values(ascending=False), columns=["count"])

    #~~~~~~~~~~~~~~PUBLIC METHODS~~~~~~~~~~~~~~#

    def overlapping_features (self, refid, start, end, strand=".", min_overlap=0, **kwargs):
        """
        Return features overlapping a given genomic interval
        * refid
            Reference sequence id
        * start
            start coordinates of the interval
        * end
            end coordinates of the interval
        * strand
            Strand of the interval
        * min_overlap
            Minimun coordinate overlap in base
        """
        try:
            df =  self.metaindex[refid][strand].query("end >= {} and start <= {}".format (start+min_overlap, end-min_overlap))
            if df.empty:
                return pd.DataFrame()
            else:
                return df
        except KeyError:
            return pd.DataFrame()

    def merge_coord (self,
        stranded = False,
        merge_gap = 0,
        counter_step = 10000,
        verbose = False,
        **kwargs): ####### Could be updated to use the metaindex
        """
        Merge the overlapping annotation features. Feature IDs are concatenated in a new feature id
        Type and all attribute features are lost after merging. if not stranded, the strand information is lost as well.
        A score containing the number of merged feature is created.
        * stranded
            If True, only the features with the same strandness are merged together. If stranded is false the strand
            of the generated features is set to "." [ DEFAULT False ]
        * merge_gap
            Merge features adjacent by less than specified
        * counter_step : int (default 10000)
            Number of sequences to update the counter
        """
        # List to collect all features
        l =[]
        # Reset strand if stranded merging not required
        if not stranded:
            self.feature_df["strand"] = "."

        if verbose: sys.stdout.write ("Merging overlapping features ")
        # Iterate per chromosome and strand
        n=0
        for (refid, strand), df in self.feature_df.groupby(["refid", "strand"]):

            # Iterate on all features of the refid/strand df
            score = 0
            for index, row in df.iterrows():

                # Counter of features processed
                n+=1
                if verbose and n%counter_step==0:
                    sys.stdout.write(".")

                # First element of a continuous interval
                if score == 0:
                    start = row.start
                    end = row.end
                    ID_list = set([row.ID])
                    score = 1

                # Central element of a continuous interval
                elif row.start <= end+merge_gap:
                    if row.end > end:
                        end = row.end
                    ID_list.add(row.ID)
                    score += 1

                # End of previous interval and begining of a new continuous interval
                else:
                    ID_str = ";".join ([str(i) for i in ID_list])
                    l.append ([refid, start, end, ID_str, score, strand])

                    start = row.start
                    end = row.end
                    ID_list = set([row.ID])
                    score = 1

            # Append last element in the main list
            ID_str = ";".join ([str(i) for i in ID_list])
            l.append([refid, start, end, ID_str, score, strand])

        if verbose: print ("\n\tFeatures upstream merging: {:,}".format(self.feature_count))
        self.feature_df=pd.DataFrame(l, columns=["refid","start","end","ID","score","strand"])
        if verbose: print ("\tFeatures after merging: {:,}".format(self.feature_count))

        if verbose: print("\tUpdate metaindex")
        self.metaindex = self._make_metaindex()

    def subtract(self,
        other_annotation,
        stranded = False,
        counter_step = 10000,
        verbose = False,
        **kwargs):
        """
        Subtract the positions of features from another annotation object.
        Features not overlapped will be kept as they are
        Features completly overlapped will be discarded.
        Features partially ovelapped can be either truncated or eventually split.
        Type and all attribute features are lost after subtraction.
        * other_annotation
            Other annotation object containing feature which positions will be removed if they overlap those in self
        * stranded
            If True, only the overlapping features with the same strandness are subtracted [ DEFAULT False ]
        * counter_step : int (default 10000)
            Number of sequences to update the counter
        """
        if verbose: sys.stdout.write ("Subtracting coordinates overlaping the reference annotation features")

        # Collection data structures
        l = []
        f = namedtuple("feature", ["refid", "start", "end", "ID", "score", "strand"])
        non_overlapped_features = overlapped_features = split_features = bases_removed = 0

        # Iterate over all the feature in the reference annotation
        for n, sf in enumerate (self.feature_df.itertuples(index=False, name="feature")):

            if verbose and n%counter_step==0:
                sys.stdout.write(".")

            # Get features from annotation_sub overlapping the current feature
            if stranded:
                overlap_df = other_annotation.overlapping_features(
                    refid=sf.refid, start=sf.start, end=sf.end, strand=sf.strand)
            else:
                overlap_df = other_annotation.overlapping_features(
                    refid=sf.refid, start=sf.start, end=sf.end)

            # If no overlap just append in the list
            if overlap_df.empty:
                l.append(f(refid=sf.refid, start=sf.start, end=sf.end, ID=sf.ID, score=sf.score, strand=sf.strand))
                non_overlapped_features+=1

            # If overlap it more complicated as we have to take into acount that the feature could be split
            else:
                overlapped_features+=1

                # First, create a Series associating position and presence of overlaping feature or not
                position_array = pd.Series(name="Position", index=np.arange(sf.start,sf.end+1), data=0)

                # Iterate over all the overlapping features
                for of in overlap_df.itertuples(index=False, name="feature"):
                    for i in np.arange(of.start, of.end+1):
                        if i in position_array.index:
                            position_array[i] = 1

                # Second, parse the position array to extract all non overlapping intervals
                found = False
                for pos, val in position_array.items():

                    # If no overlap was found for this position
                    if val == 0:
                        # If first position of a non-overlapping interval
                        if not found:
                            found = True
                            start = pos
                        # If not first position of a non-overlapping interval
                        else:
                            end = pos
                            # Special case for the last position of the reference feature
                            if pos == sf.end:
                                l.append(f(refid=sf.refid, start=start, end=end, ID=sf.ID, score=sf.score, strand=sf.strand))
                                split_features+=1

                    # If overlap add the feature to the list if it is the begining of an overlapping interval
                    if val == 1:
                        bases_removed+=1
                        if found:
                            l.append(f(refid=sf.refid, start=start, end=end, ID=sf.ID, score=sf.score, strand=sf.strand))
                            split_features+=1
                        found = False

        # Convert the list of tuple to a dataframe and update the metaindex
        if verbose: print("\n\tConvert collected features to dataframe")
        self.feature_df = pd.DataFrame(l)
        if verbose: print("\tUpdate metaindex")
        self.metaindex = self._make_metaindex()

        if verbose:
            print ("\tFeatures not overlapped: {:,}".format (non_overlapped_features))
            print ("\tFeatures overlapped {:,}".format (overlapped_features))
            print ("\tTotal Base removed {:,}".format (bases_removed))
            print ("\tFinal features {:,}\n".format (split_features+non_overlapped_features))


    def annotate (self,
        other_annotation,
        stranded = False,
        min_overlap = 0,
        fetch_field_names = [],
        counter_step = 10000,
        remove_no_overlap = False,
        remove_multi_overlap = False,
        multi_overlap_field_names = [],
        verbose = False,
        **kwargs):
        """
        Annotate the positions of self features with ovelapping features from another annotation object.
        With the default parameters features not overlapped will be kept as they are. If multiple features from the
        other annotation object overlap, a non-redundant semicolon separated list is saved. If a field name already
        exists in the self annotation it is overriden by the new one.
        * other_annotation
            An annotation object containing interval coordinates to annotate the reads of all the intervals
            can also be a list of Annotation object that will subsequently be concatenated prior annotation
        * stranded (default False)
            If True, only the overlapping features with the same strandness are subtracted [ DEFAULT False ]
        * min_overlap (default 0)
            Minimun coordinate overlap in base
        * fetch_field_names (default [])
            name(s) of the fields to fetch from the other object for all overlapping features.
            example ["gene_name", "gene_type"]
        * counter_step : int (default 10000)
            Number of sequences to update the counter
        * remove_no_overlap (default False)
            If True non overlaped feature are remove from the results
        * remove_multi_overlap (default False)
            If True features overlapped by multiple annotations are remove from the results
        * multi_overlap_field_names (default [])
            By default same as fetch_field_names. But if given will look for multiple overlap for these fields only
        """
        if type(other_annotation) == list:
            other_annotation = Annotation.concat (other_annotation, full=True, verbose=verbose)

        if verbose: sys.stdout.write ("Annotate features from reference annotation with overlapping features ")

        # Collection data structures
        l =[]
        no_overlap = single_overlap = multi_overlap = 0

        # Cast in list if check if fields exist in the other dataframe
        if type(fetch_field_names) in [int, float, str, bool]:
            fetch_field_names = [fetch_field_names]
        if type(multi_overlap_field_names) in [int, float, str, bool]:
            multi_overlap_field_names = [multi_overlap_field_names]
        for field in fetch_field_names:
            if field not in other_annotation.feature_df:
                print("\nCannot find the field '{}' in the other annotation object".format(field))
                return

        # Iterate over all the feature in the reference annotation
        for n, sf in self.feature_df.iterrows():

            if verbose and n%counter_step==0:
                sys.stdout.write(".")

            # Get features from annotation_sub overlapping the current feature
            strand = sf.strand if stranded else "."
            overlap_df = other_annotation.overlapping_features(
                refid=sf.refid,
                start=sf.start,
                end=sf.end,
                strand=strand,
                min_overlap=min_overlap)

            # No annotation hit
            if overlap_df.empty:
                no_overlap+=1
                if not remove_no_overlap:
                    l.append(sf)

            # 1 annotation hit
            elif len(overlap_df) == 1:
                single_overlap+= 1
                for field in fetch_field_names:
                    sf[field] = overlap_df.iloc[0][field]
                l.append(sf)

            # Multi-annotation but not for given features
            elif multi_overlap_field_names and remove_multi_overlap:
                valid = True
                for field in fetch_field_names:
                    unique_vals = set(overlap_df[field])
                    if field in multi_overlap_field_names and len(unique_vals) > 1:
                        valid = False
                        multi_overlap+=1
                        break
                    else:
                        sf[field] = "|".join([str(v) for v in unique_vals])
                if valid:
                    single_overlap+= 1
                    l.append(sf)

            # Multiple annotation hits
            else:
                multi_overlap+=1
                if not remove_multi_overlap:
                    for field in fetch_field_names:
                        unique_vals = set(overlap_df[field])
                        sf[field] = "|".join([str(v) for v in unique_vals])
                    l.append(sf)

        # Convert the list of tuple to a dataframe and update the metaindex
        if verbose: print("\n\tConvert collected features to dataframe")
        self.feature_df = pd.DataFrame(l)
        if verbose: print("\tUpdate metaindex")
        self.metaindex = self._make_metaindex()

        if verbose:
            print ("\tFeatures without overlap: {:,}".format (no_overlap))
            print ("\tFeatures with single overlap {:,}".format (single_overlap))
            print ("\tFeatures with multiple overlaps {:,}".format (multi_overlap))

    def extend_coord (self, left=None, right=None, verbose=False, **kwargs):
        """
        Extend coordinates left or right. Coordinates lesser that 0 are corrected to 0.
        The object does not know the upper limit of the chromosome size. Consequently it can extend the end coordinate
        to values beyond the chromosome size.
        * left
            Number of nucleotides to extend coordinates to the left of each feature
        * right
            Number of nucleotides to extend coordinates to the right of each feature
        """
        if left:
            if verbose: print("Update left coordinates by {}".format(left))
            self.feature_df["start"]-=left
            # Adjust negative start values
            for index in (self.feature_df[(self.feature_df["start"]<1)].index):
                self.feature_df.loc[index, "start"] = 1

        if right:
            if verbose: print("Update right coordinates by {}".format(left))
            self.feature_df["end"]+=right
            #### Eventually add a right side control with a list of refid length

    def resize_coord (self, size, cap_to_feature_len=True, verbose=False, **kwargs):
        """
        Resize features to a fixed length by extending or reducing coordinates centered on the middle of each feature.
        * size
            size of features after resizing
        * cap_to_feature_size
            If True the final size cannot be longer that the original feature
        """
        if verbose: print ("\tResizing features to a lenght of {}".format(size))
        unchanged = updated = 0

        for index, line in self.feature_df.iterrows():
            feature_len = line["end"]-line["start"]

            # Cap lenght to feature size if needed
            if cap_to_feature_len and size >= feature_len:
                unchanged+=1
            else :
                center = line["start"] + feature_len//2
                self.feature_df.loc [index, "start"] = center-(size//2)
                self.feature_df.loc [index, "end"] = center+(size//2)
                updated+=1

        if verbose:
            print ("\t\t{} Unchanged features".format(unchanged))
            print ("\t\t{} Updated features".format(updated))

    def select_len (self, min_len=None, max_len=None, sort=True, reindex=True, metaindex=True, verbose=False, **kwargs):
        """
        Select features longer or shorter that given values
        * min_len
            Minimal length of features
        * max_len
            Maximal length of features
        * sort
            Sort the dataframe by cordinates after filtering
        * reindex
            reset the integer index of the dataframe by cordinates after filtering
        * metaindex
            update the metaindex after filtering
        """
        if verbose: print ("\tSelecting features based on length")
        l = len(self.feature_df)

        # Filter min len
        if min_len:
            self.feature_df = self.feature_df[((self.feature_df["end"]-self.feature_df["start"]) >= min_len)]
            if verbose: print ("\t\t{} Features removed by minimal length filtering".format(l-len(self.feature_df)))
            l = len(self.feature_df)

        # Filter max len
        if max_len:
            self.feature_df = self.feature_df[((self.feature_df["end"]-self.feature_df["start"]) <= max_len)]
            if verbose: print ("\t\t{} Features removed by maximal length filtering".format(l-len(self.feature_df)))
            l = len(self.feature_df)

        # Cleanup
        if sort:
            if verbose: print("\t\tSorting by coordinates")
            self.feature_df.sort_values(by=["refid","start","end"], inplace=True)
        if reindex:
            if verbose: print("\t\tReset index")
            self.feature_df.reset_index(drop=True, inplace=True)
        if metaindex:
            if verbose: print("\t\tUpdate metaindex")
            self.metaindex = self._make_metaindex(verbose=verbose)

    def select_features (self, invert_selection=None, sort=True, reindex=True, metaindex=True, verbose=False, **kwargs):
        """
        Select annotation features can be selected based on a list of values (or a single value) for any matching
        fields name. Field names are not predefined in the function arguments and need be one of the standard
        "main fields" or of one of the "attribute fields" defined in gtf and gff3 files. Non-matching fields names
        will be ignored.
        * Main fields:
            "refid","start","end","ID","score","strand","type"
        * Example of Attribute fields:
            "gene_id", "gene_type", "gene_name", "transcript_id", "transcript_type", "transcript_name", "exon_id" ...
        * invert_selection: str or list (default None)
            Invert the selection for any of the kwargs fields listed here. For example:
            select_features(invert_selection = "type", type = "exon") will select all non exon features
        * sort
            Sort the dataframe by cordinates after filtering
        * reindex
            reset the integer index of the dataframe by cordinates after filtering
        * metaindex
            update the metaindex after filtering
        """
        l = len(self.feature_df)

        # If atomic type, cast invert_selection in a list
        if invert_selection and type(invert_selection) in [int, float, str, bool]:
            invert_selection = [invert_selection]

        for key, val in kwargs.items():
            if key in self.feature_df and val:

                # If atomic type, cast the value in a list
                if type(val) in [int, float, str, bool]:
                    val = [val]

                # Include or exclude the features matching
                if invert_selection and key in invert_selection:
                    if verbose: print ("\tFiltering out the features with {} in {}".format(key, val))
                    self.feature_df = self.feature_df[(~self.feature_df[key].isin(val))]
                else:
                    if verbose: print ("\tSelecting the features with {} in {}".format(key, val))
                    self.feature_df = self.feature_df[(self.feature_df[key].isin(val))]

                if verbose: print ("\t\t{:,} Features removed".format(l-len(self.feature_df)))
                l = len(self.feature_df)
            else:
                if verbose: print('\t"{}" field was not found in the parsed annotation file'.format(key))

        # Cleanup
        if sort:
            if verbose: print("\t\tSorting by coordinates")
            self.feature_df.sort_values (by=["refid","start","end"], inplace=True)
        if reindex:
            if verbose: print("\t\tReset index")
            self.feature_df.reset_index (drop=True, inplace=True)
        if metaindex:
            if verbose: print("\t\tUpdate metaindex")
            self.metaindex = self._make_metaindex (verbose=verbose)

    def select_longest (self, groupby="gene_id", verbose=False, **kwargs):
        """ Select the longest feature per given field group
        """
        if verbose: print ("Select longest feature per {}".format(groupby))
        assert groupby in self.feature_df, "Groupby field:{} was not found in feature_df".format(groupby)
        l = len(self.feature_df)

        longest_list = []
        for gene_id, sdf in self.feature_df.groupby (groupby):
            max_index = (sdf["end"]-sdf["start"]).idxmax()
            max_row = sdf.loc [max_index].copy()
            longest_list.append (max_row)

        self.feature_df = pd.DataFrame (longest_list)
        self.feature_df.sort_values (by=["refid","start","end"], inplace=True)
        self.feature_df.reset_index (drop=True, inplace=True)
        self.metaindex = self._make_metaindex (verbose=verbose)

        if verbose: print ("\t\t{} Features removed by longest feature filtering".format (l-len(self.feature_df)))

    def to_pickle (self, fp, verbose=False, **kwargs):
        """Store the parsed file in a pickle file for further use. Sugested extension => .pkl
        """
        if verbose: print ("Pickle Annotation '{}' to file '{}'".format (self.name, fp))
        self.feature_df.to_pickle(fp)

    def to_file (self, fp, verbose=False, **kwargs):
        """
        Write the imported and eventually modified annotation features in a file
        * fp
            Path of the file where to output the results. The ouput type will be determined based on the filename.
            The file can be gzip compressed (.gz) or not. Valid file extension are .bed (bed6 format).
            May be extended to other file format in future releases.
        """
        if verbose: print ("Save Annotation '{}' to file '{}'".format (self.name, fp))
        # Find if gziped
        if has_extension (fp, pos=-1, ext=["gz","tgz"]):
            compression=True
            ext_pos=-2
        else:
            compression=False
            ext_pos=-1

        # Find extension type
        if has_extension (fp, pos=ext_pos, ext="bed"):
            self._to_bed(fp=fp, compression=compression, verbose=verbose)

        # else if invalid
        else:
            print ("Invalid file extension. pick one of the following extensions: bed +-gz/tgz")

    def annotated_intervals_metrics (self,
        groupby = None,
        summary_fields = [],
        count_only = []):
        """
        Generate simple count metrics from annotated intervals containing given fields such as
        "gene_id", "gene_name", "end", "start" + "transcript_id" and "exon_id" if available
        """

        # Verify the presence of the requested fields
        if groupby:
            assert groupby in self.feature_df, "Groupby field:{} was not found in feature_df".format(groupby)
        for f in summary_fields:
            assert f in self.feature_df, "Summary field {} was not found in feature_df".format(f)
        for f in count_only:
            assert f in summary_fields, "Count only field {} was not found in summary_fields".format(f)

        # If data is to be grouped by a specific field
        d = OrderedDict()
        if groupby:
            for group_id, group_df in self.feature_df.groupby(groupby):
                d [group_id] = self._fields_decompose (group_df, summary_fields, count_only)
        else:
            d["overall"] = self._fields_decompose (self.feature_df, summary_fields, count_only)

        df =  pd.DataFrame.from_dict (d, orient="index")
        return df

    def reformat_id (self, groupby="gene_id", prefix=None, numbered=True, verbose=False, **kwargs):
        """
        Reformat ID field based on another field (groupby) +- a unique number per field (numbered) +- prefix (prefix)
        """
        assert groupby in self.feature_df, "{} was not found in feature_df".format(groupby)

        c=Counter()
        for line in self.feature_df.itertuples():
            index = line.Index

            # Rename ID
            group = self.feature_df.loc [index, groupby]
            id = str(group)
            if prefix:
                id = "{}_{}".format(prefix, id)
            if numbered:
                c [group] +=1
                id = "{}_{}".format(id, c[group])

            self.feature_df.loc [index, "ID"] = id

    #~~~~~~~~~~~~~~PRIVATE METHODS~~~~~~~~~~~~~~#

    def _bed_parser(self, fp, compression=None, verbose=False, **kwargs):
        """Parse a bed formated file
        """
        if verbose: print("\tUse BED parser to parse annotations")
        # try to import the file as a bed6 in a dataframe
        try:
            col_names = ["refid","start","end","ID","score","strand"]
            df = pd.read_csv(fp, sep="\t", names=col_names, index_col=False, comment="#", compression=compression)

        # else try to import as a bed12
        except IndexError as E:
            col_names = ["refid","start","end","ID","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts"]
            df = pd.read_csv(fp, sep="\t", names=col_names, index_col=False, comment="#", compression=compression)

        df = self._clean_df(df, verbose=verbose)
        return df

    def _gff3_parser(self, fp, compression=None, verbose=False, **kwargs):
        """ Parse a gff3 formated file
        """
        if verbose: print("\tUse GFF3 parser to parse annotations")
        # Define behaviour is compressed or not
        if compression == "gzip":
            open_fun, open_mode = gzip.open, "rt"
        else :
            open_fun, open_mode = open, "r"

        # Parse the file lines by lines and extract values in a list of dict
        all_features = []
        invalid_features=0
        with open_fun(fp, open_mode) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                line = line.strip().split("\t")

                # Try to extract the info from the main fields. If unsuccessful skip the feature
                d = OrderedDict()
                try:
                    d["refid"] = line[0]
                    d["source"] = line[1]
                    d["type"] = line[2]
                    d["start"] = int(line[3])
                    d["end"] = int(line[4])
                    d["score"] = line[5]
                    d["strand"] = line[6]
                    d["frame"] = line[7]
                    attribute_list = line[8].split(";")
                except (IndexError, ValueError) as E:
                    invalid_features+=1
                    continue

                # Try to extract the info from the attribute fields. If unsuccessful skip the attributes and continue
                for attribute in attribute_list:
                    try:
                        k, v = attribute.split("=")
                        d[k] = v
                    except (IndexError, ValueError) as E:
                        pass

                # Check the presence of the gene_id field and append the feature
                if "ID" not in d:
                    d["ID"] = np.nan
                all_features.append(d)

        if verbose and invalid_features:
            print("\t\tFound {:,} invalid lines during file parsing".format(invalid_features))

        # Convert in a dataframe and cleanup
        df = pd.DataFrame(all_features)
        df = self._clean_df(df, verbose=verbose)
        return df

    def _gtf_parser(self, fp, compression=None, verbose=False, **kwargs):
        """ Parse a gtf formated file
        """
        if verbose: print("\tUse GTF parser to parse annotations")
        # Define behaviour is compressed or not
        if compression == "gzip":
            open_fun, open_mode = gzip.open, "rt"
        else :
            open_fun, open_mode = open, "r"

        # Parse the file lines by lines and extract values in a list of dict
        all_features = []
        invalid_features=0
        with open_fun(fp, open_mode) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                line = line.strip().split("\t")

                # Try to extract the info from the main fields. If unsuccessful skip the feature
                d = OrderedDict()
                try:
                    d["refid"] = line[0]
                    d["source"] = line[1]
                    d["type"] = line[2]
                    d["start"] = int(line[3])
                    d["end"] = int(line[4])
                    d["score"] = line[5]
                    d["strand"] = line[6]
                    d["frame"] = line[7]
                    attribute_list = line[8].split(";")
                except (IndexError, ValueError) as E:
                    invalid_features+=1
                    continue

                # Try to extract the info from the attribute fields. If unsuccessful skip the attributes and continue
                for attribute in attribute_list:
                    try:
                        k, v = attribute.split()
                        d[k] = v[1:-1]
                    except (IndexError, ValueError) as E:
                        pass

                # Check the presence of the gene_id field and append the feature
                if "gene_id" in d:
                    d["ID"] = d["gene_id"]
                    del d["gene_id"]
                else:
                    d["ID"] = np.nan
                all_features.append(d)

        if verbose and invalid_features:
            print("\t\tFound {:,} invalid lines during file parsing".format(invalid_features))

        # Convert in a dataframe and cleanup
        df = pd.DataFrame(all_features)
        df = self._clean_df(df, verbose=verbose)
        return df

    def _clean_df(self, df, verbose=False, **kwargs):
        """ Shared cleanup function
        """
        # Drop lines with empty values in fundamental fields
        if verbose: print ("\tFiltering features containing empty values for fundamental fields")
        n_features = len(df)
        df.dropna(subset=["refid","start","end","ID","strand"], inplace=True)
        if verbose: print ("\t\t{:,} features filtered out".format(n_features-len(df)))

        # Cast the start and end field in integer
        if verbose: print("\t\tCast coordinates to integer")
        df[['start', 'end']] = df[['start', 'end']].astype(int)

        # Verify than the dataframe is not empty
        if df.empty:
            raise ValueError("No valid features imported. Is the file valid?")
        return df

    def _pickle_parser (self, fp, verbose=False, **kwargs):
        """Parse a pickle database
        """
        # Import the file in a dataframe
        if verbose: print ("\tTry to load as a pickle file")
        df = pd.read_pickle(fp)
        return df

    def _to_bed (self, fp, compression=False, verbose=False, **kwargs):
        """Output a bed file
        """
        # Select and reorder columns
        df = self.feature_df[["refid","start","end","ID","score","strand"]]
        if compression:
            if verbose: print ("\tOutput a gzip compressed bed file")
            df.to_csv(fp, sep="\t", compression="gzip", index=False, na_rep=".", header=False)

        else:
            if verbose: print ("\tOutput an uncompressed bed file")
            df.to_csv(fp, sep="\t", index=False, na_rep=".", header=False)

    def _make_metaindex (self, verbose=False, **kwargs):
        """ Create a multilevel dict to allow fast access to features by coordinates
        """
        if verbose: print("\tCreate metaindex")
        d = OrderedDict()
        for refid, df in self.feature_df.groupby("refid"):
            d[refid] = { ".": df, "+":df[df["strand"] == "+"], "-":df[df["strand"] == "-"]}
        return d

    def _fields_decompose (self, df, summary_fields, count_only):
        """
        """
        d = OrderedDict()

        # Decompose fields
        for field in summary_fields:

            # Collect non redundant values in a set
            val_set = set ()
            for val_list in df[field]:
                for val in val_list.split("|"):
                    val_set.add (val)

            # Aggregate values
            if field in count_only:
                d["n_{}".format(field)] = len (val_set)
            else:
                d[field] = "|".join(val_set)

        d["n_regions"] = len (df)
        d["n_bases"] = (df["end"]-df["start"]).sum()

        return d
