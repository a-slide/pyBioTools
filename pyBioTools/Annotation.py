# -*- coding: utf-8 -*-

def gff3_parser(self, fp, compression=None, verbose=False, **kwargs):
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
