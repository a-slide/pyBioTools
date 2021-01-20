# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import sys
import os
import inspect
from collections import *
import logging
from glob import iglob
import gzip

# Third party library imports
import pysam
import colorlog

#~~~~~~~~~~~~~~FUNCTIONS~~~~~~~~~~~~~~#

def file_readable (fn, **kwargs):
    """Check if the file is readable"""
    return os.path.isfile (fn) and os.access (fn, os.R_OK)

def dir_writable (fn, **kwargs):
    """Check if the file is readable"""
    if not os.path.isdir(fn):
        fn = os.path.dirname(fn)
    return os.path.dirname(fn) and os.access (fn, os.W_OK)

def is_gziped (fp, **kwargs):
    """ Return True if the file is Gziped else False """
    return fp[-2:].lower() == "gz"

def super_iglob (pathname, recursive=False, regex_list=[]):
    """ Same as iglob but pass multiple path regex instead of one. does not store anything in memory"""
    if type(pathname) == str:
        pathname = [pathname]

    if type(pathname) in [list, tuple, set]:
        for paths in pathname:
            for path in iglob(pathname=paths, recursive=recursive):
                if os.path.isdir(path) and regex_list:
                    for regex in regex_list:
                        regex_paths = os.path.join(path, regex)
                        for regex_path in iglob(pathname=regex_paths, recursive=recursive):
                            yield regex_path
                elif os.path.isfile(path):
                    yield path
    else:
        raise ValueError ("Invalid file type")

def mkdir (fn, exist_ok=False):
    """ Create directory recursivelly. Raise IO error if path exist or if error at creation """
    if fn not in [ "", "..", "."]:
        try:
            os.makedirs (fn, exist_ok=exist_ok)
        except:
            raise pyBioToolsError ("Error creating output folder `{}`".format(fn))

def mkbasedir (fn, exist_ok=False):
    """ Create directory for a given file recursivelly. Raise IO error if path exist or if error at creation """
    dir_fn = os.path.dirname(fn)
    if dir_fn:
        mkdir (dir_fn, exist_ok=True)

def doc_func (func):
    """Parse the function description string"""

    if inspect.isclass(func):
        func = func.__init__

    docstr_list = []
    for l in inspect.getdoc(func).split("\n"):
        l = l.strip()
        if l:
            if l.startswith("*"):
                break
            else:
                docstr_list.append(l)

    return " ".join(docstr_list)

def make_arg_dict (func):
    """Parse the arguments default value, type and doc"""

    # Init method for classes
    if inspect.isclass(func):
        func = func.__init__

    if inspect.isfunction(func) or inspect.ismethod(func):
        # Parse arguments default values and annotations
        d = OrderedDict()
        for name, p in inspect.signature(func).parameters.items():
            if not p.name in ["self","cls"]: # Object stuff. Does not make sense to include in doc
                d[name] = OrderedDict()
                if not name in ["kwargs","args"]: # Include but skip default required and type
                    # Get Annotation
                    if p.annotation != inspect._empty:
                        d[name]["type"] = p.annotation
                    # Get default value if available
                    if p.default == inspect._empty:
                        d[name]["required"] = True
                    else:
                        d[name]["default"] = p.default

        # Parse the docstring in a dict
        docstr_dict = OrderedDict()
        lab=None
        for l in inspect.getdoc(func).split("\n"):
            l = l.strip()
            if l:
                if l.startswith("*"):
                    lab = l[1:].strip()
                    docstr_dict[lab] = []
                elif lab:
                    docstr_dict[lab].append(l)

        # Concatenate and copy doc in main dict
        for name in d.keys():
            if name in docstr_dict:
                d[name]["help"] = " ".join(docstr_dict[name])
        return d

def arg_from_docstr (parser, func, arg_name, short_name=None):
    """Get options corresponding to argument name from docstring and deal with special cases"""

    if short_name:
        arg_names = ["-{}".format(short_name), "--{}".format(arg_name)]
    else:
        arg_names = ["--{}".format(arg_name)]

    arg_dict = make_arg_dict(func)[arg_name]
    if "help" in arg_dict:
        if "default" in arg_dict:
            if arg_dict["default"] == "" or arg_dict["default"] == [] :
                arg_dict["help"] += " (default: None)"
            else:
                arg_dict["help"] += " (default: %(default)s)"
        else:
            arg_dict["help"] += " (required)"

        if "type" in arg_dict:
            arg_dict["help"] += " [%(type)s]"

    # Special case for boolean args
    if arg_dict["type"] == bool:
        if arg_dict["default"] == False:
            arg_dict["action"] = 'store_true'
            del arg_dict["type"]
        elif arg_dict["default"] == True:
            arg_dict["action"] = 'store_false'
            del arg_dict["type"]

    # Special case for lists args
    elif isinstance(arg_dict["type"], list):
        arg_dict["nargs"]='*'
        arg_dict["type"]=arg_dict["type"][0]

    parser.add_argument(*arg_names, **arg_dict)

def jhelp (f:"python function or method"):
    """
    Display a Markdown pretty help message for functions and class methods (default __init__ is a class is passed)
    jhelp also display default values and type annotations if available.
    The docstring synthax should follow the same synthax as the one used for this function
    * f
        Function or method to display the help message for
    """
    # Private import as this is only needed if using jupyter
    from IPython.core.display import display, Markdown

    f_doc = doc_func(f)
    arg_doc = make_arg_dict(f)

    # Signature and function documentation
    s = "**{}** ({})\n\n{}\n\n---\n\n".format(f.__name__, ", ".join(arg_doc.keys()), f_doc)

    # Args doc
    for arg_name, arg_val in arg_doc.items():
        # Arg signature section
        s+= "* **{}**".format(arg_name)
        if "default" in arg_val:
            if arg_val["default"] == "":
                s+=" (default: \"\")".format(arg_val["default"])
            else:
                s+=" (default: {})".format(arg_val["default"])
        if "required" in arg_val:
            s+= " (required)"
        if "type" in arg_val:
            if isinstance(arg_val["type"], type):
                s+= " [{}]".format(arg_val["type"].__name__)
            elif isinstance(arg_val["type"], list):
                s+= " [list({})]".format(arg_val["type"][0].__name__)
            else:
                s+= " [{}]".format(arg_val["type"])
        s+="\n\n"
        # Arg doc section
        if "help" in arg_val:
            s+= "{}\n\n".format(arg_val["help"])

    # Display in Jupyter
    display (Markdown(s))

def get_logger (name=None, verbose=False, quiet=False):
    """Multilevel colored log using colorlog"""

    # Define conditional color formatter
    formatter = colorlog.LevelFormatter(
        fmt = {
            'DEBUG':'%(log_color)s\t[DEBUG]: %(msg)s',
            'INFO': '%(log_color)s\t%(msg)s',
            'WARNING': '%(log_color)s## %(msg)s ##',
            'ERROR': '%(log_color)sERROR: %(msg)s',
            'CRITICAL': '%(log_color)sCRITICAL: %(msg)s'},
        log_colors={
            'DEBUG': 'white',
            'INFO': 'green',
            'WARNING': 'bold_blue',
            'ERROR': 'bold_red',
            'CRITICAL': 'bold_purple'},
        reset=True)

    # Define logger with custom formatter
    logging.basicConfig(format='%(message)s')
    logging.getLogger().handlers[0].setFormatter(formatter)
    logger = logging.getLogger(name)

    # Define logging level depending on verbosity
    if verbose:
        logger.setLevel(logging.DEBUG)
    elif quiet:
        logger.setLevel(logging.WARNING)
    else:
        logger.setLevel(logging.INFO)

    return logger

def log_dict (d, logger, level=1):
    """ log a multilevel dict """

    if isinstance(d, Counter):
        for i, j in d.most_common():
            logger("{}{}: {:,}".format(" "*level, i, j))
    else:
        for i, j in d.items():
            if isinstance(j, dict):
                logger("{}{}".format(" "*level, i, j))
                log_dict(j, logger, level=level+1)
            else:
                logger("{}{}: {}".format(" "*level, i, j))

def log_list (l, logger):
    """ log a list """
    for i in l:
        logger("* {}".format(i))

def head (fp, n=10, ignore_comment_line=False, comment_char="#", max_char_line=300, sep="\t", max_char_col=30, **kwargs):
    """
    Emulate linux head cmd. Handle gziped files and bam files
    * fp
        Path to the file to be parsed. Works with text, gunziped and binary bam/sam files
    * n
        Number of lines to print starting from the begining of the file (Default 10)
    * ignore_comment_line
        Skip initial lines starting with a specific character. Pointless for bam files(Default False)
    * comment_char
        Character or string for ignore_comment_line argument (Default "#")
    * max_char_line
        Maximal number of character to print per line (Default 150)
    """
    line_list = []
    ext = fp.rpartition(".")[2].lower()

    # For bam files
    if ext in ["bam", "sam"]:
        with pysam.AlignmentFile(fp) as f:

            for line_num, read in enumerate(f):
                if line_num >= n:
                    break
                l = read.to_string()
                if sep:
                    line_list.append (l.split(sep))
                else:
                    line_list.append (l)
                line_num+=1

    # Not bam file
    else:
        # For text files
        if ext=="gz":
            open_fun = gzip.open
            open_mode =  "rt"
        else:
            open_fun = open
            open_mode =  "r"

        try:
            with open_fun(fp, open_mode) as fh:
                line_num = 0
                while (line_num < n):
                    l= next(fh).strip()
                    if ignore_comment_line and l.startswith(comment_char):
                        continue
                    if sep:
                        line_list.append (l.split(sep))
                    else:
                        line_list.append (l)
                    line_num+=1

        except StopIteration:
            print ("Only {} lines in the file".format(line_num))

    # Print lines
    if sep:
        try:
            # Find longest elem per col
            col_len_list = [0 for _ in range (len(line_list[0]))]
            for ls in line_list:
                for i in range (len(ls)):
                    len_col = len(ls[i])
                    if len_col > max_char_col:
                        col_len_list[i] = max_char_col
                    elif len_col > col_len_list[i]:
                        col_len_list[i] = len_col

            line_list_tab = []
            for ls in line_list:
                s = ""
                for i in range (len(ls)):
                    len_col = col_len_list[i]
                    len_cur_col = len(ls[i])
                    s += ls[i][0:len_col] + " "*(len_col-len_cur_col)+" "
                line_list_tab.append(s)
            line_list = line_list_tab

        # Fall back to none tabulated display
        except IndexError:
            return head (fp=fp, n=n, ignore_comment_line=ignore_comment_line, comment_char=comment_char, max_char_line=max_char_line, sep=None)

    for l in line_list:
        if max_char_line and len(l) > max_char_line:
            print (l[0:max_char_line]+"...")
        else:
            print (l)
    print()

#~~~~~~~~~~~~~~CUSTOM EXCEPTION AND WARN CLASSES~~~~~~~~~~~~~~#
class pyBioToolsError (Exception):
    """ Basic exception class for NanopolishComp package """
    pass
