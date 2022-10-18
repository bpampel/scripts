#!/usr/bin/env python3
"""
Merge columns of many files
"""

import argparse
import numpy as np

from helpers import number_format
from helpers import plumed_header


def parse_args():
    """Get cli args"""
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="+",
                        help="Path of the files to be merged")
    parser.add_argument("-o", "--outfile",
                        help="Path of the files to be merged")
    parser.add_argument("-c", "--column", type=int,
                        help="Column to be joined")
    args = parser.parse_args()
    return args


args = parse_args()
args.column -= 1 # offset for python numbering

number_string_from_file = number_format.get_string_from_file(args.files[0], args.column)
fmt_output = number_format.NumberFmt(number_string_from_file)

# parse header and replace fields with the same one
header = plumed_header.PlumedHeader()
header.parse_file(args.files[0])
header.fields = [f"{header.fields[args.column]}.{i}" for i in range(len(args.files))]

data = [np.genfromtxt(f)[:,args.column] for f in args.files]
np.savetxt(args.outfile, np.c_[data].T, fmt=fmt_output.get(), header=str(header))
