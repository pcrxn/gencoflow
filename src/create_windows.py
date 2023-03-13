#!/usr/bin/env python

"""
create_windows.py

Create gene windows using custom annotations and parsed Prokka output.

Dependencies: python=3.11.0, pandas=1.5.2
Other package versions may work but are untested.
"""

__author__ = "Liam Brown"
__email__ = "liam.brown@inspection.gc.ca"
__license__ = "Crown Copyright 2023"

import os
import sys
import argparse
import pandas as pd

#-------------------------------------------------------------------------------
# parse_arguments()
#-------------------------------------------------------------------------------

def parse_arguments():
    """
    Parse command-line arguments.
    :returns args: List of parsed arguments.
    """

    parser = argparse.ArgumentParser(
        description = """
        Create gene windows using custom annotations and parsed Prokka output.
        """)

    # Required arguments
    required_args = parser.add_argument_group('Required')
    required_args.add_argument('-p', '--prokka', type = str, required = True,
        help = """
        Path to a TSV file produced by parse_gbk.py.
        """)
    required_args.add_argument('-c', '--custom', type = str, required = True,
        help = """
        Path to a TSV file of custom annotations, such as those produced by
        an ARG identifier. Must contain the following columns in any
        order, with headers in the first row: "contig_id", "start", "end", 
        "strand". Other recognized columns (optional): "seq_id".
        """)
    required_args.add_argument('-o', '--output', type = str, required = True,
        help = """
        Path for a new TSV file of merged annotations.
        """)    

    # If no arguments provided:
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    return args

#-------------------------------------------------------------------------------
# Other functions
#-------------------------------------------------------------------------------

def parse_custom_annot(custom_annot_file):
    custom_annot_df = pd.read_csv(custom_annot_file, sep = '\t')

#-------------------------------------------------------------------------------
# main()
#-------------------------------------------------------------------------------

def main(args):
    parse_custom_annot(custom_annot_file = args.custom)

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    args = parse_arguments()
    main(args)