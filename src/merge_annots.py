#!/usr/bin/env python

"""
merge_annots.py

Merge custom gene annotations with parsed Prokka output resulting from
parse_gbk.py for downstream use in GenCoFlow.

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
        Merge custom gene annotations with parsed Prokka output resulting from
        parse_gbk.py for downstream use in GenCoFlow.
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
        an ARG or MGE identifer. Must contain the following columns in any
        order, with headers in the first row: "contig_id", "start", "end", 
        "strand".
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



#-------------------------------------------------------------------------------
# main()
#-------------------------------------------------------------------------------

def main(args):
    orfs = parse_gbk(genbank_file = args.input)
    df = orfs_to_df(orfs)
    export_tsv(df, output = args.output)

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    args = parse_arguments()
    main(args)