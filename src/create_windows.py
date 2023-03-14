#!/usr/bin/env python

"""
create_windows.py

Create gene windows using custom annotations and parsed Prokka output.

Dependencies: python=3.11.0, pandas=1.5.2, numpy=1.24.2
Other package versions may work but are untested.
"""

__author__ = "Liam Brown"
__email__ = "liam.brown@inspection.gc.ca"
__license__ = "Crown Copyright 2023"

import os
import sys
import argparse
import pandas as pd
import numpy as np

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
        an ARG identifier. The first row must contain headers for columns, and
        the following headers must be present:"contig_id", "start", "end", 
        "strand".
        """)
    required_args.add_argument('-o', '--output', type = str, required = True,
        help = """
        Path for a new TSV file of gene windows for downstream analysis.
        """)
    
    # Optional arguments
    optional_args = parser.add_argument_group('Optional')
    optional_args.add_argument('-x', '--annot_field', type = str, 
                               required = False,
        help = """
        Name of a column in -c/--custom to use for annotation.
        """)    
    optional_args.add_argument('-w', '--window_size', type = int, 
                               required = False, default = 3000,
        help = """
        Number of bp flanking the start and end positions of each custom
        annotation to include for analysis.
        Default: 3000
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
    """
    """
    custom_df = pd.read_csv(custom_annot_file, sep = '\t', index_col = False)
    custom_df['source'] = 'custom'
    custom_df['strand'].replace('+', 1, inplace = True)
    custom_df['strand'].replace('-', -1, inplace = True)
    return custom_df

def parse_prokka(parsed_gbk_file):
    """
    """
    prokka_df = pd.read_csv(parsed_gbk_file, sep = '\t', index_col = False)
    prokka_df['source'] = 'prokka'
    return prokka_df

def merge_dfs(custom_df, prokka_df, annot_field):
    """
    """
    if annot_field is None:
        df = pd.concat([
            prokka_df[['contig_id', 'source', 'start', 'end', 'strand']],
            custom_df[['contig_id', 'source', 'start', 'end', 'strand']]
        ])
    # If an annotation field was provided on the cmd line:
    else:
        df = pd.concat([
            prokka_df[['contig_id', 'source', 'start', 'end', 'strand']],
            custom_df[['contig_id', 'source', 'start', 'end', 'strand', annot_field]]
        ])

    # Fix start and end positions if necessary
    # Start should always be less than end, regardless of strand
    df['new_start'] = np.where((df['start'] < df['end']), df['start'], df['end'])
    df['new_end'] = np.where((df['start'] < df['end']), df['end'], df['start'])
    df.drop(columns = ['start', 'end'], inplace = True)
    df.rename(columns = {'new_start': 'start', 'new_end': 'end'}, inplace = True)

    return df

# def create_windows(df, window_size):
#     """
#     """
    # There should be one window per custom annotation, with the size of the
    # window determined by --window_size
    



#-------------------------------------------------------------------------------
# main()
#-------------------------------------------------------------------------------

def main(args):
    custom_df = parse_custom_annot(custom_annot_file = args.custom)
    prokka_df = parse_prokka(parsed_gbk_file = args.prokka)
    df = merge_dfs(custom_df, prokka_df, annot_field = args.annot_field)
    print(df)
    # create_windows(df, window_size = args.window_size)
    df.to_csv("combined_df.tsv", sep = "\t", index = False)

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    args = parse_arguments()
    main(args)