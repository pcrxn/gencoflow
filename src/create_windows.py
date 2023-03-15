#!/usr/bin/env python

"""
create_windows.py

Create gene windows using target annotations and parsed Prokka output.

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
        Create gene windows using target annotations and parsed Prokka output.
        """)

    # Required arguments
    required_args = parser.add_argument_group('Required')
    required_args.add_argument('-p', '--prokka', type = str, required = True,
        help = """
        Path to a TSV file produced by parse_gbk.py.
        """)
    required_args.add_argument('-t', '--targets', type = str, required = True,
        help = """
        Path to a TSV file of target annotations, such as those produced by
        an ARG identifier. The first row must contain headers for columns, and
        the following headers must be present:"contig_id", "start", "end", 
        "strand".
        """)
    required_args.add_argument('-o', '--output', type = str, required = True,
        help = """
        Path for a new TSV file of gene windows for downstream analysis.
        """)
    required_args.add_argument('-n', '--target_name', type = str, 
                               required = False,
        help = """
        Name of a column in -t/--targets to use for annotation.
        """)    
    
    # Optional arguments
    optional_args = parser.add_argument_group('Optional')
    optional_args.add_argument('-w', '--window_size', type = int, 
                               required = False, default = 3000,
        help = """
        Number of bp flanking the start and end positions of each target
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

def parse_target_annot(target_annot_file):
    """
    """
    target_df = pd.read_csv(target_annot_file, sep = '\t', index_col = False)
    target_df['source'] = 'target'
    target_df['strand'].replace('+', 1, inplace = True)
    target_df['strand'].replace('-', -1, inplace = True)

    return target_df

def parse_prokka(parsed_gbk_file):
    """
    """
    prokka_df = pd.read_csv(parsed_gbk_file, sep = '\t', index_col = False)
    prokka_df['source'] = 'prokka'

    return prokka_df

def merge_dfs(target_df, prokka_df, target_name):
    """
    """
    df = pd.concat([
        prokka_df[['contig_id', 'source', 'start', 'end', 'strand']],
        target_df[['contig_id', 'source', 'start', 'end', 'strand', target_name]]
    ])

    # Fix start and end positions if necessary
    # Start should always be less than end, regardless of strand
    df['new_start'] = np.where((df['start'] < df['end']), df['start'], df['end'])
    df['new_end'] = np.where((df['start'] < df['end']), df['end'], df['start'])
    df.drop(columns = ['start', 'end'], inplace = True)
    df.rename(columns = {'new_start': 'start', 'new_end': 'end'}, inplace = True)

    return df

def obtain_windows_list(df, window_size):
    """
    """
    windows_list = []

    # Store each target and its terminal window coordinates as a dictionary
    # within `windows_list`
    for target in df[df['source'] == 'target'].to_dict(orient = 'records'):

        window_handle = {}

        # Target info
        window_handle['window_contig_id'] = target['contig_id']
        window_handle['window_target_strand'] = target['strand']
        # window_handle['window_target_source'] = target['source']

        # window_start
        if target['start'] - window_size > 0:
            window_handle['window_start'] = target['start'] - window_size
        else: 
            window_handle['window_start'] = 0

        # window_end
        window_handle['window_end'] = target['end'] + window_size
        windows_list.append(window_handle)
    
    return windows_list

def create_wdf(df, windows_list, target_name):
    """
    Use windows_list to separate `df` by target windows.
    """

    wdf = pd.DataFrame()
    window_count = 0

    for window in windows_list:

        window_count += 1

        # Filter `df` for ORFs within the current window scope
        wdf_handle = df[(df['start'] >= window['window_start']) & 
            (df['end'] <= window['window_end']) &
            (df['contig_id'] == window['window_contig_id'])].copy()

        # Add target info to the dataframe
        wdf_handle['window_target_strand'] = window['window_target_strand']
        # wdf_handle['window_target_source'] = window['window_target_source']
        wdf_handle['orig_window_start'] = window['window_start']
        wdf_handle['orig_window_end'] = window['window_end']
        wdf_handle['window_id'] = window_count

        # Recode start and stop positions to use window-specific coordinates
        wdf_handle = wdf_handle.rename(columns = {'start': 'orig_start', 'end': 'orig_end'})
        wdf_handle['start'] = wdf_handle['orig_start'] - min(wdf_handle['orig_start'])
        wdf_handle['end'] = min(wdf_handle['orig_end']) - min(wdf_handle['orig_start']) + wdf_handle['start']

        wdf = pd.concat([wdf, wdf_handle], ignore_index = True)

    # Reorder columns
    wdf = wdf[['window_id', 'contig_id', 'source', 'start', 'end', 'strand',
                target_name, 'orig_start', 'orig_end', 'orig_window_start',
                'orig_window_end', 'window_target_strand']]

    return wdf

def save_wdf(wdf, output):
    """
    Save dataframe of windows for downstream analysis.
    """
    wdf.to_csv(output, sep = '\t', index = False)

#-------------------------------------------------------------------------------
# main()
#-------------------------------------------------------------------------------

def main(args):
    target_df = parse_target_annot(target_annot_file = args.targets)
    prokka_df = parse_prokka(parsed_gbk_file = args.prokka)
    df = merge_dfs(target_df, prokka_df, target_name = args.target_name)
    # print(df)
    windows_list = obtain_windows_list(df, window_size = args.window_size)
    wdf = create_wdf(df, windows_list, target_name = args.target_name)
    # print(wdf)
    save_wdf(wdf, output = args.output)

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    args = parse_arguments()
    main(args)