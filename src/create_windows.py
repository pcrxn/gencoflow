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
import re
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
                               required = True,
        help = """
        Name of a column -t/--targets file to use for annotation of the targets.
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
    optional_args.add_argument('-s', '--sample_id', type = str, 
                            required = False,
        help = """
        Column name in -t/--targets file which contains sample identifiers. If
        none provided, GenCoFlow will run in single-sample mode and generate a
        single plot of the gene windows. If a column name is provided, GenCoFlow
        will run in multi-sample mode and generate a plot of gene windows for 
        each unique sample identifier.
        """)
    optional_args.add_argument('-c', '--contig_id', type = str, 
                        required = False, default = 'contig_id',
        help = """
        Column name in -t/--targets file which contains unique contig
        identifiers for each sample. Each unique contig identifier corresponds
        to a different gene track in the generated plot(s).
        Default: 'contig_id'
        """)
    optional_args.add_argument('-1', '--start', type = str, 
                        required = False, default = 'start',
        help = """
        Column name in -t/--targets which contains the start position of each
        target.
        Default: 'start'
        """)
    optional_args.add_argument('-2', '--end', type = str, 
                        required = False, default = 'end',
        help = """
        Column name in -t/--targets which contains the end position of each
        target.
        Default: 'end'
        """)
    optional_args.add_argument('-x', '--strand', type = str, 
                        required = False, default = 'strand',
        help = """
        Column name in -t/--targets which contains the strand of each target
        (-1, 1, +, or -).
        Default: 'strand'
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

def determine_run_mode(sample_id_name):
    """
    Determine whether GenCoFlow should be run in single-sample or multi-sample
    mode.
    """
    run_mode = None

    if sample_id_name is None:
        run_mode = 'single'
    else:
        run_mode = 'multi'

    return run_mode

def read_target_annot(run_mode, target_annot_file, sample_id_name, 
                      contig_id_name, start_name, end_name, strand_name,
                      target_name):
    """

    """
    target_df = pd.read_csv(target_annot_file, sep = '\t', index_col = False)

    if run_mode == 'single':
        target_df.rename(columns = {contig_id_name: 'contig_id',
                                    start_name: 'start',
                                    end_name: 'end',
                                    strand_name: 'strand',
                                    target_name: 'target_name'},
                                    inplace = True)
    elif run_mode == 'multi':
        target_df.rename(columns = {contig_id_name: 'contig_id',
                                    start_name: 'start',
                                    end_name: 'end',
                                    strand_name: 'strand',
                                    target_name: 'target_name',
                                    sample_id_name: 'sample_id'},
                                    inplace = True)

    target_df['source'] = 'target'
    try:
        target_df['strand'].replace('+', 1, inplace = True)
        target_df['strand'].replace('-', -1, inplace = True)
    except KeyError:
        raise KeyError("""'strand' not found in column headers of -t/--targets
        file. If the strand column has a different name, you can specify this
        using -x/--strand.
        """)

    return target_df

def read_prokka(parsed_gbk_file):
    """
    Read a TSV file of parsed Prokka annotations and store as a dataframe.

    :param parsed_gbk_file: Path to a TSV file produced by parse_gbk.py.
    :type parsed_gbk_file: str
    :returns prokka_df: A dataframe containing Prokka annotation information and
    coordinates.
    :rtype prokka_df: <class 'pandas.core.frame.DataFrame'>
    """
    prokka_df = pd.read_csv(parsed_gbk_file, sep = '\t', index_col = False)
    # Add prefixes to some column names to avoid conflicts with merge_dfs()
    prokka_df.rename(columns = {'topology': 'prokka_topology',
                                'type': 'prokka_type',
                                'orf_id': 'prokka_orf_id',
                                'inference': 'prokka_inference',
                                'codon_start': 'prokka_codon_start',
                                'transl_table': 'prokka_transl_table',
                                'product': 'prokka_product',
                                'translation': 'prokka_translation',
                                'db_xref': 'prokka_db_xref',
                                'gene': 'prokka_gene',
                                'ec_num': 'prokka_ec_num',
                                'note': 'prokka_note'})
    prokka_df['source'] = 'prokka'

    return prokka_df

def merge_dfs(target_df, prokka_df, run_mode):
    """
    Combine Prokka and target annotations into one dataframe, and fix the start
    and end coordinates if necessary.

    :param target_df: A dataframe containing target annotation information and
    coordinates.
    :param prokka_df:
    :param target_name:
    """

    if run_mode == 'single':
        df = pd.concat([
            prokka_df,
            target_df[['contig_id', 'source', 'start', 'end', 'strand', 
                    'target_name']]
        ])
    elif run_mode == 'multi':
        df = pd.concat([
        prokka_df,
        target_df[['sample_id', 'contig_id', 'source', 'start', 'end', 'strand', 
                'target_name']]
    ])

    # Fix start and end positions if necessary
    # Start should always be less than end, regardless of strand
    df['new_start'] = np.where((df['start'] < df['end']), df['start'], 
                               df['end'])
    df['new_end'] = np.where((df['start'] < df['end']), df['end'], df['start'])
    df.drop(columns = ['start', 'end'], inplace = True)
    df.rename(columns = {'new_start': 'start', 'new_end': 'end'},
              inplace = True)

    return df

def obtain_windows_list(df, run_mode, window_size):
    """
    """
    windows_list = []

    # Store each target and its terminal window coordinates as a dictionary
    # within `windows_list`
    for target in df[df['source'] == 'target'].to_dict(orient = 'records'):

        window_handle = {}

        # Target info
        if run_mode == 'multi':
            window_handle['window_sample_id'] = target['sample_id']
        elif run_mode == 'single':
            continue

        # More target info
        window_handle['window_contig_id'] = target['contig_id']
        window_handle['window_target_strand'] = target['strand']

        # window_start
        if target['start'] - window_size > 0:
            window_handle['window_start'] = target['start'] - window_size
        else: 
            window_handle['window_start'] = 0

        # window_end
        window_handle['window_end'] = target['end'] + window_size
        windows_list.append(window_handle)
    
    return windows_list

def create_wdf(df, windows_list, run_mode):
    """
    Use windows_list to separate `df` by target windows.
    """

    wdf = pd.DataFrame()
    window_count = 0

    for window in windows_list:

        window_count += 1

        # Filter `df` for ORFs within the current window scope
        if run_mode == 'single':
            wdf_handle = df[(df['start'] >= window['window_start']) & 
                (df['end'] <= window['window_end']) &
                (df['contig_id'] == window['window_contig_id'])].copy()
        elif run_mode == 'multi':
            wdf_handle = df[(df['start'] >= window['window_start']) & 
                (df['end'] <= window['window_end']) &
                (df['contig_id'] == window['window_contig_id']) &
                (df['sample_id'] == window['window_sample_id'])].copy()

        # Add target info to the dataframe
        wdf_handle['window_target_strand'] = window['window_target_strand']
        wdf_handle['orig_window_start'] = window['window_start']
        wdf_handle['orig_window_end'] = window['window_end']
        wdf_handle['window_id'] = window_count

        # Recode start and stop positions to use window-specific coordinates
        wdf_handle = wdf_handle.rename(columns = {'start': 'orig_start',
                                                  'end': 'orig_end'})
        wdf_handle['start'] = wdf_handle['orig_start'] - min(
            wdf_handle['orig_start'])
        wdf_handle['end'] = min(wdf_handle['orig_end']) - min(
            wdf_handle['orig_start']) + wdf_handle['start']

        wdf = pd.concat([wdf, wdf_handle], ignore_index = True)

    # Reorder columns
    if run_mode == 'single':
        wdf = wdf[['window_id', 'contig_id', 'source', 'start', 'end', 'strand',
                    'target_name', 'orf_id', 'inference', 'product', 'db_xref', 
                    'gene', 'note', 'orig_start', 'orig_end', 'orig_window_start',
                    'orig_window_end', 'window_target_strand']]
    if run_mode == 'multi':
        wdf = wdf[['window_id', 'sample_id', 'contig_id', 'source', 'start',
                    'end', 'strand', 'target_name', 'orf_id', 'inference',
                     'product', 'db_xref', 'gene', 'note', 'orig_start',
                      'orig_end', 'orig_window_start', 'orig_window_end',
                      'window_target_strand']]

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
    run_mode = determine_run_mode(sample_id_name = args.sample_id)
    target_df = read_target_annot(run_mode,
                                  target_annot_file = args.targets,
                                  sample_id_name = args.sample_id,
                                  contig_id_name = args.contig_id,
                                  start_name = args.start,
                                  end_name = args.end,
                                  strand_name = args.strand,
                                  target_name = args.target_name)
    prokka_df = read_prokka(parsed_gbk_file = args.prokka)
    df = merge_dfs(target_df, prokka_df, run_mode)
    # print(df)
    windows_list = obtain_windows_list(df, run_mode,
                                       window_size = args.window_size)
    wdf = create_wdf(df, windows_list, run_mode)
    # print(wdf)
    save_wdf(wdf, output = args.output)

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    args = parse_arguments()
    main(args)