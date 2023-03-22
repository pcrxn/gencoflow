#!/usr/bin/env python

"""
parse_gbk.py

Parse one or many Prokka-outputted GenBank files of genomic features into a TSV
file containing only ORF annotations for downstream use in GenCoFlow.

Dependencies: python=3.11.0, biopython=1.80, pandas=1.5.2
Other package versions may work but are untested.
"""

__author__ = "Liam Brown"
__email__ = "liam.brown@inspection.gc.ca"
__license__ = "Crown Copyright 2023"

import os
import sys
import argparse
import pandas as pd
from pathlib import Path
from Bio import SeqIO

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
        Parse one or many Prokka-outputted GenBank files of genomic features
        into a TSV file containing only ORF annotations for downstream use in
        GenCoFlow.       
        """)

    # Required arguments
    required_args = parser.add_argument_group('Required')
    required_args.add_argument('-o', '--output', type = str, required = True,
        help = """
        Path for a new TSV file.
        """)
    
    # Either arguments
    either_args = parser.add_mutually_exclusive_group(required = True)
    either_args.add_argument('-i', '--input', type = str,
        help = """
        Path to a GenBank file to parse (single-sample mode).
        """)
    either_args.add_argument('-d', '--dir', type = str,
        help = """
        Path to a directory to recursively search within for GenBank files
        (multi-sample mode).
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

def determine_run_mode(input, dir):
    """
    Determine whether GenCoFlow should be run in single-sample or multi-sample
    mode.
    """
    run_mode = None

    if input is not None and dir is None:
        run_mode = 'single'
    elif input is None and dir is not None:
        run_mode = 'multi'

    return run_mode

def obtain_genbank_files(run_mode, input, dir):
    """
    Store the path(s) of the GenBank file(s) provided on the command-line to
    --input or --dir in a list.
    """

    genbank_files = []
    if run_mode == 'single':
        genbank_files.append(os.path.abspath(input))
    elif run_mode == 'multi':
        for path in Path(dir).rglob("*.gbk"):
            genbank_files.append(os.path.abspath(path))

    return genbank_files

def parse_gbk(genbank_files):
    """

    """

    orfs = []
    for file in genbank_files:
        for record in SeqIO.parse(file, 'genbank'):
            sample_id = os.path.splitext(os.path.basename(file))[0]
            for feature in record.features:
                # Extract ORF info only
                if feature.type == 'CDS':
                    # Create a dict to hold ORF info and minimal parent contig info
                    # Populate parent contig info
                    orf_dict = dict()
                    orf_dict['sample_id'] = sample_id
                    orf_dict['contig_id'] = record.id
                    orf_dict['topology'] = record.annotations['topology']
                    # Populate ORF info
                    orf_dict['type'] = feature.type
                    orf_dict['start'] = int(feature.location.start)
                    orf_dict['end'] = int(feature.location.end)
                    orf_dict['strand'] = feature.strand
                    # orf_dict['id'] = feature.id
                    orf_dict['orf_id'] = feature.qualifiers.get('locus_tag', 'NA')
                    orf_dict['inference'] = feature.qualifiers.get('inference', 'NA')
                    orf_dict['codon_start'] = int(feature.qualifiers.get('codon_start', 'NA')[0])
                    orf_dict['transl_table'] = int(feature.qualifiers.get('transl_table', 'NA')[0])
                    orf_dict['product'] = feature.qualifiers.get('product', 'NA')
                    orf_dict['translation'] = feature.qualifiers.get('translation', 'NA')
                    orf_dict['db_xref'] = feature.qualifiers.get('db_xref', 'NA')
                    orf_dict['gene'] = feature.qualifiers.get('gene', 'NA')
                    orf_dict['ec_num'] = feature.qualifiers.get('EC_number', 'NA')
                    orf_dict['note'] = feature.qualifiers.get('note', 'NA')
                    # If any of the values in orf_dict are lists of length 1, 
                    # unlist:
                    for key, value in orf_dict.items():
                        if isinstance(value, list):
                            if len(value) == 1:
                                orf_dict[key] = value[0]
                    orfs.append(orf_dict)

    return orfs

def orfs_to_df(orfs):
    """

    """

    df = pd.DataFrame(orfs)
    return df

def export_tsv(df, output):
    """

    """
    df.to_csv(output, sep = '\t', na_rep = 'NA', index = False)

#-------------------------------------------------------------------------------
# main()
#-------------------------------------------------------------------------------

def main(args):
    run_mode = determine_run_mode(input = args.input, dir = args.dir)
    genbank_files = obtain_genbank_files(run_mode, input = args.input,
                                         dir = args.dir)
    orfs = parse_gbk(genbank_files)
    df = orfs_to_df(orfs)
    export_tsv(df, output = args.output)

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    args = parse_arguments()
    main(args)