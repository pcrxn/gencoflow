#!/usr/bin/env python

"""
parse_gbk.py

Parse a Prokka-outputted GenBank file of genomic features into a TSV file 
containing only ORF annotations for downstream use in GenCoFlow.

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
        Parse a Prokka-outputted GenBank file of genomic features into a TSV
        file containing only ORF annotations for downstream use in GenCoFlow.       
        """)

    # Required arguments
    required_args = parser.add_argument_group('Required')
    required_args.add_argument('-i', '--input', type = str, required = True,
        help = """
        Path to a GenBank file to parse.
        """)
    required_args.add_argument('-o', '--output', type = str, required = True,
        help = """
        Path for a new TSV file.
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

def parse_gbk(genbank_file):
    """
    Parse the contents of a GenBank file and extract only CDS feature info.
    :param genbank_file: Path to a GenBank file to parse.
    :type genbank_file: str
    :returns orfs: List of dictionaries, where each dictionary holds details for
    a single ORF.
    :rtype orfs: list
    """
    orfs = []
    for record in SeqIO.parse(genbank_file, 'genbank'):
        for feature in record.features:
            # Extract ORF info only
            if feature.type == 'CDS':
                # Create a dict to hold ORF info and minimal parent contig info
                # Populate parent contig info
                orf_dict = dict()
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
    Convert list of ORFs to a Pandas dataframe.
    :param orfs: List of dictionaries, where each dictionary holds details for
    a single ORF.
    :type orfs: list
    :returns df: Dataframe of ORF information.
    :rtype df: <class 'pandas.core.frame.DataFrame'>
    """

    df = pd.DataFrame(orfs)
    return df

def export_tsv(df, output):
    """
    Export the Pandas dataframe of ORF details to a TSV file.
    :param df: Dataframe of ORF information.
    :type df: <class 'pandas.core.frame.DataFrame'>
    :param output: Path for a new TSV file.
    :type output: str
    """
    df.to_csv(output, sep = '\t', na_rep = 'NA', index = False)

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