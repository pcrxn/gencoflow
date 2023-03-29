#!/usr/bin/env python

"""
gencoflow.py

Wrapper script for 'parse_gbk.py', 'create_windows.py' and 'vis_windows.R'.

Parse a Prokka-outputted GenBank file of genomic features into a TSV file 
containing only ORF annotations for downstream use in GenCoFlow.

Create target windows using target annotations and parsed Prokka output.

# TODO: better description
Create output figure

# TODO: add R dependencies
Dependencies: python=3.11.0, biopython=1.80, pandas=1.5.2, 
Other package versions may work but are untested.
"""

__author__ = "Liam Brown"
__email__ = "liam.brown@inspection.gc.ca"
__license__ = "Crown Copyright 2023"

import os
import sys
import logging
import argparse
from subprocess import call

from parse_gbk import main as parse_gbk_main
from create_windows import main as create_windows_main


def parse_arguments():
    """
    Parse command-line arguments.

    :return list args: List of parsed arguments.
    """

    parser = argparse.ArgumentParser(
        description = """
        Create target windows using target annotations and parsed Prokka output.
        """)
        

    # Required arguments
    required_args = parser.add_argument_group('Required')
    # Parse GBK file
    required_args.add_argument('-i', '--input', type = str, required = True,
        help = """
        Path to a GenBank file to parse.
        """)
    required_args.add_argument('-r', '--report_file', type = str, required = True,
        help = """
        Path for a new TSV file.
        """)
    # Create windows
    required_args.add_argument('-p', '--prokka', type = str, required = True,
        help = """
        Path to a TSV file produced by parse_gbk.py.
        """)
    required_args.add_argument('-t', '--targets', type = str, required = True,
        help = """
        Path to a TSV file of target annotations, such as those produced by
        an ARG identifier. The first row must contain headers for columns, and
        the following headers must be present: "contig_id", "start", "end", 
        "strand".
        """)
    required_args.add_argument('-o', '--output', type = str, required = True,
        help = """
        Path for exporting target windows and annotations in TSV format.
        """)
    required_args.add_argument('-n', '--target_name', type = str, 
                               required = True,
        help = """
        Name of a column -t/--targets file to use for annotation of the targets.
        """)    
    
    # Optional arguments
    optional_args = parser.add_argument_group('Optional')
    optional_args.add_argument('-s', '--sample_id', type = str, 
                            required = False,
        help = """
        Column name in -t/--targets file which contains sample identifiers. If
        none provided, GenCoFlow will run in single-sample mode and generate a
        single plot of the target windows. If a column name is provided, 
        GenCoFlow will run in multi-sample mode and generate a plot of gene 
        windows for each unique sample identifier.
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
    optional_args.add_argument('-w', '--window_size', type = int, 
                               required = False, default = 3000,
        help = """
        Number of bp flanking the start and end positions of each target
        annotation to include for analysis.
        Default: 3000
        """)
    optional_args.add_argument('-v', '--verbosity',
        choices=['debug', 'info', 'warning', 'error', 'critical'],
        metavar='verbosity',
        default='info',
        help = """
        Set the logging level. Options are debug, info, warning, error, and
        critical.
        Default: 'info'
        """
    )
    
    
    # Vis windows
    optional_args.add_argument('-pile', '--pile', action='store_true',
        help = """
        Save figure in strandpile format
        """)


    # If no arguments provided:
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    # Get the arguments into an object
    args = parser.parse_args()
    
    logging.basicConfig(level=args.verbosity.upper())

    return args


def main(args):
    """
    Run the three scripts in the workflow

    :param args: Argument parser object with required arguments
    """
    logging.info('Validating arguments')
    # Ensure that the supplied arguments are valid
    parse_gbk_check(args)

    # Parse the supplied GBK file
    parse_gbk_main(args)

    # Create windows
    create_windows_main(args)
    
    # Create figure
    call('vis_windows.R')

def parse_gbk_check(args):
    """
    Perform checks to ensure that the supplied arguments are valid

    :param args: Argument parser object with required arguments
    """
    # List to store all the errors encountered
    errors = []
    # Input
    if not os.path.isfile(args.input):
        errors.append(
                f'Could not locate the supplied GBK file, {args.input}. Please ensure that it '
                'exists.'
            )
    report_path = os.path.basename(args.report_file)
    if not os.path.isdir(report_path):
        try:
            os.makedirs(report_path)
            logging.warning(
                'Directory for supplied report file, %s did not exist. It has been created for '
                'you.', args.report_file
                )
        except IOError as exc:
            errors.append(
                f'Could not create the necessary folder into which the report is to be written: ' 
                f' {report_path}. Please ensure that you have the necessary permissions.\n'
                f'Error: {exc}'
            )
    if errors:
        error_writing(errors)

def error_writing(errors):
    """
    Write the grammatically correct error message to terminal

    :param errors: List of errors encountered with the arguments
    """
    error_string = '\n'.join(errors)
    was_were = 'was' if len(errors) == 1 else 'were'
    correct = 'error' if len(errors) == 1 else 'errors'
    logging.error(
        'There %s %s %s when attempting to run your command: \n%s', 
        was_were, len(errors), correct, error_string)
    raise SystemExit

if __name__ == '__main__':
    args = parse_arguments()
    main(args)