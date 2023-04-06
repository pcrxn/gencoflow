#!/usr/bin/env python

"""
gencoflow.py

Wrapper script for 'parse_gbk.py', 'create_windows.py' and 'vis_windows.R'.

Create visualizations of windows surrounding gene targets of interest using
Prokka output and custom target annotations.

Dependencies: 
    - python=3.11.0
        - biopython=1.80, pandas=1.5.2
    - R=4.2.3
        - ggtree=3.6.2
        - IRanges=2.32.0
        - thacklr=0.0.0.9000
        - gggenomes=0.9.5.9000

Other package versions may work but are untested.
"""

__author__ = "Liam Brown"
__email__ = "liam.brown@inspection.gc.ca"
__license__ = "Crown Copyright 2023"

import os
import sys
import logging
import argparse
from pathlib import Path
from subprocess import call

from src.version import __version__
from parse_gbk import main as parse_gbk_main
from create_windows import main as create_windows_main

#-------------------------------------------------------------------------------
# parse_arguments()
#-------------------------------------------------------------------------------

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
    required_args.add_argument('-r', '--report_file', type = str, required = True,
        help = """
        Path for exporting data as a new TSV report.
        Example: 'path/new.tsv'
        """)
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
    
    # Either arguments
    either_args = parser.add_mutually_exclusive_group(required = True)
    either_args.add_argument('-i', '--input', type = str,
        help = """
        Path to a GenBank file to parse (single-sample mode).
        Example: 'path/prokka_annots.gbk'
        """)
    either_args.add_argument('-d', '--dir', type = str,
        help = """
        Path to a directory to recursively search within for GenBank files
        (multi-sample mode). The base file name of each GenBank file will become
        the sample_id.
        Example: 'path/to/gbks/'
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
    optional_args.add_argument(
        '-version', '--version',
        action='version',
        version=f'%(prog)s commit {__version__}'
        )
    optional_args.add_argument('-pile', '--pile', action='store_true',
        help = """
        Save figure in strandpile format
        """)


    # If no arguments provided:
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    
    logging.basicConfig(level=args.verbosity.upper())

    return args

#-------------------------------------------------------------------------------
# Check functions
#-------------------------------------------------------------------------------

def parse_gbk_check(args):
    """
    Perform checks to ensure that the supplied arguments are valid.

    :param args: Argument parser object with required arguments.
    """

    # List to store all the errors encountered
    errors = []

    # Input
    if args.input:
        if not os.path.isfile(args.input):
            errors.append(f"""
            Could not locate the supplied GBK file, {args.input}. Please ensure
            that it exists.            
            """)
    elif args.dir:
        if not os.path.isdir(args.dir):
            errors.append(f"""
            Could not locate the supplied directory of GBK files, {args.dir}.
            Please ensure that it exists.
            """)
        else:
            gbk_files = []
            for path in Path(args.dir).rglob('*.gbk'):
                gbk_files.append(path)
            if not gbk_files:
                errors.append(f"""
                Could not locate any GBK files in the supplied directory, 
                {args.dir}. Please ensure that you have supplied the correct 
                folder name.
                """)

    # Report file
    report_path = os.path.basename(args.report)
    if not os.path.isdir(report_path):
        try:
            os.makedirs(report_path)
            logging.warning(f"""
            Directory for supplied report file, %s did not exist. It has been
            created for you: {args.report}
            """)
        except IOError as exc:
            errors.append(f"""
            Could not create the necessary folder into which the report is to be
            written: {report_path}. Please ensure that you have the necessary 
            permissions.\n
            """)
            logging.debug(exc)
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
    logging.error(f"""
    There {was_were} {len(errors)} {correct} when attempting to run your
    command: \n{error_string}
    """)

    raise SystemExit

#-------------------------------------------------------------------------------
# main()
#-------------------------------------------------------------------------------

def main(args):
    """
    Run the three scripts in the workflow.

    :param args: Argument parser object with required arguments.
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

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    args = parse_arguments()
    main(args)