#!/usr/bin/env python

"""
Tests for parse_gbk.py
"""

__author__ = "Liam Brown"
__email__ = "liam.brown@inspection.gc.ca"
__license__ = "Crown Copyright 2023"

import argparse
from unittest.mock import patch

import pytest

from src.parse_gbk import \
    determine_run_mode, \
    parse_arguments

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------

# determine_run_mode()
@pytest.mark.parametrize('input,directory,expected_result',
                         [(True, None, 'single'),
                          (None, True, 'multi'),
                          (None, None, None),
                          (True, True, None),
                          ('path/prokka_annots.gbk', None, 'single'),
                          (None, 'path/to/gbks/', 'multi')])
def test_determine_run_mode(input, directory, expected_result):
    run_mode = determine_run_mode(
        input=input,
        dir=directory
    )
    assert run_mode == expected_result

#-------------------------------------------------------------------------------
# Argument parsing
#-------------------------------------------------------------------------------

# # --empty args
# @patch('argparse.ArgumentParser.parse_args')
# def test_arguments_empty_args(mock_args):
#     mock_args.return_value = argparse.Namespace()
#     args = parse_arguments()
#     # Ensure that args is empty
#     assert not args
    

# --no report
@patch('argparse.ArgumentParser.parse_args')
def test_arguments_no_report(mock_args):
    mock_args.return_value = argparse.Namespace(
        report=None,
        input='test'
    )
    args = parse_arguments()