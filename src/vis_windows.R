#!/usr/bin/env/R

#-------------------------------------------------------------------------------
# Install gggenomes dependencies
#-------------------------------------------------------------------------------
if (!requireNamespace("ggtree", quietly = TRUE))
  BiocManager::install("ggtree")
if (!requireNamespace("thacklr", quietly = TRUE))
  devtools::install_github("thackl/thacklr")
if (!requireNamespace("gggenomes", quietly = TRUE))
  devtools::install_github("thackl/gggenomes")

#-------------------------------------------------------------------------------
# Load libraries
#-------------------------------------------------------------------------------
library(optparse)
library(tidyverse)
library(gggenomes)

#-------------------------------------------------------------------------------
# Parse command-line arguments
#-------------------------------------------------------------------------------
parser = OptionParser()
parser = add_option(parser,
                    c("-i", "--input"),
                    type = "character",
                    # metavar = "path/to/gene_windows.tsv",
                    help = "Path to input TSV file of gene windows produced by
                    `create_windows.py`.")
parser = add_option(parser,
                    c("-o", "--outdir"),
                    type = "character",
                    metavar = "path/to/outdir/",
                    help = "Path to directory for saving plots.")

parse_args(parser)

ifelse(!is.na(parser$input),
       input_tsv <- parser$input,
       stop("Input TSV file must be provided. See script usage (--help)"))

ifelse(!is.na(parser$outdir),
       input_tsv <- parser$outdir,
       stop("Input TSV file must be provided. See script usage (--help)"))