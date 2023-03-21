#!/usr/bin/env Rscript

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

option_list <- list( 
  make_option(c("-i", "--input"),
              type = "character",
              help = "Path to input TSV file of gene windows produced by 
              `create_windows.py`."),
  make_option(c("-o", "--outdir"), 
              type = "character",
              metavar = "path/to/outdir/",
              help = "Path to directory for saving plots.")
)

# parser = parse_args(OptionParser(option_list = option_list))

# TESTING
parser = parse_args(OptionParser(option_list=option_list),
           args = c("--input=/home/brownli/Desktop/gencoflow/example/gene-windows.tsv",
                    "--outdir=/home/brownli/Desktop/gencoflow/example/"))

if (!is.null(parser$input)) {
  input_tsv <- parser$input
} else {
  stop("Input TSV file must be provided. See script usage (--help)")
}

if (!is.null(parser$outdir)) {
  output_dir <- parser$outdir
} else {
  stop("Output directory path must be provided. See script usage (--help)")
}

#-------------------------------------------------------------------------------
# Import gene windows
#-------------------------------------------------------------------------------

windows = read_tsv(input_tsv)

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
