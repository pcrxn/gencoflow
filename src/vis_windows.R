#!/usr/bin/env Rscript

#-------------------------------------------------------------------------------
# Install gggenomes dependencies
#-------------------------------------------------------------------------------

if (!requireNamespace("ggtree", quietly = TRUE))
  BiocManager::install("ggtree")
if (!requireNamespace("IRanges", quietly = TRUE))
  BiocManager::install("IRanges")
if (!requireNamespace("thacklr", quietly = TRUE))
  devtools::install_github("thackl/thacklr")
if (!requireNamespace("gggenomes", quietly = TRUE))
  devtools::install_github("thackl/gggenomes")

#-------------------------------------------------------------------------------
# Load libraries
#-------------------------------------------------------------------------------

library.path = .libPaths()
# library("optparse", lib.loc = library.path)
# library("tidyverse", lib.loc = library.path)
# library("gggenomes", lib.loc = library.path)
suppressMessages(library("optparse", lib.loc = library.path))
suppressMessages(library("tidyverse", lib.loc = library.path))
suppressMessages(library("gggenomes", lib.loc = library.path))

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
              help = "Path to directory for saving plots."),
  make_option(c("-p", "--pile"),
              action = "store_true",
              default = FALSE,
              help = "Save figure in strandpile format.
              Default: %default")
)

parser = parse_args(OptionParser(option_list = option_list))

# TESTING
# parser = parse_args(OptionParser(option_list=option_list),
#            args = c("--input=gencoflow/example/multi-sample/gene-windows.tsv",
#                     "--outdir=gencoflow/example/multi-sample/",
#                     "--pile"))
# parser = parse_args(OptionParser(option_list=option_list),
#                     args = c("--input=gencoflow/example/multi-sample/gene-windows.tsv",
#                              "--outdir=gencoflow/example/multi-sample/"))

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

windows = read_tsv(input_tsv) %>%
  # Rename the window_id column for compatibility with gggenomes
  rename(seq_id = window_id)

#-------------------------------------------------------------------------------
# Create figures
#-------------------------------------------------------------------------------

# In the future, make these into functions for brevity
for (i in 1:length(unique(windows$sample_id))) {
  num_windows = windows %>% 
    filter(sample_id == unique(windows$sample_id)[i]) %>%
    select(seq_id) %>% 
    unique() %>% 
    nrow()
  
  seqs = windows %>% 
    filter(sample_id == unique(windows$sample_id)[i]) %>% 
    group_by(seq_id) %>% 
    mutate(length = max(end) - min(start)) %>% 
    select(seq_id, contig_id, length) %>% 
    unique()
  
  prokka = windows %>% 
    filter(sample_id == unique(windows$sample_id)[i]) %>% 
    filter(source == 'prokka')
  
  targets = windows %>% 
    filter(sample_id == unique(windows$sample_id)[i]) %>% 
    filter(source == 'target')
  
  # --pile
  if (parser$pile == TRUE) {
    gg = gggenomes(seqs = seqs, genes = prokka, feats = targets) +
      geom_seq() +
        geom_seq_label(aes(label = paste0(seq_id, "_", contig_id)), vjust = -10) +
        geom_gene(fill = 'grey80', position = position_strandpile(offset = 0.1)) +
        geom_gene_note(aes(label = gene), position = position_strandpile(offset = 0.1), vjust = 0, hjust = -0.3, size = 1.75) +
        geom_feat(linewidth = 1, colour = 'dodgerblue', position = position_strandpile(offset = 0.2)) +
        geom_feat_note(aes(label = target_name), position = position_strandpile(offset = 0.2))
    
    # Save plots
    if (num_windows < 2){
      ggsave(paste0(file.path(output_dir), unique(windows$sample_id)[i], '.pdf'),
             plot = gg,
             height = 2,
             width = 12,
             limitsize = FALSE)
      ggsave(paste0(file.path(output_dir), unique(windows$sample_id)[i], '.png'),
             plot = gg,
             height = 2,
             width = 12,
             dpi = 300,
             limitsize = FALSE)
    } else {
      ggsave(paste0(file.path(output_dir), unique(windows$sample_id)[i], '.pdf'),
             plot = gg,
             height = as.integer(num_windows[i]) * 1.3,
             width = 12,
             limitsize = FALSE)
      ggsave(paste0(file.path(output_dir), unique(windows$sample_id)[i], '.png'),
             plot = gg,
             height = as.integer(num_windows[i]) * 1.3,
             width = 12,
             dpi = 300,
             limitsize = FALSE)
    }
  # No --pile
  } else {
    gg = gggenomes(seqs = seqs, genes = prokka, feats = targets) + 
      geom_seq() +
      geom_seq_label(aes(label = paste0(seq_id, "_", contig_id)), vjust = -5) +
      geom_gene(fill = 'grey80') +
      geom_gene_note(aes(label = gene), vjust = 0, hjust = -0.3, size = 1.75) +
      geom_feat(linewidth = 1, colour = 'dodgerblue', position = position_strandpile(offset = 0.15)) +
      geom_feat_note(aes(label = target_name), position = position_strandpile(offset = 0.15))
    
    # Save plots
    if (num_windows < 2){
      # pdf
      ggsave(paste0(file.path(output_dir), unique(windows$sample_id)[i], '.pdf'),
             plot = gg,
             height = 2,
             width = 12,
             limitsize = FALSE)
      # png
      ggsave(paste0(file.path(output_dir), unique(windows$sample_id)[i], '.png'),
             plot = gg,
             height = 2,
             width = 12,
             dpi = 300,
             limitsize = FALSE)
    } else {
      # pdf
      ggsave(paste0(file.path(output_dir), unique(windows$sample_id)[i], '.pdf'),
             plot = gg,
             height = as.integer(num_windows) * 0.8,
             width = 12,
             limitsize = FALSE)
      # png
      ggsave(paste0(file.path(output_dir), unique(windows$sample_id)[i], '.png'),
             plot = gg,
             height = as.integer(num_windows) * 0.8,
             width = 12,
             dpi = 300,
             limitsize = FALSE)
    }
  }
}
