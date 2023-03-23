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

# parser = parse_args(OptionParser(option_list = option_list))

# TESTING
parser = parse_args(OptionParser(option_list=option_list),
           args = c("--input=/home/brownli/Desktop/gencoflow/example/multi-sample/gene-windows.tsv",
                    "--outdir=/home/brownli/Desktop/gencoflow/example/multi-sample/"))

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

##################
# Store ggplot objects in a dataframe
test = windows %>% 
  group_nest(sample_id) %>%
  # Used for determining height of the resulting figure
  mutate(num_windows = map(data, .f = \(x){
    x %>% 
      select(seq_id) %>% 
      unique() %>% 
      nrow()
  })) %>% 
  mutate(gg = map(data, .f = \(x){
    gggenomes(seqs = x %>%
                group_by(seq_id) %>% 
                mutate(length = max(end) - min(start)) %>% 
                select(seq_id, contig_id, length) %>% 
                unique(), 
              genes = x %>%
                filter(source == 'prokka'), 
              feats = x %>% 
                filter(source == 'target')) +
      geom_seq() +
      geom_seq_label(aes(label = paste0(seq_id, "_", contig_id)), vjust = -10) +
      geom_gene(fill = 'grey80', position = position_strandpile(offset = 0.1)) +
      geom_gene_note(aes(label = gene), position = position_strandpile(offset = 0.1), vjust = 0, hjust = -0.3, size = 1.75) +
      geom_feat(linewidth = 1, colour = 'dodgerblue', position = position_strandpile(offset = 0.2)) +
      geom_feat_note(aes(label = target_name), position = position_strandpile(offset = 0.2))
  }))

# # Plot ggplot plots
# test$gg

for (i in 1:length(unique(test$sample_id))) {
  test$gg[i] 
  ggsave(filename = paste0(file.path(output_dir), test$sample_id[i], '.pdf'),
         height = as.integer(test$num_windows[i]) * 1.3,
         width = 12,
         limitsize = FALSE)
}



###############

































# Used for determining height of the resulting figure
num_windows = windows %>%
  select(seq_id) %>% 
  unique() %>% 
  nrow()

# Gene track
s0 = windows %>%
  group_by(seq_id) %>% 
  mutate(length = max(end) - min(start)) %>% 
  select(seq_id, contig_id, length) %>% 
  unique()

# Prokka annotations
p0 = windows %>% 
  filter(source == 'prokka')

# Target annotations
t0 = windows %>% 
  filter(source == 'target')

# If the --pile option was inputted, use position_strandpile for geom_gene
if (parser$pile == FALSE) {
  gg = gggenomes(seqs = s0, genes = p0, feats = t0) +
    geom_seq() +
    geom_seq_label(aes(label = paste0(seq_id, "_", contig_id)), vjust = -5) +
    geom_gene(fill = 'grey80') +
    geom_gene_note(aes(label = gene), vjust = 0, hjust = -0.3, size = 1.75) +
    geom_feat(linewidth = 1, colour = 'dodgerblue', position = position_strandpile(offset = 0.15)) +
    geom_feat_note(aes(label = arg), position = position_strandpile(offset = 0.15))

  ggsave(paste0(file.path(output_dir), "gene-windows.pdf"),
         plot = gg,
         height = (num_windows * 0.8),
         width = 12,
         limitsize = FALSE)
} else {
  gg = gggenomes(seqs = s0, genes = p0, feats = t0) +
    geom_seq() +
    geom_seq_label(aes(label = paste0(seq_id, "_", contig_id)), vjust = -10) +
    geom_gene(fill = 'grey80', position = position_strandpile(offset = 0.1)) +
    geom_gene_note(aes(label = gene), position = position_strandpile(offset = 0.1), vjust = 0, hjust = -0.3, size = 1.75) +
    geom_feat(linewidth = 1, colour = 'dodgerblue', position = position_strandpile(offset = 0.2)) +
    geom_feat_note(aes(label = arg), position = position_strandpile(offset = 0.2))

  ggsave(paste0(file.path(output_dir), "gene-windows.pdf"),
         plot = gg,
         height = (num_windows * 1.3),
         width = 12,
         limitsize = FALSE)
}