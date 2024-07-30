suppressPackageStartupMessages(suppressWarnings(library(optparse)))
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))
suppressPackageStartupMessages(suppressWarnings(library(scales)))
suppressPackageStartupMessages(suppressWarnings(library(ggpubr)))
suppressPackageStartupMessages(suppressWarnings(library(ggsci)))

############################################################
# OPTIONS                                                  #
############################################################

option_list = list(
  make_option(c("-d", "--dir"), type = "character",
              help = "full path to repository", metavar = "character"),
  make_option(c("--bassik"), type = "character",
              help = "full path to Bassik results", metavar = "character"),
  make_option(c("--bassik_binary"), type = "character",
              help = "full path to Bassik binary results", metavar = "character"),
  make_option(c("-m", "--mapping"), type = "character",
              help = "full path to sample mapping", metavar = "character"),
  make_option(c("-r", "--rds"), type = "character",
              help = "full path to RDS directory", metavar = "character"),
  make_option(c("-o", "--out"), type = "character",
              help = "full path to output directory", metavar = "character"))

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

############################################################
# GENERAL                                                  #
############################################################

# Set top level directory
repo_path <- opt$dir

# Check top level directory exists
if (!dir.exists(repo_path)) {
  stop(smessagef("Repository directory not exist: %s", repo_path))
}

# Add helper functions
helper_path <- ifelse(is.null(opt$helper), file.path(repo_path, 'SCRIPTS', 'postprocessing', 'helper.R'), opt$helper)
if (!file.exists(helper_path)){
  stop(paste('Helper file does not exist:', helper_path))
}
source(helper_path)

# Check output path exists
check_dir_exists(opt$out)

# Set RDS path
rds_path <- ifelse(is.null(opt$rds), file.path(repo_path, 'DATA', 'RDS', 'QC'), opt$rds)
check_dir_exists(rds_path)

############################################################
# Bassik results matrix                                    #
############################################################

message('Reading in Bassik binary matrix...')
# Read in binary matrix
bassik_hits <- read.delim(opt$bassik, header = T, sep = "\t", check.names = F)

############################################################
# Bassik binary matrix                                     #
############################################################

message('Reading in Bassik binary matrix...')
# Read in binary matrix
bassik_binary_hits <- read.delim(opt$bassik_binary, header = T, sep = "\t", check.names = F)

############################################################
# Sample mapping                                           #
############################################################

check_file_exists(opt$mapping)

message(paste("Reading sample annotations from:", opt$mapping))

# Read in sample mapping
sample_mapping <- read_sample_metadata(opt$mapping)

# Cancer type palette
cancer_type_pal <- pal_npg(alpha = 0.3)(3)
names(cancer_type_pal) <- c(levels(sample_mapping$cancer_type))

############################################################
# Number of Bassik hits per cell line                      #
############################################################

message('Plotting number of Bassik hits per cell line...')

bassik_hits_per_cell_line <- bassik_binary_hits %>%
  select(sorted_gene_pair, 'A-375':'SU.86.86') %>%
  gather(cell_line_label, is_bassik_hit, -sorted_gene_pair) %>%
  filter(is_bassik_hit == 1) %>%
  group_by(cell_line_label) %>%
  summarise(n = sum(is_bassik_hit), .groups = 'keep') %>%
  left_join(sample_mapping %>% select(cell_line_label, cancer_type), by = 'cell_line_label')

n_bassik_hits_per_cell_line_plot <- plot_n_bassik_hits_per_cell_line(bassik_hits_per_cell_line, cancer_type_pal)
n_bassik_hits_per_cell_line_plot_path <- file.path(opt$out, 'number_of_bassik_hits_per_cell_line.png')
ggsave(filename = n_bassik_hits_per_cell_line_plot_path, plot = n_bassik_hits_per_cell_line_plot, device = 'png', dpi = 300 , width = 4000, height = 3000, units = 'px')

message(paste("Number of Bassik hits per cell line written to:", n_bassik_hits_per_cell_line_plot_path))

############################################################
# Number of Bassik hits per cell line                      #
############################################################

message('Plotting Bassik hits (number of cell lines by tumour type)...')

bassik_hits_mt5 <- bassik_binary_hits %>%
  filter(bassik__total >= 5) %>%
  select(sorted_gene_pair, 'Melanoma' = bassik__Melanoma, 'Pancreas' = bassik__Pancreas, 'Lung NSCLC' = bassik__Lung_NSCLC) %>%
  gather(cancer_type, n, -sorted_gene_pair)

stacked_bassik_hits_per_cell_line_plot <- plot_stacked_bassik_hits_per_cell_line(bassik_hits_mt5, cancer_type_pal)
stacked_bassik_hits_per_cell_line_plot_path <- file.path(opt$out, 'bassik_hits_cell_line_stacked_barplot.png')
ggsave(filename = stacked_bassik_hits_per_cell_line_plot_path, plot = stacked_bassik_hits_per_cell_line_plot, device = 'png', dpi = 300 , width = 3000, height = 4000, units = 'px')

message(paste("Stacked barplot of Bassik hits by cell line written to:", stacked_bassik_hits_per_cell_line_plot_path))

############################################################
# Bassik mean gi                                           #
############################################################

message('Plotting Bassik hits by cell line (mean normalised GI)...')

bassik_strict_hits <- bassik_binary_hits %>%
  filter(bassik__total >= 5) %>%
  select(sorted_gene_pair, bassik__total)

bassik_mean_gi <- bassik_hits %>%
  filter(sorted_gene_pair %in% bassik_strict_hits$sorted_gene_pair) %>%
  filter(is_bassik_hit == 1) %>%
  select(sorted_gene_pair, cell_line_label, mean_norm_gi) %>%
  left_join(bassik_strict_hits, by = 'sorted_gene_pair') %>%
  left_join(sample_mapping %>% select(cell_line_label, cancer_type) %>% unique, by = 'cell_line_label')

bassik_hits_mean_gi_plot <- plot_bassik_hits_mean_gi_per_cell_line(bassik_mean_gi, cancer_type_pal)
bassik_hits_mean_gi_plot_path <- file.path(opt$out, 'bassik_hits_mean_gi_per_cell_line.png')
ggsave(filename = bassik_hits_mean_gi_plot_path, plot = bassik_hits_mean_gi_plot, device = 'png', dpi = 300 , width = 3000, height = 4000, units = 'px')

message(paste("Mean normalised GI of Bassik hits by cell line written to:", bassik_hits_mean_gi_plot_path))

message('Done.')