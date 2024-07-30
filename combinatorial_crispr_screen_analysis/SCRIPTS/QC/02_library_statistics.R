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
  make_option(c("-c", "--counts"), type = "character",
              help = "full path to count matrix", metavar = "character"),
  make_option(c("-m", "--mapping"), type = "character",
              help = "full path to sample mapping", metavar = "character"),
  make_option(c("--helper"), type = "character",
              help = "full path to helper functions", metavar = "character"),
  make_option(c("--controls"), type = "character",
              help = "full path of file with control sample names", metavar = "character"),
  make_option(c("--low"), type = "integer", default = 20,
              help = "threshold for low counts", metavar = "integer"),
  make_option(c("--annotations"), type = "integer",
              help = "number of annotation columns", metavar = "integer"),
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
helper_path <- ifelse(is.null(opt$helper), file.path(repo_path, 'SCRIPTS', 'QC', 'helper.R'), opt$helper)
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
# Sample mapping                                           #
############################################################

check_file_exists(opt$mapping)

message(paste("Reading sample annotations from:", opt$mapping))

# Read in sample mapping
sample_mapping <- read_sample_metadata(opt$mapping)

############################################################
# Count matrix                                             #
############################################################

check_file_exists(opt$counts)

message(paste("Reading count matrix from:", opt$counts))

# Read in count matrix
count_matrix <- read.delim(opt$counts, sep = "\t", header = T, check.names = F)


############################################################
# Control samples                                          #
############################################################

message("Identifying control samples...")

# Get user-defined control sample labels
control_sample_labels <- scan(opt$controls, what = 'character', sep = "\n")

# Check control samples are in sample mapping
for (cs in control_sample_labels) {
  if (!cs %in% sample_mapping$sample_label) {
    stop(print(paste('Control sample not found in sample annotations:', cs)))
  }
}

# Identify indices of control samples (except control_mean)
control_sample_indices <- which(colnames(count_matrix) %in% setdiff(control_sample_labels, c('control_mean')))

# Identify index of control_mean
control_mean_index <- which(colnames(count_matrix) == 'control_mean')

############################################################
# Calculate library statistics                             #
############################################################

message("Gathering count matrix...")

# Gather count matrix
count_matrix.narrow <- count_matrix %>%
  gather(sample_label, norm_count, -all_of(colnames(count_matrix)[1:opt$annotations]))

message("Calculating library statistics...")

# Calculate library statistics
library_statistics <- get_library_statistics(count_matrix.narrow)

message("Adding sample metadata to library statistics...")

# Add sample annotations to library statistics
library_statistics <- library_statistics %>%
  left_join(sample_mapping, by = 'sample_label')

############################################################
# Plot library statistics                                  #
############################################################

message("Plotting median counts per sample...")

# Plot median counts per sample as bar plot
median_count_barplot <- plot_library_statistic(library_statistics %>% filter(sample_label != 'control_mean'), 'median', 'Number of reads')
median_count_barplot_path <- file.path(opt$out, 'median_normalised_counts_per_sample.png')
ggsave(filename = median_count_barplot_path, plot = median_count_barplot, device = 'png', dpi = 300 , width = 4000, height = 2000, units = 'px')

message(paste("Median counts per sample barplot written to:", median_count_barplot_path))

message("Plotting low counts per sample...")

# Plot low counts per sample as bar plot
low_count_barplot <- plot_library_statistic(library_statistics %>% filter(sample_label != 'control_mean'), 'low_counts', 'Number of guides')
low_count_barplot_path <- file.path(opt$out, 'low_normalised_counts_per_sample.png')
ggsave(filename = low_count_barplot_path, plot = low_count_barplot, device = 'png', dpi = 300 , width = 4000, height = 2000, units = 'px')

message(paste("Low counts per sample barplot written to:", low_count_barplot_path))

message("Plotting Gini index per sample...")

# Plot Gini index per sample as bar plot
gini_index_barplot <- plot_library_statistic(library_statistics %>% filter(sample_label != 'control_mean'), 'gini_index', 'Gini index')
gini_index_barplot_path <- file.path(opt$out, 'gini_index_per_sample.png')
ggsave(filename = gini_index_barplot_path, plot = gini_index_barplot, device = 'png', dpi = 300 , width = 4000, height = 2000, units = 'px')

message(paste("Gini index per sample barplot written to:", gini_index_barplot_path))

############################################################
# Write library statistics                                 #
############################################################

message("Writing library statistics...")

lib_stats_path <- file.path(opt$out, 'library_statistics.tsv')
write.table(library_statistics, lib_stats_path, row.names = F, quote = F, sep = "\t")

message(paste("Library statistics written to:", lib_stats_path))

############################################################
# RDS                                                      #
############################################################

message("Saving R objects to file...")

library_statistics.rds.path <- file.path(rds_path, 'library_statistics.Rdata')
save(count_matrix.narrow, control_sample_indices, control_mean_index, library_statistics, median_count_barplot,low_count_barplot, gini_index_barplot,  file = library_statistics.rds.path )

message(paste("Saved objects to:", library_statistics.rds.path))

message('Done.')
