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
  make_option(c("-p", "--pycroquet"), type = "character",
              help = "full path to pyCROQUET outputs", metavar = "character"),
  make_option(c("-m", "--mapping"), type = "character",
              help = "full path to sample mapping", metavar = "character"),
  make_option(c("--helper"), type = "character",
              help = "full path to helper functions", metavar = "character"),
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

# Check pyCROQUET outputs path exists
check_dir_exists(opt$pycroquet)

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
# Read in pyCROQUET statistics                             #
############################################################

message(paste("Reading pyCROQUET JSON statistics from:", opt$pycroquet))

# Read pyCROQUET JSON files into a data frame
pycroquet_statistics <- collate_pyc_stats_json(opt$pycroquet, sample_mapping)

############################################################
# Get sample mapping statistics                            #
############################################################

message("Calculating sample mapping statistics...")

# Combine pyCROQUET run statistics to get sample totals
pycroquet_sample_statistics <- get_pyc_sample_stats(pycroquet_statistics)
pycroquet_sample_statistics$sample_label <- factor(pycroquet_sample_statistics$sample_label, levels = sample_mapping$sample_label)

############################################################
# Plot read totals per pair.                               #
############################################################

message("Plotting total number of read pairs per sample...")

# Plot total read pairs per sample (coloured by cancer_type e.g. Melanoma, Lung, Pancreas)
total_pairs_barplot <- plot_total_read_pairs(pycroquet_sample_statistics)
total_pairs_barplot_path <- file.path(opt$out, 'total_read_pairs_per_sample.png')
ggsave(filename = total_pairs_barplot_path, plot = total_pairs_barplot, device = 'png', dpi = 300 , width = 4000, height = 2000, units = 'px')

message(paste("Total read pairs barplot written to:", total_pairs_barplot_path))

############################################################
# Plot mapping statistics                                  #
############################################################

message("Preparing mapping statistics...")

# Prepare mapping statistics
mapping_statistics <- prepare_mapping_statistics(pycroquet_sample_statistics, sample_mapping)

message("Plotting mapping rates of read pairs per sample...")

# Plot number of read pairs mapped
mapped_read_pairs <- plot_mapped_read_pair(mapping_statistics)
mapped_read_pairs_path <- file.path(opt$out, 'read_pairs_mapped_per_sample.png')
ggsave(filename = mapped_read_pairs_path, plot = mapped_read_pairs, device = 'png', dpi = 300 , width = 4000, height = 2000, units = 'px')

message(paste("Number of reads mapping barplot written to:", mapped_read_pairs_path))

# Plot proportion of read pairs mapped
prop_mapped_read_pairs <-plot_mapped_read_pair_proportion(mapping_statistics)
prop_mapped_read_pairs_path <- file.path(opt$out, 'proportion_of_read_pairs_mapped_per_sample.png')
ggsave(filename = prop_mapped_read_pairs_path, plot = prop_mapped_read_pairs, device = 'png', dpi = 300 , width = 4000, height = 2000, units = 'px')

message(paste("Proportion of reads mapping barplot written to:", prop_mapped_read_pairs_path))

############################################################
# RDS.                                                     #
############################################################

message("Saving R objects to file...")

pycroquet_statistics.rds.path <- file.path(rds_path, 'pycroquet_statistics.Rdata')
save(pycroquet_statistics, pycroquet_sample_statistics, mapping_statistics, total_pairs_barplot, mapped_read_pairs, prop_mapped_read_pairs, file = pycroquet_statistics.rds.path )

message(paste("Saved objects to:", pycroquet_statistics.rds.path))

message('Done.')
