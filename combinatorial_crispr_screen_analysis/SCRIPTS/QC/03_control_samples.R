suppressPackageStartupMessages(suppressWarnings(library(optparse)))
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))
suppressPackageStartupMessages(suppressWarnings(library(scales)))
suppressPackageStartupMessages(suppressWarnings(library(ggpubr)))
suppressPackageStartupMessages(suppressWarnings(library(ggsci)))
suppressPackageStartupMessages(suppressWarnings(library(GGally)))

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

# Cancer type palette
cancer_type_pal <- pal_npg(alpha = 0.3)(3)
names(cancer_type_pal) <- levels(sample_mapping$cancer_type)

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

# Subset control counts
control_counts <- count_matrix[,c(control_mean_index,control_sample_indices, 1:opt$annotation)]

# Gather control counts
control_counts.narrow <- control_counts %>%
  gather(sample_label, counts, -all_of(colnames(count_matrix)[1:opt$annotations])) %>%
  left_join(sample_mapping, by = 'sample_label')

# Sort control samples
control_counts.narrow$sample_label <- factor(control_counts.narrow$sample_label, levels = c('control_mean', control_sample_labels))

############################################################
# Correlation of control samples                           #
############################################################

message("Building pairwise correlation plot...")

# Pairwise correlation of controls
n_controls <- length(c(control_mean_index,control_sample_indices))
control_cor_path <- file.path(opt$out, 'control_correlation_normalised_counts.png')
control_counts_pairwise_correlation <- plot_control_correlation(control_counts, control_cor_path)

message(paste("Correlation plot saved to:", control_cor_path))

############################################################
# Count distribution of control samples                    #
############################################################

message("Building count distribution plot...")

count_dist_path <- file.path(opt$out, 'normalised_control_count_distributions.png')
control_counts_violin <- plot_control_violin(control_counts.narrow, cancer_type_pal)
ggsave(filename = count_dist_path, plot = control_counts_violin, device = 'png', dpi = 300 , width = 4000, height = 2000, units = 'px')

message(paste("Count distribution plot saved to:", count_dist_path))

############################################################
# Essential distribution of control samples                #
############################################################

message("Building control essentiality plot...")

# Set levels for sgrna_group
control_counts.narrow$sgrna_group <-
  factor(control_counts.narrow$sgrna_group,
         levels = c('Essential', 'Non-essential', 'Paralogue',
                    'Achilles CNV', 'CRISPR RNA', 'Safe-targeting control'))

# sgRNA group palette
sgrna_group_pal <- pal_nejm(alpha = 0.3)(6)
names(sgrna_group_pal) <- levels(control_counts.narrow$sgrna_group)

control_ess_path <- file.path(opt$out, 'normalised_control_essential_distribution.png')
control_ess_density <- plot_control_essentials(control_counts.narrow, sgrna_group_pal)
ggsave(filename = control_ess_path, plot = control_ess_density, device = 'png', dpi = 300 , width = 3800, height = 4000, units = 'px')

message(paste("Control essentiality plot saved to:", control_ess_path))

############################################################
# RDS                                                      #
############################################################

message("Saving R objects to file...")

control_samples.rds.path <- file.path(rds_path, 'control_samples.Rdata')
save(count_matrix, sample_mapping, control_sample_indices, control_mean_index, control_counts, control_counts_pairwise_correlation, control_counts_violin, control_ess_density, file = control_samples.rds.path )

message(paste("Saved objects to:", control_samples.rds.path))

message('Done.')
