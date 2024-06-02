suppressPackageStartupMessages(suppressWarnings(library(optparse)))
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))
suppressPackageStartupMessages(suppressWarnings(library(scales)))
suppressPackageStartupMessages(suppressWarnings(library(ggpubr)))
suppressPackageStartupMessages(suppressWarnings(library(ggsci)))
suppressPackageStartupMessages(suppressWarnings(library(GGally)))
suppressPackageStartupMessages(suppressWarnings(library(ggrepel)))

############################################################
# OPTIONS                                                  #
############################################################

option_list = list(
  make_option(c("-d", "--dir"), type = "character",
              help = "full path to repository", metavar = "character"),
  make_option(c("-c", "--counts"), type = "character",
              help = "full path to count matrix", metavar = "character"),
  make_option(c("-l", "--lfc"), type = "character",
              help = "full path to lfc matrix", metavar = "character"),
  make_option(c("-m", "--mapping"), type = "character",
              help = "full path to sample mapping", metavar = "character"),
  make_option(c("--helper"), type = "character",
              help = "full path to helper functions", metavar = "character"),
  make_option(c("--samples"), type = "character",
              help = "comma-delimited list of samples to remove", metavar = "character"),
  make_option(c("--annotations"), type = "integer",
              help = "number of annotation columns", metavar = "integer"),
  make_option(c("-r", "--rds"), type = "character",
              help = "full path to RDS directory", metavar = "character"),
  make_option(c("--plots_out"), type = "character",
              help = "full path to output directory", metavar = "character"),
  make_option(c("--counts_out"), type = "character",
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
check_dir_exists(opt$counts_out)
check_dir_exists(opt$plots_out)

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
cancer_type_pal <- pal_npg(alpha = 0.3)(5)
names(cancer_type_pal) <- c(levels(sample_mapping$cancer_type), 'Control', 'Summary')

############################################################
# Count matrix                                             #
############################################################

check_file_exists(opt$counts)

message(paste("Reading count matrix from:", opt$counts))

# Read in count matrix
count_matrix <- read.delim(opt$counts, sep = "\t", header = T, check.names = F)

############################################################
# LFC matrix                                               #
############################################################

check_file_exists(opt$lfc)

message(paste("Reading LFC matrix from:", opt$lfc))

# Read in LFC matrix
lfc_matrix <- read.delim(opt$lfc, sep = "\t", header = T, check.names = F)

############################################################
# Samples to remove                                        #
############################################################

message("Identifying user-defined samples to remove...")

# Get user-defined sample labels to remove
samples_to_remove <- scan(opt$samples, what = 'character', sep = "\n")

message(paste("Samples to remove:", paste(samples_to_remove, collapse = ', ')))

# Check samples are in sample mapping
for (es in samples_to_remove) {
  if (!es %in% sample_mapping$sample_label) {
    stop(print(paste('Sample not found in sample annotations:', es)))
  }
}

############################################################
# Remove samples from counts                               #
############################################################

message("Removing user-defined samples from counts...")

# Remove user-defined samples from count matrix
count_matrix.removed <- count_matrix %>%
  select(-all_of(samples_to_remove))

# Check column names have been removed from counts
if (ncol(count_matrix) - ncol(count_matrix.removed) != length(samples_to_remove)) {
  print(paste('Pre-filter column number:', ncol(count_matrix)))
  print(paste('Post-filter column number:', ncol(count_matrix.removed)))
  stop('Incorrect number of columns post-filtering.')
}

# Save updated count matrix
message("Saving count matrix with user-defined samples removed...")
excluded_sample_counts_path <- file.path(opt$counts_out, "count_matrix.norm.samples_removed.tsv")
write.table(count_matrix.removed, excluded_sample_counts_path, sep = "\t", row.names = F, quote = F)
message(paste("Count matrix with user-defined samples removed written to:", excluded_sample_counts_path))

############################################################
# Remove samples from fold changes                         #
############################################################

message("Removing user-defined samples from fold changes...")

# Remove user-defined samples from lfc matrix
lfc_matrix.removed <- lfc_matrix %>%
  select(-all_of(samples_to_remove))

# Check column names have been removed from fold changes
if (ncol(lfc_matrix) - ncol(lfc_matrix.removed) != length(samples_to_remove)) {
  print(paste('Pre-filter column number:', ncol(lfc_matrix)))
  print(paste('Post-filter column number:', ncol(lfc_matrix.removed)))
  stop('Incorrect number of columns post-filtering.')
}

# Save updated fold change  matrix
message("Saving LFC matrix with user-defined samples removed...")
excluded_sample_lfc_path <- file.path(opt$counts_out, "lfc_matrix.unscaled.samples_removed.tsv")
write.table(lfc_matrix.removed, excluded_sample_lfc_path, sep = "\t", row.names = F, quote = F)
message(paste("LFCmatrix with user-defined samples removed written to:", excluded_sample_lfc_path))

############################################################
# Sample correlations                                      #
############################################################

message("Calculating sample count correlations...")

# Get Spearman's correlation of controls
sample_cor <- cor(count_matrix.removed[(opt$annotation + 1):ncol(count_matrix.removed)], method = 'spearman')

############################################################
# Sample normalised count PCA                              #
############################################################

message("Principal component analysis...")

# PCA
pca <- prcomp(sample_cor)

# Get explained variance for plot axis titles
var_explained <- pca$sdev^2/sum(pca$sdev^2)

# Add sample mapping to principal components matrix
scores <- as.data.frame(pca$x) %>%
  rownames_to_column('sample_label') %>%
  filter(sample_label != 'control_mean') %>%
  left_join(sample_mapping, by = 'sample_label')

# Plot PCA with all samples
pca_plot <- plot_pca(scores, outliers = F)
pca_plot_path <- file.path(opt$plots_out, 'normalised_sample_count_pca.samples_removed.png')
ggsave(filename = pca_plot_path, plot = pca_plot, device = 'png', dpi = 300 , width = 4000, height = 4000, units = 'px')

message(paste("Sample PCA plot saved to:", pca_plot_path))


message('Done.')
