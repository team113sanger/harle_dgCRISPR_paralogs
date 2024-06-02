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
  make_option(c("--annotations"), type = "integer",
              help = "number of annotation columns", metavar = "integer"),
  make_option(c("-m", "--mapping"), type = "character",
              help = "full path to sample mapping", metavar = "character"),
  make_option(c("--controls"), type = "character",
              help = "full path of file with control sample names", metavar = "character"),
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
# Sample correlations                                      #
############################################################

message("Calculating sample count correlations...")

# Get Spearman's correlation of controls
sample_cor <- cor(count_matrix[(opt$annotation + 1):ncol(count_matrix)], method = 'spearman')

message("Preparing sample count correlations...")

# Collapse control correlation
ind <- which(upper.tri(sample_cor, diag = F) , arr.ind = TRUE)
sample_cor.narrow <- data.frame(sample_label.x = dimnames(sample_cor)[[2]][ind[,2]],
                                 sample_label.y = dimnames(sample_cor)[[1]][ind[,1]],
                                 r = sample_cor[ind])

# Add cell line labels to data frame
sample_cor.narrow <- sample_cor.narrow %>%
  left_join(sample_mapping %>% select(sample_label, cell_line_label), by = c('sample_label.x' = 'sample_label'), relationship = "many-to-many") %>%
  rename('cell_line_label.x' = 'cell_line_label') %>%
  left_join(sample_mapping %>% select(sample_label, cell_line_label), by = c('sample_label.y' = 'sample_label'), relationship = "many-to-many") %>%
  rename('cell_line_label.y' = 'cell_line_label')

# Prepare correlation data for box plot
boxplot_cor_data <- prepare_correlation_data_for_boxplot(sample_cor.narrow)

message("Plotting sample count correlations...")

# Plot correlation data
sample_count_cor_boxplot <- plot_sample_correlation(boxplot_cor_data)
sample_count_cor_boxplot_path <- file.path(opt$out, 'normalised_sample_count_correlation_boxplot.png')
ggsave(filename = sample_count_cor_boxplot_path, plot = sample_count_cor_boxplot, device = 'png', dpi = 300 , width = 4000, height = 3200, units = 'px')

message(paste("Sample correlation plot saved to:", sample_count_cor_boxplot_path))

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
pca_plot <- plot_pca(scores)
pca_plot_path <- file.path(opt$out, 'normalised_sample_count_pca.png')
ggsave(filename = pca_plot_path, plot = pca_plot, device = 'png', dpi = 300 , width = 4000, height = 4000, units = 'px')

message(paste("Sample PCA plot saved to:", pca_plot_path))

############################################################
# RDS                                                      #
############################################################

message("Saving R objects to file...")

sample_correlation_and_pca.rds.path <- file.path(rds_path, 'sample_correlation_and_pca.Rdata')
save(count_matrix, sample_mapping, control_sample_indices, control_mean_index, sample_cor, sample_cor.narrow, boxplot_cor_data, sample_count_cor_boxplot, pca, scores, pca_plot, file = sample_correlation_and_pca.rds.path )

message(paste("Saved objects to:", sample_correlation_and_pca.rds.path))

message('Done.')
