suppressPackageStartupMessages(suppressWarnings(library(optparse)))
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))

############################################################
# OPTIONS                                                  #
############################################################

option_list = list(
  make_option(c("-d", "--dir"), type = "character",
              help = "full path to repository", metavar = "character"),
  make_option(c("-m", "--mapping"), type = "character",
              help = "full path to sample mapping", metavar = "character"),
  make_option(c("-c", "--counts"), type = "character",
              help = "full path to count matrix", metavar = "character"),
  make_option(c("--annotations"), type = "integer",
              help = "number of annotation columns", metavar = "integer"),
  make_option(c("--samples"), type = "character",
              help = "comma-delimited list of control sample names", metavar = "character"),
  make_option(c("--helper"), type = "character",
              help = "full path to helper functions", metavar = "character"),
  make_option(c("-r", "--rds"), type = "character",
              help = "full path to RDS directory", metavar = "character"),
  make_option(c("-o", "--out"), type = "character",
              help = "full path to output directory", metavar = "character")
);

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

############################################################
# GENERAL                                                  #
############################################################

# Set top level directory
repo_path <- opt$dir

# Check top level directory exists
if (!dir.exists(repo_path)) {
  stop(sprintf("Repository directory not exist: %s", repo_path))
}

# Add helper functions
helper_path <- ifelse(is.null(opt$helper), file.path(repo_path, 'SCRIPTS', 'preprocessing', 'helper.R'), opt$helper)
if (!file.exists(helper_path)){
  stop(paste('Helper file does not exist:', helper_path))
}
source(helper_path)

# Set RDS path
rds_path <- ifelse(is.null(opt$rds), file.path(repo_path, 'DATA', 'RDS', 'preprocessing'), opt$rds)
check_dir_exists(rds_path)

# Set output path
output_path <- ifelse(is.null(opt$output), file.path(repo_path, 'DATA', 'preprocessing'), opt$output)
check_dir_exists(output_path)

############################################################
# Count matrix                                             #
############################################################

check_file_exists(opt$counts)

message(paste("Reading count matrix:", opt$counts))

# Read in count matrix
count_matrix <- read.delim(opt$counts, sep = "\t", header = T, check.names = F)

############################################################
# Sample mapping                                           #
############################################################

check_file_exists(opt$mapping)

message(paste("Reading sample annotations from:", opt$mapping))

# Read in sample mapping
sample_mapping <- read.delim(opt$mapping, header = T, sep = "\t")

############################################################
# Identify control samples                                 #
############################################################

message("Identifying control samples...")

# Get user-defined control sample labels
control_sample_labels <-  str_split(opt$samples, ',', simplify = T)

# Check control samples are in sample mapping
for (cs in control_sample_labels) {
  if (!cs %in% sample_mapping$sample_label) {
    stop(print(paste('Control sample not found in sample annotations:', cs)))
  }
}

# Identify indices of control samples
control_sample_indices <- which(colnames(count_matrix) %in% control_sample_labels)

# Check that there are the right number of control samples
if (length(control_sample_labels) != length(control_sample_indices)) {
  stop(print(paste('Number of control sample labels does not match number of control sample indices:', 
                  length(control_sample_labels), length(control_sample_indices))))
}

############################################################
# Calculate control mean                                   #
############################################################

message('Calculating control mean...')

# Get control mean
count_matrix.mean_control <- count_matrix %>%
  mutate('control_mean' = rowMeans(count_matrix[,control_sample_indices])) 

# Write filtered counts with mean control to file
count_matrix.mean_control.path <- file.path(output_path, 'count_matrix.filt.control_mean.tsv')
write.table(count_matrix.mean_control, count_matrix.mean_control.path, sep = "\t", quote = F, row.names = F)

message(paste("Count matrix with control mean written to:", count_matrix.mean_control.path))

############################################################
# Calculate log fold changes with respect to control mean  #
############################################################

message("Calculating unscaled LFCs...")

# Set annotation column names
annotation_colnames <- colnames(count_matrix)[1:opt$annotations]

# Count matrix with individual controls removed
count_matrix.mean_control <- count_matrix.mean_control %>%
  select(everything(), -contains('control'), control_mean)

# Calculate log fold changes
lfc_matrix <- suppressMessages(calculate_lfc(count_matrix.mean_control,
                            control_indices = which(colnames(count_matrix.mean_control) == 'control_mean'),
                            treatment_indices = which(! colnames(count_matrix.mean_control) %in% c(annotation_colnames, 'control_mean')),
                            pseudocount =  0.5))

# Add annotations back into matrix
lfc_matrix.ann <- lfc_matrix %>%
 left_join(count_matrix.mean_control[,1:opt$annotations], by = c('id', 'sgrna_ids')) %>%
 relocate(annotation_colnames[3:length(annotation_colnames)], .after = 'sgrna_ids')

# Check number of rows in count and LFC matrix are the same
if (nrow(count_matrix) != nrow(lfc_matrix.ann)) {
  message(paste('Number of count matrix rows:', nrow(count_matrix)))
  message(paste('Number of LFC matrix rows:', nrow(lfc_matrix.ann)))
  stop('Number of rows is not consistent.')
}

############################################################
# Outputs                                                  #
############################################################

# Write normalised/filtered log2 fold changes to file
lfc_matrix.path <- file.path(output_path, 'lfc_matrix.unscaled.tsv')
write.table(lfc_matrix.ann, lfc_matrix.path, sep = "\t", quote = F, row.names = F)

message(paste("Unscaled LFC matrix written to:", lfc_matrix.path))

# Save as RDS
lfc_matrix.rds.path <- file.path(rds_path, 'lfc_matrix.unscaled.rds')
saveRDS(lfc_matrix.ann, file = lfc_matrix.rds.path )

message(paste("Unscaled LFC matrix RDS written to:", lfc_matrix.rds.path ))

message("DONE.")