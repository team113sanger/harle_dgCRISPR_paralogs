suppressPackageStartupMessages(suppressWarnings(library(optparse)))
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))

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
  make_option(c("--filter"), type = "integer", default = 0,
              help = "minimum value for filter", metavar = "integer"),
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
# Filter low count guides (in controls)                    #
############################################################

message('Identifying low count guides...')

# Identify low count guides
filt_guides <- get_guides_failing_filter(count_matrix,
                                         id_column = 1,
                                         count_column = 5:ncol(count_matrix),
                                         filter_indices = control_sample_indices,
                                         filter_method = 'mean',
                                         min_reads = opt$filter)

message(paste('Number of low count guides ids:', length(filt_guides)))

message('Removing low count guides...')

# Remove low count guides from normalised counts (76 guides)
count_matrix.filt <- count_matrix %>%
  filter(!id %in% filt_guides)

# Check that the correct number of guides have been removed
num_guide_to_remove <- length(filt_guides)
message(paste('Count matrix rows:', nrow(count_matrix)))
message(paste('Filtered count matrix rows:', nrow(count_matrix.filt)))
if (nrow(count_matrix) - nrow(count_matrix.filt) != num_guide_to_remove) {
  message('Unexpected number of guides removed')
}

############################################################
# Outputs                                                  #
############################################################

# Write filtered guides to file
filt_guides.path <- file.path(output_path, 'filtered_guides.txt')
write.table(filt_guides, filt_guides.path, sep = "\t", quote = F, row.names = F, col.names = F)

message(paste("Filtered guide ids written to:", filt_guides.path))

# Write filtered counts to file
count_matrix.filt.path <- file.path(output_path, 'count_matrix.filt.tsv')
write.table(count_matrix.filt, count_matrix.filt.path, sep = "\t", quote = F, row.names = F)

message(paste("Filtered count matrix written to:", count_matrix.filt.path))

# Save as RDS
count_matrix.filt.rds.path <- file.path(rds_path, 'count_matrix.filt.rds')
saveRDS(count_matrix.filt, file = count_matrix.filt.rds.path )

message(paste("Filtered count matrix RDS written to:", count_matrix.filt.rds.path))

message("DONE.")