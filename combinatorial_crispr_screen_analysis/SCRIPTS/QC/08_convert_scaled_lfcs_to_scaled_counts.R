suppressPackageStartupMessages(suppressWarnings(library(optparse)))
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))

############################################################
# OPTIONS                                                  #
############################################################

option_list = list(
  make_option(c("-d", "--dir"), type = "character",
              help = "full path to repository", metavar = "character"),
  make_option(c("-f", "--fc" ), type = "character",
               help = "scaled fold change matrix file", metavar = "character"),
  make_option(c("-c", "--counts" ), type = "character",
               help = "unscaled counts file", metavar = "character"),
  make_option(c("-s", "--suffix"), type = "character",
              help = "file name suffix", metavar = "character"),
  make_option(c("--annotations"), type = "integer",
              help = "number of annotation columns", metavar = "integer"),
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

# Check LFC matrix file exists 
check_file_exists(opt$fc)

# Check LFC matrix file exists 
check_file_exists(opt$fc)

# Check unscaled count matrix file exists 
check_file_exists(opt$counts)

# Determine file extensions with suffix
if (is.null(opt$suffix)) {
  tsv_ext <- 'tsv'
  rds_ext <- 'rds'
} else {
  tsv_ext <- paste0(opt$suffix, '.tsv')
  rds_ext <- paste0(opt$suffix, '.rds')
}

############################################################
# Scaled LFC matrix                                        #
############################################################

# Set LFC matrix path
lfc_matrix_file <- opt$fc

message(paste("Reading LFC matrix:", lfc_matrix_file))

# Check LFC matrix file exists
if (!file.exists(lfc_matrix_file)) {
  stop(sprintf("LFC matrix file does not exist: %s", lfc_matrix_file))
}

# Read in LFC matrix
lfc_matrix <- read.delim(lfc_matrix_file, sep = "\t", header = T, check.names = F)

############################################################
# Unscaled count matrix                                    #
############################################################

# Set count matrix path
count_matrix_file <- opt$counts

message(paste("Reading count matrix:", count_matrix_file))

# Check count matrix file exists
if (!file.exists(count_matrix_file)) {
  stop(sprintf("Count matrix file does not exist: %s", count_matrix_file))
}

# Read in count matrix
count_matrix <- read.delim(count_matrix_file, sep = "\t", header = T, check.names = F)
  
# Set annotation column names
annotation_colnames <- colnames(count_matrix)[1:opt$annotations]

############################################################
# Gather matrices                                          #
############################################################

message("Gathering scaled LFC matrix...")
lfc_matrix.narrow <- lfc_matrix %>% 
  gather(sample_label, scaled_lfc, -all_of(annotation_colnames))

message("Gathering unscaled count matrix...")
count_matrix.narrow <- count_matrix %>%
  gather(sample_label, unscaled_counts, -all_of(annotation_colnames))

############################################################
# Get unscaled control counts                              #
############################################################

message("Getting unscaled control counts...")
unscaled_control_counts <- count_matrix.narrow %>%
  filter(sample_label == "control_mean") %>% # Get only control counts
  rename('unscaled_control_counts' = 'unscaled_counts') %>% # Rename counts column
  select(-sample_label) 

############################################################
# Remove unwanted sample_labels                            #
############################################################

message("Removing unwanted samples...")
count_matrix.narrow <- count_matrix.narrow %>%
  filter(sample_label %in% lfc_matrix.narrow$sample_label)

############################################################
# Process matrices                                         #
############################################################

# Add unscaled control counts to gathered unscaled count matrix
message("Adding unscaled control counts to unscaled count matrix...")
count_matrix.narrow <- count_matrix.narrow %>%
  left_join(unscaled_control_counts, by = c(annotation_colnames))

# Combine unscaled counts and scaled LFCs into new dataframe
message("Combining unscaled counts and scaled log fold changes...")
combined_matrix <- count_matrix.narrow %>%
  left_join(lfc_matrix.narrow, by = c('sample_label', annotation_colnames))

############################################################
# Get scaled counts                                        #
############################################################

# Revert scaled LFCs to scaled counts
message("Reverting scaled log fold changes to scaled counts...")
combined_matrix <- combined_matrix %>%
  mutate('scaled_counts' = unscaled_control_counts * (2 ^ scaled_lfc))

# Spread scaled count matrix 
message("Generating scaled count matrix...")
scaled_count_matrix <- combined_matrix %>%
  select(all_of(annotation_colnames), sample_label, scaled_counts) %>%
  spread(sample_label, scaled_counts)

# Add unscaled control counts into scaled count matrix
message("Adding control_mean into scaled count matrix...")
scaled_count_matrix <- scaled_count_matrix %>%
  left_join(count_matrix %>% select(all_of(annotation_colnames), control_mean), 
            by = c(annotation_colnames))

# Check that the new matrix is the correct length
if (nrow(count_matrix) != nrow(scaled_count_matrix)) {
  message(paste('Number of rows in count matrix pre-scaling:', nrow(count_matrix)))
  message(paste('Number of rows in count matrix post-scaling:', nrow(scaled_count_matrix)))
  stop("Incorrect number of rows in count matrix after scaling.")
}

############################################################
# Outputs                                                  #
############################################################

# Write scaled counts to RDS
scaled_count_matrix.rds.path <- file.path(rds_path, paste0('count_matrix.scaled.', rds_ext))
saveRDS(scaled_count_matrix, file = scaled_count_matrix.rds.path)
message(paste("Scaled count matrix RDS written to:", scaled_count_matrix.rds.path))

# Write scaled counts to TSV
scaled_count_matrix.path <- file.path(output_path, paste0('count_matrix.scaled.', tsv_ext))
write.table(scaled_count_matrix, scaled_count_matrix.path, sep = "\t", quote = F, row.names = F)
message(paste("Scaled count matrix TSV written to:", scaled_count_matrix.path))

message('Done.')