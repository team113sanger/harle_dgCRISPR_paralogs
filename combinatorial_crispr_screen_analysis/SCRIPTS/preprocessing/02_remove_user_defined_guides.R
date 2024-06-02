suppressPackageStartupMessages(suppressWarnings(library(optparse)))
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))

############################################################
# OPTIONS                                                  #
############################################################

option_list = list(
  make_option(c("-d", "--dir"), type = "character",
              help = "full path to repository", metavar = "character"),
  make_option(c("-g", "--guide"), type = "character",
              help = "full path to file with list of guide ids to remove", metavar = "character"),
  make_option(c("-c", "--counts"), type = "character",
              help = "full path to count matrix", metavar = "character"),
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
# User-defined guides                                      #
############################################################

# Check file with guide ids exists
check_file_exists(opt$guide)

# Read in guide ids to remove
guides_to_remove <- scan(opt$guide, what = character(), strip.white = T, quiet = T)

# Check we have one or more ids
if (is.null(guides_to_remove) || length(guides_to_remove) == 0) {
  print(paste('Error reading guides to remove:', opt$guide))
  stop()
}

############################################################
# Count matrix                                             #
############################################################

check_file_exists(opt$counts)

print(paste("Reading count matrix:", opt$counts))

# Read in count matrix
count_matrix <- read.delim(opt$counts, sep = "\t", header = T, check.names = F)

############################################################
# User-defined guide pairs to remove                       #
############################################################

print("Removing user defined guides.")

count_matrix.guides_removed <- count_matrix %>%
  filter(!id %in% guides_to_remove)

print(paste('Number of guides removed:', length(guides_to_remove)))

# Check that the new matrix is the correct length
if (nrow(count_matrix) - nrow(count_matrix.guides_removed) != length(guides_to_remove)) {
  print(paste('Number of rows pre-removal:', nrow(count_matrix)))
  print(paste('Number of rows post-removal:', nrow(count_matrix.guides_removed)))
  stop("Incorrect number of rows after removing guides.")
}

############################################################
# Outputs                                                  #
############################################################

# Write count matrix with user-defined guides removed to file
count_matrix.guides_removed.path <- file.path(output_path, 'count_matrix.guides_removed.tsv')
write.table(count_matrix.guides_removed, count_matrix.guides_removed.path, sep = "\t", quote = F, row.names = F)

print(paste("Count matrix with user-defined guides removed written to:", count_matrix.guides_removed.path))

# Save as RDS
count_matrix.guides_removed.rds.path <- file.path(rds_path, 'count_matrix.guides_removed.rds')
saveRDS(count_matrix.guides_removed, file = count_matrix.guides_removed.rds.path )

print(paste("Count matrix with user-defined guides removed RDS written to:", count_matrix.guides_removed.rds.path))

print('DONE.')