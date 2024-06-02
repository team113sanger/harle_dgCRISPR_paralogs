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
  make_option(c("--annotations"), type = "integer",
              help = "number of annotation columns", metavar = "integer"),
  make_option(c("-p", "--pseudocount"), type = "integer", default = 5,
              help = "pseudocount", metavar = "integer"),
  make_option(c("-s", "--scalingfactor"), type = "integer", default = 10000000,
              help = "scaling factor", metavar = "integer"),
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
  stop(smessagef("Repository directory not exist: %s", repo_path))
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
# Normalise using BAGEL method                             #
############################################################

message("Adding pseudocount...")

# Add a user-defined pseudocount 
count_matrix.pseudo <- count_matrix %>% mutate(across(c((opt$annotations + 1):ncol(count_matrix)), ~ . + opt$pseudocount))

message("Total normalisation with scaling factor...")

# Total normalisation with user-defined scaling factor 
count_matrix.norm <- count_matrix.pseudo %>%
  mutate(across((opt$annotations + 1):ncol(count_matrix.pseudo), ~ (. / sum(.)) * opt$scalingfactor))

# Check same number of guides as we started with
if (nrow(count_matrix) != nrow(count_matrix.norm)) {
  message(paste('Number of rows pre-normalisation:', nrow(count_matrix)))
  message(paste('Number of rows post-normalisation:', nrow(count_matrix.norm)))
  stop('Number of rows is not consistent.')
}

############################################################
# Outputs                                                  #
############################################################

# Write normalised counts to file
count_matrix.norm.path <- file.path(output_path, 'count_matrix.norm.tsv')
write.table(count_matrix.norm, count_matrix.norm.path, sep = "\t", quote = F, row.names = F)

message(paste("Normalised count matrix written to:", count_matrix.norm.path))

# Save as RDS
count_matrix.norm.rds.path <- file.path(rds_path, 'count_matrix.norm.rds')
saveRDS(count_matrix.norm, file = count_matrix.norm.rds.path )

message(paste("Normalised count matrix RDS written to:", count_matrix.norm.rds.path))

message("DONE.")