suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))

############################################################
# OPTIONS                                                  #
############################################################

option_list = list(
  make_option(c("-d", "--dir"), type = "character",
              help = "full path to repository", metavar = "character"),
  make_option(c("-f", "--fc"), type = "character",
              help = "full path to fold change matrix file", metavar = "character"),
  make_option(c("--doubles_guide_matrix"), type = "character", default = 'norm',
              help = "full path to doubles guide matrix file", metavar = "character"),
  make_option(c("-i", "--chunk_index"), type="integer", default = NULL,
              help="which dual matrix chunk to process", metavar="integer"),
  make_option(c("-n", "--chunk_size"), type="integer", default=NULL,
              help="dual matrix chunk size", metavar="integer"),
  make_option(c("--annotations"), type = "integer",
              help = "number of annotation columns", metavar = "integer"),
  make_option(c("-o", "--out"), type = "character", default = '.',
              help = "output directory [Default: . ]", metavar = "character"),
  make_option(c("--helper"), type = "character",
              help = "full path to helper functions", metavar = "character")
);

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

############################################################
# VALIDATION                                               #
############################################################

# Set top level directory
repo_path <- opt$dir

# Fold change matrix file
if (is.null(opt$fc)) {
  print_help(opt_parser)
  stop("Please provide a fold change matrix file", call.=FALSE)
}

if (!file.exists(opt$fc)) {
  print_help(opt_parser)
  stop(paste("Fold change matrix file does not exist:", opt$fc), call.=FALSE)
}

# Doubles guide matrix file
if (is.null(opt$doubles_guide_matrix)) {
  print_help(opt_parser)
  stop("Please provide a doubles guide matrix file", call.=FALSE)
}

if (!file.exists(opt$doubles_guide_matrix)) {
  print_help(opt_parser)
  stop(paste("Doubles guide matrix file does not exist:", opt$doubles_guide_matrix), call.=FALSE)
}

############################################################
# FUNCTIONS                                                #
############################################################

# Add helper functions
helper_path <- ifelse(is.null(opt$helper), file.path(repo_path, 'SCRIPTS', 'dual_guide', 'helper.R'), opt$helper)
if (!file.exists(helper_path)){
  stop(paste('Helper file does not exist:', helper_path))
}
source(helper_path)

############################################################
# MAIN SCRIPT                                              #
############################################################

# Read in dual guide matrix (maps dual guide to its single components)
message(paste("Reading in dual guide matrix:", opt$doubles_guide_matrix))
doubles_guide_matrix <- read.delim(file = opt$doubles_guide_matrix , sep = "\t", header = T, check.names = F)

# Read in fold change matrix
message(paste("Reading in fold change matrix:", opt$fc))
fc <- read.delim(file = opt$fc, sep = "\t", header = T, check.names = F)

# Set annotation column names
message('Setting annotation column names...')
annotation_colnames <- colnames(fc)[1:opt$annotations]

# Narrow the FC matrix
message('Gathering fold change matrix...')
fc.narrow <- fc %>% gather(sample, fc, -all_of(annotation_colnames))

# Get sample names
message('Getting sample names...')
samples <- unique(fc.narrow$sample)
message(paste('Number of samples found:', length(samples)))

# Get chunks for double guide matrix
message('Getting chunks from dual guide matrix...')
chunk_list <- split(c(1:nrow(doubles_guide_matrix)), ceiling(seq_along(c(1:nrow(doubles_guide_matrix))) / opt$chunk_size))
message(paste('Number of chunks:', length(names(chunk_list))))

# Getting user-defined chunk
doubles_chunk_indexes <- unlist(chunk_list[opt$chunk_index])
chunk_doubles_guide_matrix <- doubles_guide_matrix[doubles_chunk_indexes, ]
message(paste('Chunk selected for processing:', opt$chunk_index))

# Get observed and predicted fold changes for doubles (per sample lists)
message('Building observed and predicted fold changes...')
results <- get_pred_vs_obs_y12_for_all_samples(samples, chunk_doubles_guide_matrix, fc.narrow)

# Bring together results into single data frame
message('Collating results...')

pred_vs_obs_y12 <- data.frame()
missing_data <- data.frame()

for (sn in samples) {
  print(paste("Merging results for:", sn))
  if (nrow(pred_vs_obs_y12) == 0) {
    pred_vs_obs_y12 <- results[['pred_vs_obs_y12']][[sn]]
    missing_data <- results[['missing_data']][[sn]]
  } else {
    pred_vs_obs_y12 <- rbind(pred_vs_obs_y12, as.data.frame(results[['pred_vs_obs_y12']][[sn]]))
    missing_data <- rbind(missing_data, results[['missing_data']][[sn]])
  }
}

# Save pred vs obs for chunk
message('Saving pred_vs_obs_y12...')
pred_vs_obs_y12_filename = file.path(opt$out, paste("pred_vs_obs_y12", opt$chunk_index, "tsv", sep = "."))
write.table(pred_vs_obs_y12, file = pred_vs_obs_y12_filename, row.names = F, sep = "\t", quote = F)

# Save missing data for chunk
message('Saving missing data...')
missing_data_filename = file.path(opt$out, paste( "missing_data", opt$chunk_index, "tsv", sep = "."))
write.table(missing_data, file = missing_data_filename, row.names = F, sep = "\t", quote = F)

message('Done.')