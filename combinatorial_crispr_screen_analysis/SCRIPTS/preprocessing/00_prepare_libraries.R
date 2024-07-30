suppressPackageStartupMessages(suppressWarnings(library(optparse)))
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))
suppressPackageStartupMessages(suppressWarnings(library(stringi)))

############################################################
# OPTIONS                                                  #
############################################################

option_list = list(
  make_option(c("-d", "--dir"), type = "character",
              help = "full path to repository", metavar = "character"),
  make_option(c("-m", "--mapping"), type = "character",
              help = "full path to sample mapping", metavar = "character"),
  make_option(c("-a", "--annotation"), type = "character",
              help = "full path to annotation", metavar = "character"),
  make_option(c("--helper"), type = "character",
              help = "full path to helper functions", metavar = "character"),
  make_option(c("-r", "--rds"), type = "character",
              help = "full path to RDS directory", metavar = "character"))

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

############################################################
# Sample mapping                                           #
############################################################

check_file_exists(opt$mapping)

message(paste("Reading sample annotations from:", opt$mapping))

# Read in sample mapping
sample_mapping <- read.delim(opt$mapping, header = T, sep = "\t")

############################################################
# Library annotation                                       #
############################################################

check_file_exists(opt$annotation)

message(paste("Reading library annotations:", opt$annotation))

# Read in library
library_annotation <- read.delim(opt$annotation, sep = "\t", header = T)

# Save library as RDS
library_annotation.rds.path <- file.path(rds_path, 'paralog_library.rds')
saveRDS(library_annotation, file = library_annotation.rds.path )
message(paste("Library RDS written to:", library_annotation.rds.path))

############################################################
# pyCROQUET library                                        #
############################################################

message("Formatting library to use with pyCROQUET...")

# Set path for formatted library
pyc_library.path <- file.path(repo_path, 'METADATA', 'libraries', 'paralog_library.pycroquet.tsv')

# Write header lines
write(paste0("##library-type: dual"), file = pyc_library.path)
write(paste0("##library-name: dgCRISPR_paralogs"), file = pyc_library.path, append = TRUE)
write(paste0("##dual-orientation: R2_R1"), file = pyc_library.path, append = TRUE)
write(paste("#id", "sgrna_ids", "sgrna_seqs", "gene_pair_id", sep = "\t"), file = pyc_library.path, append = TRUE)

# Format library for quantification with pyCROQUET
# Note: library needs to be non-redundant for the sum of the pair classifications to be the same as the total pairs processed
pyc_library <- library_annotation %>% 
  mutate('gene_pair_id' = paste(sgrnaA, sgrnaB, sep = '|')) %>%
  select(sgrna_ids, sgrna_seqs, gene_pair_id) %>%
  unique()

# Loop over library dataframe 
for (i in 1:nrow(pyc_library)) {
  write(paste(pyc_library$sgrna_ids[i], 
              pyc_library$sgrna_ids[i], 
              pyc_library$sgrna_seqs[i], 
              pyc_library$gene_pair_id[i], 
              sep = "\t" ), file = pyc_library.path, append = TRUE)
}

message(paste("pyCROQUET-formatted library written to:", pyc_library.path))

# Save as RDS
pyc_library.rds.path <- file.path(rds_path, 'dual_guide_matrix.rds')
saveRDS(pyc_library, file = pyc_library.rds.path )

message(paste("pyCROQUET library RDS written to:", pyc_library.rds.path))

############################################################
# Build dual guide matrix                                  #
############################################################

# For each dual (gene|gene) guide we need to know the corresponding single guides (e.g. safe|gene and gene|safe)
# Function is in helper.R
dual_guide_matrix <- get_guide_matrix(library_annotation, 
                                      control_guide_file = file.path(repo_path, 'METADATA', 'safe_targeting_guides.txt'))

# Write guide matrix to output file
dual_guide_matrix.path <- file.path(repo_path, 'METADATA', 'libraries', 'dual_guide_matrix.tsv')
write.table(dual_guide_matrix, dual_guide_matrix.path, sep = "\t", quote = F, row.names = F)

message(paste("Guide matrix written to:", dual_guide_matrix.path))

# Save as RDS
dual_guide_matrix.rds.path <- file.path(rds_path, 'dual_guide_matrix.rds')
saveRDS(dual_guide_matrix, file = dual_guide_matrix.rds.path )

message(paste("Guide matrix RDS written to:", dual_guide_matrix.rds.path))

message('Done.')