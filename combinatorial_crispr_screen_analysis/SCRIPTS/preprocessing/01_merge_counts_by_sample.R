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
# Library                                                  #
############################################################

check_file_exists(opt$annotation)

message(paste("Reading library:", opt$annotation))

# Read in library
lb <- read.delim(opt$annotation, sep = "\t", header = T)

############################################################
# Sample mapping                                           #
############################################################

check_file_exists(opt$mapping)

message(paste("Reading sample annotations from:", opt$mapping))

# Read in sample mapping
sample_mapping <- read.delim(opt$mapping, header = T, sep = "\t")

############################################################
# pyCROQUET counts per lane                                #
############################################################

# Get list of count file paths
pycroquet_dir <- file.path(repo_path, 'DATA', 'pyCROQUET')
pycroquet_count_files <- list.files(path = pycroquet_dir, pattern = '.counts.tsv', full.names = T)

message(paste("Reading pyCROQUET count files from:", pycroquet_dir))
message(paste("Number of count files:", length(pycroquet_count_files)))

# Read in count files. Organise by sample and lane. Takes a while to run.
pycroquet_counts_per_lane <- data.frame()

for (pcf in pycroquet_count_files) {
  rlt <- basename(pcf) %>% str_replace('.counts.tsv.gz', '')
  sn <- str_split(read_lines(pcf, n_max=1, skip = 2), '\t')[[1]][6]
  sn <- str_replace(sn, 'reads_', '')
  
  tmp_counts <- read.delim(pcf, header = F, sep = "\t", skip = 3, col.names = c('id', 'sgrna_ids', 'sgrna_seqs', 'gene_pair_id', 'unique_guide', 'counts'), colClasses = c(rep('character', 4), rep('integer', 2)))
  tmp_counts <- tmp_counts %>% mutate('sample_name' = sn, 'run_info' = rlt)
  
  if (nrow(pycroquet_counts_per_lane) == 0) {
    pycroquet_counts_per_lane <- tmp_counts
  } else {
    pycroquet_counts_per_lane <- bind_rows(pycroquet_counts_per_lane, tmp_counts)
  }
}

message(paste("Number of samples found:", length(unique(pycroquet_counts_per_lane$sample_name))))

############################################################
# pyCROQUET counts per sample                              #
############################################################

message("Calculating total counts per sample...")

# Get total counts per sample
pycroquet_counts <- pycroquet_counts_per_lane %>%
  separate(gene_pair_id, into = c('sgrnaA', 'sgrnaB'), sep = "\\|") %>%
  group_by(sample_name, sgrna_ids, sgrnaA, sgrnaB) %>%
  summarise('sample_counts' = sum(counts), .groups = 'keep') %>%
  ungroup()

message("Building count matrix...")

# Spread into sample count matrix
pycroquet_counts.wide <- sample_mapping %>% select(sanger_sample_name, sample_label) %>%
  left_join(pycroquet_counts, by = c('sanger_sample_name' = 'sample_name')) %>%
  select(-sanger_sample_name) %>%
  spread(sample_label, sample_counts, fill = 0) %>%
  filter(!is.na(sgrna_ids))

message("Expanding with redundant library...")
# Add annotations into count matrix
# Expands counts across non-redundant library
pycroquet_counts.wide.ann <- lb %>% 
   left_join(pycroquet_counts.wide, by = c('sgrna_ids', 'sgrnaA', 'sgrnaB'))

message(paste("Number of guides in count matrix:", nrow(pycroquet_counts.wide.ann)))

# Sort count matrix
pycroquet_counts.wide.ann <- pycroquet_counts.wide.ann %>%
  arrange(sorted_gene_pair, id)

############################################################
# Outputs                                                  #
############################################################

# Write raw counts per sample to output file
pycroquet_counts.path <- file.path(output_path, 'count_matrix.tsv')
write.table(pycroquet_counts.wide.ann, pycroquet_counts.path, sep = "\t", quote = F, row.names = F)

message(paste("Raw count matrix written to:", pycroquet_counts.path))

# Save as RDS
pycroquet_counts.rds.path <- file.path(rds_path, 'count_matrix.rds')
saveRDS(pycroquet_counts.wide.ann, file = pycroquet_counts.rds.path )

message(paste("Raw count matrix RDS written to:", pycroquet_counts.rds.path))

message('Done.')