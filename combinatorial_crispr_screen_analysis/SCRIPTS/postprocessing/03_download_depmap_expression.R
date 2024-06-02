suppressPackageStartupMessages(suppressWarnings(library(optparse)))
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))

############################################################
# OPTIONS                                                  #
############################################################

option_list = list(
  make_option(c("-d", "--dir"), type = "character",
              help = "repository directory path", metavar = "character"),
  make_option(c("-m", "--mapping"), type = "character",
              help = "full path to sample mapping", metavar = "character"),
  make_option(c("-a", "--annotation"), type = "character",
              help = "full path to annotation", metavar = "character")
);

opt_parser <- OptionParser( option_list = option_list );
opt <- parse_args( opt_parser );

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
source(file.path(repo_path, 'SCRIPTS', 'postprocessing', 'helper.R'))

############################################################
# Sample mapping                                           #
############################################################

check_file_exists(opt$mapping)

message(paste("Reading sample annotations from:", opt$mapping))

# Read in sample mapping
sample_mapping <- read.delim(opt$mapping, header = T, sep = "\t")

############################################################
# Expanded library                                         #
############################################################

check_file_exists(opt$annotation)

message(paste("Reading expanded library:", opt$annotation))

# Read in library
expanded_library <- read.delim(opt$annotation, sep = "\t", header = T)

############################################################
# Download DepMap expression                               #
############################################################

message('Downloading DepMap gene expression...')

# Read in DepMap 22Q2 gene expression per cell line
# RNAseq TPM gene expression data for just protein coding genes using RSEM. Log2 transformed, using a pseudo-count of 1.  
tmpfile <- tempfile()
options(timeout=360)
download.file(url = 'https://ndownloader.figshare.com/files/34989919', destfile = tmpfile, quiet = T)
depmap_expression <- read.delim(tmpfile, sep = ",", header = T)
colnames(depmap_expression)[1] <- 'depMapID'

############################################################
# Filter expression for screened lines                     #
############################################################

message('Filtering to remove cell lines not screened...')

# Filter to get only the cell lines screened
depmap_expression.filt <- depmap_expression %>% 
  filter(depMapID %in% sample_mapping$depMapID)

message('Reshaping gene expression data frame...')

# Get genes per cell line by row
depmap_expression.filt <- depmap_expression.filt %>% 
  gather(gene_label, tpm_log2ps1, -depMapID) %>%
  separate(gene_label, sep = "\\.\\.", into = c('gene', 'entrez_id')) %>%
  mutate(entrez_id = str_replace_all(entrez_id, "\\.", ""))

############################################################
# Add to expanded library.                                 #
############################################################

message('Get gene pairs from library...')

# Reduce guide level library down to just gene pairs
expanded_library <- expanded_library %>%
  select(sorted_gene_pair, targetA, targetB) %>%
  unique()

############################################################
# Expand by library metadata and cell line                 #
############################################################

message('Expanding library by cell line...')

depmap_expression.ann <- expand.grid(expanded_library$sorted_gene_pair, sample_mapping$depMapID) %>%
  unique() %>%
  select('sorted_gene_pair' = 'Var1', 'depMapID' = 'Var2' ) %>%
  filter(!is.na(depMapID))

message('Adding additional library metadata...')

depmap_expression.ann <- depmap_expression.ann %>%
  left_join(expanded_library %>% select(sorted_gene_pair, targetA, targetB), by = 'sorted_gene_pair', relationship = 'many-to-many') %>%
  unique()

message('Adding copy number for targetA...')

depmap_expression.ann <- depmap_expression.ann %>%
  left_join(depmap_expression.filt %>% 
              select(-entrez_id) %>%
              rename_at(vars(-gene, -depMapID), ~ paste('targetA', ., sep = '__')), by = c('targetA' = 'gene', 'depMapID')) %>%
  unique()

message('Adding copy number for targetB...')

depmap_expression.ann <- depmap_expression.ann %>%
  left_join(depmap_expression.filt %>% 
              select(-entrez_id) %>%
              rename_at(vars(-gene, -depMapID), ~ paste('targetB', ., sep = '__')), by = c('targetB' = 'gene', 'depMapID')) %>%
  unique()

message('Adding sample mapping...')

depmap_expression.ann <- depmap_expression.ann %>%
  left_join(sample_mapping %>% select(cell_line_label, cancer_type, depMapID), by = 'depMapID', relationship = 'many-to-many') %>%
  unique()
              
message('Postprocessong...')

depmap_expression.ann <- depmap_expression.ann %>%
  select(depMapID, cell_line_label, cancer_type, 
         sorted_gene_pair, 
         targetA, targetB, 
         targetA__tpm_log2ps1, 
         targetB__tpm_log2ps1)

############################################################
# Write outputs to file                                    #
############################################################

message('Writing to file...')

# Save as a TSV
depmap_expression.filepath <- file.path(opt$dir, 'DATA', 'postprocessing', 'intermediate_tables', 'depmap_expression.tsv')
write.table(depmap_expression.ann, file = depmap_expression.filepath, sep = "\t", row.names = F, quote = F)
message(paste0('DepMap expression TSV written to:', depmap_expression.filepath))

# Save as an Rdata object
depmap_expression.rds.filepath <- file.path(opt$dir, 'DATA', 'RDS', 'postprocessing', 'depmap_expression.rds')
saveRDS(depmap_expression.ann, file = depmap_expression.rds.filepath)
message(paste0('DepMap expression RDS written to:', depmap_expression.rds.filepath))

message('Done.')