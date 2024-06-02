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
# Download CMP model metadata                       #
############################################################

message('Downloading CMP model metadata...')

# Download and read in CMP model metadata
tmpfile <- tempfile()
download.file(url = 'https://cog.sanger.ac.uk/cmp/download/model_list_20220628.csv', destfile = tmpfile, quiet = T)
cmp_models <- read.csv(tmpfile, check.names = F)

message('Filtering CMP models...')

cmp_models <- cmp_models %>% 
  rename('depMapID' = BROAD_ID) %>%
  filter(depMapID %in% unique(sample_mapping$depMapID)) %>%
  select(model_id, depMapID)

############################################################
# Download CMP mutations                                   #
############################################################

message('Downloading CMP mutations...')

# Download and read in CMP mutations 
tmpfile <- tempfile()
download.file(url = 'https://cog.sanger.ac.uk/cmp/download/mutations_summary_20220509.zip', destfile = tmpfile, quiet = T)
conn <- unz(tmpfile, 'mutations_summary_20220509.csv')
cmp_mutations <- read.csv(conn)
unlink(tmpfile)

############################################################
# Process CMP mutations                                    #
############################################################

message('Filtering CCLE mutations for screened cell lines...')

cmp_mutations.filt <- cmp_mutations %>%
  filter(model_id %in% cmp_models$model_id) %>% # Only need screened lines
  left_join(cmp_models, by = 'model_id')        # Add DepMap IDs

message('Filtering CCLE mutations for screened genes...')

cmp_mutations.filt <- cmp_mutations.filt %>%
  filter(gene_symbol %in% unique(expanded_library$targetA) | gene_symbol %in% unique(expanded_library$targetB))  # Only need screened genes

message('Filtering CCLE mutations for cancer drivers...')

cmp_mutations.filt <- cmp_mutations.filt %>%
  filter(cancer_driver == 'True') %>%
  mutate(cancer_driver = 'TRUE')

message('Selecting columns from CCLE mutations...')

cmp_mutations.ann <- sample_mapping %>%
  select(depMapID, cell_line_label, cancer_type) %>%
  unique() %>%
  left_join(cmp_mutations.filt, by = 'depMapID') %>%
  select(depMapID, 'CMP_model_id' = model_id, 
         cell_line_label, cancer_type,
         gene_symbol, protein_mutation, rna_mutation, cdna_mutation, effect, cancer_driver)

message('Removing NAs...')

cmp_mutations.ann <- cmp_mutations.ann %>%
  filter(!is.na(depMapID))

############################################################
# Write outputs to file                                    #
############################################################

message('Writing to file...')

# Save as a TSV
cmp_mutations.filepath <- file.path(opt$dir, 'DATA', 'postprocessing', 'intermediate_tables', 'cmp_cancer_driver_mutations.tsv')
write.table(cmp_mutations.ann, file = cmp_mutations.filepath, sep = "\t", row.names = F, quote = F)
message(paste0('CMP cancer driver mutations TSV written to:', cmp_mutations.filepath))

# Save as an Rdata object
cmp_mutations.rds.filepath <- file.path(opt$dir, 'DATA', 'RDS', 'postprocessing', 'cmp_cancer_driver_mutations.rds')
saveRDS(cmp_mutations.ann, file = cmp_mutations.rds.filepath)
message(paste0('CMP cancer driver mutations RDS written to:', cmp_mutations.rds.filepath))

message('Done.')