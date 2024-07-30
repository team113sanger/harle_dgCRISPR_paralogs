suppressPackageStartupMessages(suppressWarnings(library(optparse)))
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))

options(timeout=200)

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
# Download DepMap depletion probabilities                  #
############################################################

message('Downloading DepMap depletion probabilities...')

# Download and read in CRISPR gene dependency (DepMap 22Q2, post_chronos)
tmpfile <- tempfile()
download.file(url = 'https://ndownloader.figshare.com/files/34990033', destfile = tmpfile, quiet = T)
depmap_gene_dependency <- read.csv(tmpfile, check.names = F)
colnames(depmap_gene_dependency)[1] <- 'depMapID'

############################################################
# Calculate cell line essentiality totals                  #
############################################################

message('Calculating essentiality frequency per gene...')

# Calculate number of cell lines in which gene is essential
depmap_gene_dependency <- depmap_gene_dependency %>% 
  gather(gene_label, depmap_depletion_probability, -depMapID) %>%
  separate(gene_label, sep = " ", into = c('gene', 'entrez_id')) %>%
  mutate('entrez_id' = str_replace_all(entrez_id, "[()]", "")) %>% 
  mutate('is_depmap_dependent' = ifelse(depmap_depletion_probability > 0.5, 1, 0))

message('Calculating number of cell lines per gene...')

# For each gene, calculate total number of cell lines screened and number of cell lines in which it is depleted
depmap_gene_dependency.summary <- depmap_gene_dependency %>%
  drop_na() %>%
  group_by(gene) %>%
  summarise('n_cell_lines' = n(), 'n_cell_lines_dependent' = sum(is_depmap_dependent), .groups = 'keep') %>%
  unite('n_depmap_dependent_cell_lines', n_cell_lines_dependent, n_cell_lines, sep = '/')

message('Summarising essentiality statistics...')

# Combine frequency and totals
depmap_gene_dependency <- depmap_gene_dependency %>%
  left_join(depmap_gene_dependency.summary, by = 'gene')

############################################################
# Filter for screened lines and genes                      #
############################################################

message('Filter to include only screened genes and cell lines...')

# Add DepMap dependencies
depmap_gene_dependency.filt <- depmap_gene_dependency %>%
  filter(gene %in% c(expanded_library$targetA, expanded_library$targetB)) %>% # Filter for screened genes
  filter(depMapID %in% unique(sample_mapping$depMapID)) # Filter for screened cell lines
  
############################################################
# Expand by library metadata and cell line                 #
############################################################

message('Expanding library by cell line...')

depmap_gene_dependency.ann <- expand.grid(expanded_library$sorted_gene_pair, sample_mapping$depMapID) %>%
  unique() %>%
  select('sorted_gene_pair' = 'Var1', 'depMapID' = 'Var2' ) %>%
  filter(!is.na(depMapID))

message('Adding additional library metadata...')

depmap_gene_dependency.ann <- depmap_gene_dependency.ann %>%
  left_join(expanded_library %>% select(sorted_gene_pair, targetA, targetB), by = 'sorted_gene_pair', relationship = 'many-to-many') %>%
  unique()

message('Adding depletion probabilities for targetA...')

depmap_gene_dependency.ann <- depmap_gene_dependency.ann %>%
  left_join(depmap_gene_dependency.filt %>% 
              select(-entrez_id) %>%
              rename_at(vars(-gene, -depMapID), ~ paste('targetA', ., sep = '__')), by = c('targetA' = 'gene', 'depMapID')) %>%
  unique()

message('Adding depletion probabilities for targetB...')

depmap_gene_dependency.ann <- depmap_gene_dependency.ann %>%
  left_join(depmap_gene_dependency.filt %>% 
              select(-entrez_id) %>%
              rename_at(vars(-gene, -depMapID), ~ paste('targetB', ., sep = '__')), by = c('targetB' = 'gene', 'depMapID')) %>%
  unique()

message('Adding sample mapping...')

depmap_gene_dependency.ann <- depmap_gene_dependency.ann %>%
  left_join(sample_mapping %>% select(cell_line_label, cancer_type, depMapID), by = 'depMapID', relationship = 'many-to-many') %>%
  unique()
              
message('Postprocessong...')

depmap_gene_dependency.ann <- depmap_gene_dependency.ann %>%
  select(depMapID, cell_line_label, cancer_type, 
         sorted_gene_pair, 
         targetA, targetB, 
         targetA__is_depmap_dependent, targetA__n_depmap_dependent_cell_lines,
         targetB__is_depmap_dependent, targetB__n_depmap_dependent_cell_lines)

############################################################
# Write outputs to file                                    #
############################################################

message('Writing to file...')

# Save as a TSV
depmap_gene_dependency.filepath <- file.path(opt$dir, 'DATA', 'postprocessing', 'intermediate_tables', 'depmap_depletion_probability.tsv')
write.table(depmap_gene_dependency.ann, file = depmap_gene_dependency.filepath, sep = "\t", row.names = F, quote = F)
message(paste0('DepMap depletion probability TSV written to:', depmap_gene_dependency.filepath))

# Save as an Rdata object
depmap_gene_dependency.rds.filepath <- file.path(opt$dir, 'DATA', 'RDS', 'postprocessing', 'depmap_depletion_probability.rds')
saveRDS(depmap_gene_dependency.ann, file = depmap_gene_dependency.rds.filepath)
message(paste0('DepMap depletion probability RDS written to:', depmap_gene_dependency.rds.filepath))

message('Done.')