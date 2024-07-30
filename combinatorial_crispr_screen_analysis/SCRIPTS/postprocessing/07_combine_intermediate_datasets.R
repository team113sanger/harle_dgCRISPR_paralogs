suppressPackageStartupMessages(suppressWarnings(library(optparse)))
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))

## NOTE: THIS SCRIPT REQUIRES STEPS 1-6 TO HAVE ALREADY BEEN RUN!

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
# Read in required objects                                 #
############################################################

# Bassik results
screen_bassik <- readRDS(file.path(repo_path, 'DATA', 'RDS', 'dual_guide', 'Bassik', 'summary_by_gene.rds'))

# MAGeCK and BAGEL single essentiality for screened cell lines (binary)
screen_singles <- readRDS(file.path(repo_path, 'DATA', 'RDS', 'postprocessing', paste0('combined_singles_results.binary.rds')))

# Cell Model Passport Total Copy Number
cmp_copy_number <- readRDS(file.path(repo_path, 'DATA', 'RDS', 'postprocessing', 'cmp_total_copy_number.rds'))

# Cell Model Passport Cancer Driver Mutations
cmp_mutations <- readRDS(file.path(repo_path, 'DATA', 'RDS', 'postprocessing', 'cmp_cancer_driver_mutations.rds'))

# DepMap expression
depmap_expression <- readRDS(file.path(repo_path, 'DATA', 'RDS', 'postprocessing', 'depmap_expression.rds'))

# DepMap depletion probability
depmap_depletion <- readRDS(file.path(repo_path, 'DATA', 'RDS', 'postprocessing', 'depmap_depletion_probability.rds'))

# Common essentials and non-essentials
common_essentials <- readRDS(file.path(repo_path, 'DATA', 'RDS', 'postprocessing', 'depmap_and_bagel_common_genes.rds'))

############################################################
# Set up predefined cell lines and gene pairs              #
############################################################

# Limit gene pairs to those that made it through the dual Bassik analysis
cell_lines <- screen_bassik %>% 
  pull(cell_line_label) %>%
  unique()
gene_pairs <- unique(screen_bassik$sorted_gene_pair)

message('Creating gene pairs x cell lines matrix...')

data <- expand.grid(gene_pairs, cell_lines) %>%
  select('sorted_gene_pair' = 'Var1', 'qc_group' = 'Var2')

message('Adding additional library metadata...')

data <- data %>%
  left_join(expanded_library %>% 
              filter(guide_type == 'gene|gene') %>%
              select(sorted_gene_pair, targetA, targetB, sgrna_group) %>%
              unique(), by = 'sorted_gene_pair') %>%
  mutate(sgrna_group = str_replace(sgrna_group, "\\|.*", "")) 

message('Adding additional sample metadata...')
data <- data %>%
  left_join(sample_mapping %>% 
              select(qc_group, cell_line_label, cancer_type, depMapID) %>%
              unique(), by = 'qc_group') %>%
  relocate(qc_group, .after = 'cell_line_label')

############################################################
# Add screen Bassik results                                #
############################################################

message('Adding screened Bassik results...')

data <- data %>% 
  left_join(screen_bassik, by = c('sorted_gene_pair', 'targetA', 'targetB', 'cell_line_label', 'cancer_type'))

############################################################
# Add screen single essentiality results                   #
############################################################

message('Adding screened single essentiality results...')

data <- data %>% 
  left_join(screen_singles, by = c('sorted_gene_pair', 'targetA', 'targetB', 'cell_line_label', 'cancer_type'))

############################################################
# Add common essentials                                    #
############################################################

message('Adding BAGEL and DepMap common essentials...')

data <- data %>% 
  left_join(common_essentials, by = c('sorted_gene_pair', 'targetA', 'targetB'))

############################################################
# Add DepMap depletion probability                         #
############################################################

message('Adding DepMap depletion probabilities...')

data <- data %>% 
  left_join(depmap_depletion, by = c('sorted_gene_pair', 'targetA', 'targetB', 'cell_line_label', 'cancer_type', 'depMapID'))

############################################################
# Add DepMap expression                                    #
############################################################

message('Adding DepMap gene expression...')

data <- data %>% 
  left_join(depmap_expression, by = c('sorted_gene_pair', 'targetA', 'targetB', 'cell_line_label', 'cancer_type', 'depMapID'))

############################################################
# Add Cell Model Passport copy number                      #
############################################################

message('Adding CMP copy number...')

data <- data %>% 
  left_join(cmp_copy_number, by = c('sorted_gene_pair', 'targetA', 'targetB', 'cell_line_label', 'cancer_type', 'depMapID')) %>%
  relocate('CMP_model_id', .after = depMapID)

############################################################
# Add binary for Bassik                                    #
############################################################

message('Parsing Bassik results...')

# Set Bassik hit as residual < -0.5 and fdr < 0.05
data <- data %>%
  mutate('is_bassik_hit' = ifelse(mean_norm_gi < -0.5 & fdr < 0.01 & targetA__is_single_depleted == 0 & targetB__is_single_depleted == 0, 1, 0), .after = 'fdr')

# Add % of significant guides
data <- data %>% 
  rowwise() %>% 
  mutate('pct_sig_guide_pairs' = round((n_sig_guide_pairs / n_guide_pairs) * 100, 2), .after = 'n_sig_guide_pairs')

############################################################
# Add binary for expression                                #
############################################################

message('Parsing DepMap expression...')

# DepMap TPM = log2(TPM+1)
# Threshold of TPM = 1 or log2(1 + 1) = 1
data <- data %>% 
  rowwise() %>% 
  mutate('targetA__is_expressed' = ifelse(targetA__tpm_log2ps1 >= 1, 1, 0), .after = 'is_bassik_hit') %>% 
  mutate('targetB__is_expressed' = ifelse(targetB__tpm_log2ps1 >= 1, 1, 0), .after = 'targetA__is_expressed') 

############################################################
# Add binary for single essentiality                       #
############################################################

message('Parsing singles essentiality results...')

# For gene to be essential, it must be identified by both BAGEL and MAGeCK
data <- data %>% 
  rowwise() %>% 
  mutate('targetA__is_single_depleted' = ifelse(targetA__is_depleted_bagel == 1 & targetA__is_depleted_mageck == 1, 1, 0), .after = 'is_bassik_hit') %>% 
  mutate('targetB__is_single_depleted' = ifelse(targetB__is_depleted_bagel == 1 & targetB__is_depleted_mageck == 1, 1, 0), .after = 'targetA__is_single_depleted') 

############################################################
# Add Cell Model Passport cancer driver mutations          #
############################################################

message('Adding CMP cancer driver mutations...')

cmp_mutations <- cmp_mutations %>%
  select(-cancer_driver) %>%
  group_by(depMapID, CMP_model_id, cell_line_label, cancer_type) %>% 
  summarise_at(vars(-group_cols()), ~ paste(na.omit(.), collapse = "|")) %>%
  rename_at(vars(-cell_line_label, -cancer_type, -depMapID, -CMP_model_id), ~ paste('CMP_mutation', ., sep = '__'))

data <- data %>% 
  left_join(cmp_mutations, by = c('cell_line_label', 'cancer_type', 'depMapID', 'CMP_model_id'))

############################################################
# Write outputs to file                                    #
############################################################

message('Writing to file...')

# Save as a TSV
data.filepath <- file.path(opt$dir, 'DATA', 'postprocessing', paste0('combined_gene_level_results.tsv'))
write.table(data, file = data.filepath, sep = "\t", row.names = F, quote = F)
message(paste0('Combined gene level results TSV written to:', data.filepath))

# Save as an Rdata object
data.rds.filepath <- file.path(opt$dir, 'DATA', 'RDS', 'postprocessing', paste0('combined_gene_level_results.rds'))
saveRDS(data, file = data.rds.filepath)
message(paste0('Combined gene level results RDS written to:', data.rds.filepath))

message('Done.')
