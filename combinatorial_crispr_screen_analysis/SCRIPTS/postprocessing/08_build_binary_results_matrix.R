suppressPackageStartupMessages(suppressWarnings(library(optparse)))
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))

## NOTE: THIS SCRIPT REQUIRES STEPS 1-7 TO HAVE ALREADY BEEN RUN!

############################################################
# OPTIONS                                                  #
############################################################

option_list = list(
  make_option(c("-d", "--dir"), type = "character",
              help = "repository directory path", metavar = "character")
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
# Read in combined gene-level results                      #
############################################################

message('Reading combined results...')

# Combined gene level results
combined_results <- readRDS(file.path(repo_path, 'DATA', 'RDS', 'postprocessing', paste0('combined_gene_level_results.rds')))

# Replace spaces in cancer_type with underscores
combined_results <- combined_results %>%
  mutate(cancer_type = str_replace(cancer_type, ' ', '_'))

############################################################
# Summarise number of Bassik hits by cancer type           #
############################################################

message('Summarising Bassik hits...')

cancer_type_summary <- combined_results %>%
  group_by(sorted_gene_pair, cancer_type) %>%
  summarise(n_cell_lines = sum(is_bassik_hit), .groups = 'keep') %>%
  spread(cancer_type, n_cell_lines) %>%
  ungroup() %>%
  rowwise() %>% 
  mutate('total' = sum(c_across(`Lung_NSCLC`:Pancreas))) %>%
  rename_at(vars(-sorted_gene_pair), ~ paste('bassik', ., sep = '__'))

############################################################
# Summarise expression by cancer type                      #
############################################################

message('Summarising targetA expression...')

targetA_expression_summary <- combined_results %>%
  select(sorted_gene_pair, cancer_type, targetA__is_expressed) %>%
  drop_na() %>%
  group_by(sorted_gene_pair, cancer_type) %>%
  summarise(n_cell_lines = sum(targetA__is_expressed), .groups = 'keep') %>%
  spread(cancer_type, n_cell_lines) %>%
  ungroup() %>%
  rowwise() %>% 
  mutate('total' = sum(c_across(`Lung_NSCLC`:Pancreas))) %>%
  rename_at(vars(-sorted_gene_pair), ~ paste('targetA_expression', ., sep = '__'))

message('Summarising targetB expression...')

targetB_expression_summary <- combined_results %>%
  select(sorted_gene_pair, cancer_type, targetB__is_expressed) %>%
  drop_na() %>%
  group_by(sorted_gene_pair, cancer_type) %>%
  summarise(n_cell_lines = sum(targetB__is_expressed), .groups = 'keep') %>%
  spread(cancer_type, n_cell_lines) %>%
  ungroup() %>%
  rowwise() %>% 
  mutate('total' = sum(c_across(`Lung_NSCLC`:Pancreas))) %>%
  rename_at(vars(-sorted_gene_pair), ~ paste('targetB_expression', ., sep = '__'))

message('Averaging targetA expression...')

targetA_expression_avg <- combined_results %>%
  select(sorted_gene_pair, cancer_type, targetA__tpm_log2ps1) %>%
  drop_na() %>%
  group_by(sorted_gene_pair, cancer_type) %>%
  summarise(avg_expression = round(mean(targetA__tpm_log2ps1), 2), .groups = 'keep') %>%
  spread(cancer_type, avg_expression) %>%
  ungroup() %>%
  rename_at(vars(-sorted_gene_pair), ~ paste('targetA_avg_tpm_log2ps1', ., sep = '__'))

message('Averaging targetB expression...')

targetB_expression_avg <- combined_results %>%
  select(sorted_gene_pair, cancer_type, targetB__tpm_log2ps1) %>%
  drop_na() %>%
  group_by(sorted_gene_pair, cancer_type) %>%
  summarise(avg_expression = round(mean(targetB__tpm_log2ps1), 2), .groups = 'keep') %>%
  spread(cancer_type, avg_expression) %>%
  ungroup() %>%
  rename_at(vars(-sorted_gene_pair), ~ paste('targetB_avg_tpm_log2ps1', ., sep = '__'))

############################################################
# Summarise singles depletion by cancer type               #
############################################################

message('Summarising targetA singles depletion...')

targetA_depletion_summary <- combined_results %>%
  select(sorted_gene_pair, cancer_type, targetA__is_single_depleted) %>%
  drop_na() %>%
  group_by(sorted_gene_pair, cancer_type) %>%
  summarise(n_cell_lines = sum(targetA__is_single_depleted), .groups = 'keep') %>%
  spread(cancer_type, n_cell_lines) %>%
  ungroup() %>%
  rowwise() %>% 
  mutate('total' = sum(c_across(`Lung_NSCLC`:Pancreas))) %>%
  rename_at(vars(-sorted_gene_pair), ~ paste('targetA_single_depleted', ., sep = '__'))

message('Summarising targetB singles depletion...')

targetB_depletion_summary <- combined_results %>%
  select(sorted_gene_pair, cancer_type, targetB__is_single_depleted) %>%
  drop_na() %>%
  group_by(sorted_gene_pair, cancer_type) %>%
  summarise(n_cell_lines = sum(targetB__is_single_depleted), .groups = 'keep') %>%
  spread(cancer_type, n_cell_lines) %>%
  ungroup() %>%
  rowwise() %>% 
  mutate('total' = sum(c_across(`Lung_NSCLC`:Pancreas))) %>%
  rename_at(vars(-sorted_gene_pair), ~ paste('targetB_single_depleted', ., sep = '__'))

############################################################
# Spread Bassik binary results.                            #
############################################################

message('Building binary matrix...')

# Spread Bassik binary results by cell line label 
data <- combined_results %>% 
  select(sorted_gene_pair, cell_line_label, is_bassik_hit) %>%
  spread(cell_line_label, is_bassik_hit)

# Add Bassik cancer type summary
data <- data %>%
  left_join(cancer_type_summary, by = 'sorted_gene_pair')
  
# Add expression and depletion
data <- data %>%
  left_join(targetA_expression_summary, by = 'sorted_gene_pair') %>%
  left_join(targetA_expression_avg, by = 'sorted_gene_pair') %>%
  left_join(targetB_expression_summary, by = 'sorted_gene_pair') %>%
  left_join(targetB_expression_avg, by = 'sorted_gene_pair') %>%
  left_join(targetA_depletion_summary, by = 'sorted_gene_pair') %>%
  left_join(targetB_depletion_summary, by = 'sorted_gene_pair')

# Add common essentiality
data <- data %>%
  left_join(combined_results %>%
              select(sorted_gene_pair,
                     targetA__is_common_essential_bagel2:targetB__is_common_nonessential_depmap_ceres, 
                     targetA__n_depmap_dependent_cell_lines, targetB__n_depmap_dependent_cell_lines) %>%
              unique() %>%
              filter(!is.na(targetA__n_depmap_dependent_cell_lines) & !is.na(targetB__n_depmap_dependent_cell_lines)), by = 'sorted_gene_pair')

############################################################
# Write outputs to file                                    #
############################################################

message('Writing to file...')

# Save as a TSV
data.filepath <- file.path(opt$dir, 'DATA', 'postprocessing', paste0('combined_gene_level_results.binary.tsv'))
write.table(data, file = data.filepath, sep = "\t", row.names = F, quote = F)
message(paste0('Combined gene level binary results TSV written to:', data.filepath))

# Save as an Rdata object
data.rds.filepath <- file.path(opt$dir, 'DATA', 'RDS', 'postprocessing', paste0('combined_gene_level_results.binary.rds'))
saveRDS(data, file = data.rds.filepath)
message(paste0('Combined gene level binary results RDS written to:', data.rds.filepath))

message('Done.')