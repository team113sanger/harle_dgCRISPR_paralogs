suppressPackageStartupMessages(suppressWarnings(library(optparse)))
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))
suppressPackageStartupMessages(suppressWarnings(library(ggsci)))
suppressPackageStartupMessages(suppressWarnings(library(ggpubr)))
suppressPackageStartupMessages(suppressWarnings(library(pROC)))
suppressPackageStartupMessages(suppressWarnings(library(scales)))

############################################################
# OPTIONS                                                  #
############################################################

option_list = list(
  make_option(c("-d", "--dir"), type = "character",
              help = "full path to repository", metavar = "character"),
  make_option(c("-m", "--mapping"), type = "character",
              help = "full path to sample mapping", metavar = "character"),
  make_option(c("-f", "--fc"), type = "character",
              help = "full path to fold change matrix", metavar = "character"),
  make_option(c("--annotations"), type = "character",
              help = "number of annotation columns", metavar = "character"),
  make_option(c("--mageck"), type = "character",
              help = "full path to MAGeCK results", metavar = "character"),
  make_option(c("--bagel"), type = "character",
              help = "full path to BAGEL2 results", metavar = "character"),
  make_option(c("--ess"), type = "character",
              help = "full path to essential gene file", metavar = "character"),
  make_option(c("--noness"), type = "character",
              help = "full path to non-essesntial gene file", metavar = "character"),
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
helper_path <- ifelse(is.null(opt$helper), file.path(repo_path, 'SCRIPTS', 'single_guide', 'helper.R'), opt$helper)
if (!file.exists(helper_path)){
  stop(paste('Helper file does not exist:', helper_path))
}

source(helper_path)

# Check output path exists
check_dir_exists(opt$out)

# Set RDS path
rds_path <- ifelse(is.null(opt$rds), file.path(repo_path, 'DATA', 'RDS', 'single_guide'), opt$rds)
check_dir_exists(rds_path)

############################################################
# Sample mapping                                           #
############################################################

check_file_exists(opt$mapping)

message(paste("Reading sample annotations from:", opt$mapping))

# Read in sample mapping
sample_mapping <- read_sample_metadata(opt$mapping)

# Get stripped cell lines and cell line labels
cell_lines <- list.dirs(opt$mageck, full.names = F, recursive = F)
cell_line_mapping <- sample_mapping %>%
  filter(stripped_cell_line_name %in% cell_lines) %>%
  select(stripped_cell_line_name, cell_line_label, cancer_type) %>%
  unique()

# Cancer type palette
cancer_type_pal <- pal_npg(alpha = 0.3)(5)
names(cancer_type_pal) <- c(levels(sample_mapping$cancer_type), 'Control', 'Summary')

############################################################
# Scaled LFC matrix                                        #
############################################################

# Set LFC matrix path
lfc_matrix_file <- opt$fc

message(paste("Reading LFC matrix:", lfc_matrix_file))

# Check LFC matrix file exists
if (!file.exists(lfc_matrix_file)) {
  stop(sprintf("LFC matrix file does not exist: %s", lfc_matrix_file))
}

# Read in LFC matrix
lfc_matrix <- read.delim(lfc_matrix_file, sep = "\t", header = T, check.names = F)

# Set annotation column names
annotation_colnames <- colnames(lfc_matrix)[1:opt$annotations]

############################################################
# MAGeCK gene-level results                                #
############################################################

mageck_gene_results <- collate_mageck_gene_results(opt$mageck, sample_mapping)

############################################################
# BAGEL2 gene-level results                                #
############################################################

# Note: these are collated from the guide-level results
# Average BF across replicates
# Average BF across guides
bagel_gene_results <- collate_bagel_gene_results(opt$bagel, sample_mapping)

# Read in essential genes
message('Reading in essential and non-essential genes...')
ess <- read.delim(opt$ess, header = T, check.names = F, sep = "\t")
noness <- read.delim(opt$noness, header = T, check.names = F, sep = "\t")

# Add classifications to BAGEL2 results
message('Classifying BAGEL results...')
bagel_gene_results <- bagel_gene_results %>%
  mutate('classification' = case_when(GENE %in% ess$GENE ~ 'Essential',
                                      GENE %in% noness$GENE ~ 'Non-essential',
                                      TRUE ~ 'unknown'))

# Get ROC metrics
message('Getting ROC metrics...')
roc_list <- roc_metrics(bagel_gene_results)

# Annotate BAGEL pass/fail
message('Determining BAGEL2 pass/fail...')
bagel_gene_results.ann <- bagel_gene_results %>%
  left_join(roc_list[['thresholds']] %>% select(dataset, cell_line_label, 'threshold' = thresholds), by = c('dataset', 'cell_line_label')) %>%
  left_join(sample_mapping %>% select(cell_line_label, cancer_type) %>% unique(), by = 'cell_line_label') %>%
  relocate(cancer_type, .after = 'cell_line_label') %>%
  relocate(classification, .after = 'GENE') %>%
  mutate('scaled_BF' = BF - threshold) %>%
  mutate('is_depleted_bagel' = ifelse(scaled_BF > 0, 1, 0))

# Prepare data for ROC plot
message('Preparing BAGEL2 ROC...')
roc_coords <- roc_list[['roc_coords']] %>%
  left_join(sample_mapping %>% select(cell_line_label, cancer_type) %>% unique(), by = 'cell_line_label') %>%
  relocate(cancer_type, .after = 'cell_line_label')

# Plot ROC
message('Plotting BAGEL2 ROC...')
bagel_roc_plot <- plot_roc(roc_coords, cancer_type_pal)
bagel_roc_plot_path <- file.path(opt$out, 'scaled_bagel2_roc.png')
ggsave(filename = bagel_roc_plot_path, plot = bagel_roc_plot, device = 'png', dpi = 300 , width = 3000, height = 3000, units = 'px')
message(paste('BAGEL2 ROC plot saved to:', bagel_roc_plot_path))

############################################################
# Collate binary matrix.                                   #
############################################################

# Take only combined data set forward
message('Combining MAGeCK and BAGEL datasets (combined)...')
mageck_gene_results.combined <- mageck_gene_results %>% filter(dataset == 'combined')
bagel_gene_results.combined <- bagel_gene_results.ann %>% filter(dataset == 'combined')
combined_data <- mageck_gene_results.combined %>%
  left_join(bagel_gene_results.combined, by = c('dataset', 'cell_line_label', 'gene' = 'GENE')) %>%
  select(-dataset) %>%
  relocate(cancer_type, .after = 'cell_line_label') %>%
  relocate(classification, .after = 'gene') %>%
  relocate(is_enriched_mageck, is_depleted_mageck, .before = 'is_depleted_bagel')

message('Preparing binary matrix (combined)...')
gp <- lfc_matrix %>% filter(guide_type == 'gene|gene') %>% pull(sorted_gene_pair) %>% unique()
cl <- combined_data %>% pull(cell_line_label) %>% unique()
binary_matrix <- expand.grid(gp, cl) %>%
  select('sorted_gene_pair' = 'Var1', 'cell_line_label' = 'Var2') %>%
  left_join(lfc_matrix %>% select(sorted_gene_pair, targetA, targetB, sgrna_group) %>% unique(), by = 'sorted_gene_pair') %>%
  left_join(sample_mapping %>% select(cell_line_label, cancer_type) %>% unique(), by = 'cell_line_label') %>%
  relocate(cell_line_label, cancer_type, .before = 'sorted_gene_pair')

message('Adding single essentiality for targetA (combined)...')
binary_matrix.targetA <- binary_matrix %>%
  left_join(combined_data %>% select(cell_line_label, gene, is_enriched_mageck, is_depleted_mageck, is_depleted_bagel), by = c('cell_line_label', 'targetA' = 'gene')) %>%
  rename_at(vars(-sgrna_group, -cancer_type, -sorted_gene_pair, -cell_line_label, -targetA, -targetB), ~ paste('targetA', ., sep = '__'))


message('Adding single essentiality for targetB (combined)...')
binary_matrix.targetB <- binary_matrix %>%
  left_join(combined_data %>% select(cell_line_label, gene, is_enriched_mageck, is_depleted_mageck, is_depleted_bagel), by = c('cell_line_label', 'targetB' = 'gene')) %>%
  rename_at(vars(-sgrna_group, -cancer_type, -sorted_gene_pair, -cell_line_label, -targetA, -targetB), ~ paste('targetB', ., sep = '__'))

message('Combining targetA and targetB binary matrices (combined)...')
binary_matrix <- binary_matrix %>%
  left_join(binary_matrix.targetA, by = c('cell_line_label', 'cancer_type', 'sorted_gene_pair', 'targetA', 'targetB', 'sgrna_group')) %>%
  left_join(binary_matrix.targetB, by = c('cell_line_label', 'cancer_type', 'sorted_gene_pair', 'targetA', 'targetB', 'sgrna_group')) %>%
  mutate('targetA__is_single_depleted' = ifelse(targetA__is_depleted_bagel == 1 & targetA__is_depleted_mageck == 1, 1, 0)) %>%
  mutate('targetB__is_single_depleted' = ifelse(targetB__is_depleted_bagel == 1 & targetB__is_depleted_mageck == 1, 1, 0))

message('Plotting number of target A/B genes depleted by cell line...')
n_genes_binary_depleted_plot <- plot_n_genes_binary_depleted(binary_matrix, cancer_type_pal)
n_genes_binary_depleted_plot_path <- file.path(opt$out, 'scaled_bagel2_n_genes_binary_depleted.png')
ggsave(filename = n_genes_binary_depleted_plot_path, plot = n_genes_binary_depleted_plot, device = 'png', dpi = 300 , width = 4000, height = 3000, units = 'px')
message(paste('Binary depletion n genes plot saved to:', n_genes_binary_depleted_plot_path))

############################################################
# Write tables                                             #
############################################################

# MAGeCK gene results
message("Writing MAGeCK gene results...")
mageck_gene_results_path <- file.path(opt$out, 'MAGeCK_gene_results.tsv')
write.table(mageck_gene_results, mageck_gene_results_path, row.names = F, sep = "\t", quote = F)
message(paste('MAGeCK gene results written to:', mageck_gene_results_path))

# BAGEL gene results
message("Writing BAGEL2 gene results...")
bagel_gene_results_path <- file.path(opt$out, 'BAGEL2_gene_results.tsv')
write.table(bagel_gene_results, bagel_gene_results_path, row.names = F, sep = "\t", quote = F)
message(paste('BAGEL gene results written to:', bagel_gene_results_path))

# ROC coordinates
message("Writing BAGEL2 ROC coordinates...")
bagel_roc_coords_path <- file.path(opt$out, 'BAGEL2_ROC_coordinates.tsv')
write.table(roc_coords, bagel_roc_coords_path, row.names = F, sep = "\t", quote = F)
message(paste('BAGEL ROC coordinates written to:', bagel_roc_coords_path))

# ROC threshold
message("Writing BAGEL2 ROC thresholds")
bagel_roc_thresholds_path <- file.path(opt$out, 'BAGEL2_ROC_thresholds.tsv')
write.table(roc_list[['thresholds']], bagel_roc_thresholds_path, row.names = F, sep = "\t", quote = F)
message(paste('BAGEL ROC thresholds written to:', bagel_roc_thresholds_path))

# Binary matrix
message("Writing binary results...")
binary_matrix_path <- file.path(opt$out, 'Singles_analysis_binary_matrix.tsv')
write.table(binary_matrix, binary_matrix_path, row.names = F, sep = "\t", quote = F)
message(paste('Binary results written to:', binary_matrix_path))

# Save binary results as RDS
binary_results.rds.path <- file.path(repo_path, 'DATA', 'RDS', 'postprocessing', paste0('combined_singles_results.binary.rds'))
saveRDS(binary_matrix, file = binary_results.rds.path)
message(paste("Combined MAGeCK and BAGEL binary results RDSwritten to:", binary_results.rds.path))