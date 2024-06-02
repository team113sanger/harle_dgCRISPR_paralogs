suppressPackageStartupMessages(suppressWarnings(library(optparse)))
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))
suppressPackageStartupMessages(suppressWarnings(library(scales)))
suppressPackageStartupMessages(suppressWarnings(library(ggpubr)))
suppressPackageStartupMessages(suppressWarnings(library(ggsci)))
suppressPackageStartupMessages(suppressWarnings(library(GGally)))
suppressPackageStartupMessages(suppressWarnings(library(ggrepel)))

############################################################
# OPTIONS                                                  #
############################################################

option_list = list(
  make_option(c("-d", "--dir"), type = "character",
              help = "full path to repository", metavar = "character"),
  make_option(c("-c", "--counts"), type = "character",
              help = "full path to count matrix", metavar = "character"),
  make_option(c("-l", "--lfc"), type = "character",
              help = "full path to lfc matrix", metavar = "character"),
  make_option(c("-m", "--mapping"), type = "character",
              help = "full path to sample mapping", metavar = "character"),
  make_option(c("--annotations"), type = "integer",
              help = "number of annotation columns", metavar = "integer"),
  make_option(c("--helper"), type = "character",
              help = "full path to helper functions", metavar = "character"),
  make_option(c("-r", "--rds"), type = "character",
              help = "full path to RDS directory", metavar = "character"),
  make_option(c("--out"), type = "character",
              help = "full path to output directory", metavar = "character"))

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
helper_path <- ifelse(is.null(opt$helper), file.path(repo_path, 'SCRIPTS', 'QC', 'helper.R'), opt$helper)
if (!file.exists(helper_path)){
  stop(paste('Helper file does not exist:', helper_path))
}
source(helper_path)

# Check output path exists
check_dir_exists(opt$out)

# Set RDS path
rds_path <- ifelse(is.null(opt$rds), file.path(repo_path, 'DATA', 'RDS', 'QC'), opt$rds)
check_dir_exists(rds_path)

############################################################
# Sample mapping                                           #
############################################################

check_file_exists(opt$mapping)

message(paste("Reading sample annotations from:", opt$mapping))

# Read in sample mapping
sample_mapping <- read_sample_metadata(opt$mapping)

# Cancer type palette
cancer_type_pal <- pal_npg(alpha = 0.3)(5)
names(cancer_type_pal) <- c(levels(sample_mapping$cancer_type), 'Control', 'Summary')

############################################################
# Count matrix                                             #
############################################################

check_file_exists(opt$counts)

message(paste("Reading count matrix from:", opt$counts))

# Read in count matrix
count_matrix <- read.delim(opt$counts, sep = "\t", header = T, check.names = F)

############################################################
# LFC matrix.                                              #
############################################################

check_file_exists(opt$lfc)

message(paste("Reading fold change matrix from:", opt$lfc))

# Read in lfc matrix
lfc_matrix <- read.delim(opt$lfc, sep = "\t", header = T, check.names = F)

############################################################
# Essential LFC dropout per cancer type                    #
############################################################

message("Gathering LFC matrix...")

# Gather log fold change matrix
lfc_matrix.narrow <- lfc_matrix %>%
  gather(sample_label, lfc, -all_of(colnames(lfc_matrix)[1:opt$annotations])) %>%
  left_join(sample_mapping, by = 'sample_label')

# Set guide group palette
sgrna_group_pal <- pal_d3(alpha = 0.3)(2)
names(sgrna_group_pal) <- c('Non-essential', 'Essential')

# Loop over cancer types
for (ct in levels(sample_mapping$cancer_type)) {
  message(paste("Plotting essential dropout:", ct))
  # Subset the data by cancer type
  d <- lfc_matrix.narrow %>%
    filter(cancer_type == ct & sgrna_group %in% c('Essential', 'Non-essential'))
  # Plot Essential and Non-Essential density by cell line and replicate
  ess_noness_density <- plot_ess_noness_lfc_density(d)
  # Save plot
  ess_noness_density_path <- file.path(opt$out, paste0('normalised_sample_lfc_', str_replace(ct, ' ', '_'), '.png'))
  ggsave(filename = ess_noness_density_path, plot = ess_noness_density, device = 'png', dpi = 300 , width = 3000, height = 4000, units = 'px')
  message(paste("Essential dropout plot written to:", ess_noness_density_path))
}

############################################################
# NNMD                                                     #
############################################################

message("Calculating NNMD...")

# Calculate NNMD (separation of Essential and Non-essential)
nnmd_results <- calculate_nnmd(lfc_matrix.narrow) %>%
  left_join(sample_mapping, by = 'sample_label')
nnmd_results$sample_label <- factor(nnmd_results$sample_label, levels = sample_mapping$sample_label)

# Write NNMD results to TSV
nnmd_results_path <- file.path(opt$out, 'normalised_LFC_NNMD.tsv')
write.table(nnmd_results, nnmd_results_path, row.names = F, quote = F, sep = "\t")
message(paste("NNMD results written to:", nnmd_results_path))

# Plot NMMD
message("Plotting NNMD...")
nnmd_plot_path <- file.path(opt$out, 'normalised_LFC_NNMD.png')
nnmd_plot <- plot_nnmd(nnmd_results)
ggsave(filename = nnmd_plot_path, plot = nnmd_plot, device = 'png', dpi = 300 , width = 4000, height = 2000, units = 'px')
message(paste("NNMD plot written to:", nnmd_plot_path))

# Get samples with NNMD > -2
nnmd_fail <- nnmd_results %>%
  filter(NNMD > -2) %>%
  pull(sample_label)

message(paste("Samples with NNMD > -2:", paste(nnmd_fail, collapse = ', ')))

# Write samples to exclude to file
nnmd_fail_path <- file.path(opt$out, 'samples_to_remove.txt')
write_lines(nnmd_fail, nnmd_fail_path)
message(paste("Samples to exclude written to:", nnmd_fail_path))

############################################################
# RDS outputs                                              #
############################################################

# Determine file extensions
tsv_ext <- 'tsv'
rds_ext <- 'rds'

# Write NNMD to RDS
nnmd.rds.path <- file.path(rds_path, paste0('normalised_lfc_nnmd', rds_ext))
saveRDS(nnmd_results, file = nnmd.rds.path)
message(paste("NNMD results RDS written to:", nnmd.rds.path))

message('Done.')

