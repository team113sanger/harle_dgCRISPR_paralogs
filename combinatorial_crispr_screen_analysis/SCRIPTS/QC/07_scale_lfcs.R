suppressPackageStartupMessages(suppressWarnings(library(optparse)))
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))
suppressPackageStartupMessages(suppressWarnings(library(ggsci)))
suppressPackageStartupMessages(suppressWarnings(library(ggpubr)))
suppressPackageStartupMessages(suppressWarnings(library(scales)))
suppressPackageStartupMessages(suppressWarnings(library(gridExtra)))

############################################################
# OPTIONS                                                  #
############################################################

option_list = list(
  make_option(c("-d", "--dir"), type = "character",
              help = "full path to repository", metavar = "character"),
  make_option(c("-f", "--fc"), type = "character",
              help = "full path to lfc matrix", metavar = "character"),
  make_option(c("-m", "--mapping"), type = "character",
              help = "full path to sample mapping", metavar = "character"),
  make_option(c("-s", "--suffix"), type = "character",
              help = "file name suffix", metavar = "character"),
  make_option(c("--annotations"), type = "integer",
              help = "number of annotation columns", metavar = "integer"),
  make_option(c("--helper"), type = "character",
              help = "full path to helper functions", metavar = "character"),
  make_option(c("-r", "--rds"), type = "character",
              help = "full path to RDS directory", metavar = "character"),
  make_option(c("--plots_out"), type = "character",
              help = "full path to output directory", metavar = "character"),
  make_option(c("--lfc_out"), type = "character",
              help = "full path to output directory", metavar = "character"));

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
check_dir_exists(opt$lfc_out)
check_dir_exists(opt$plots_out)

# Set RDS path
rds_path <- ifelse(is.null(opt$rds), file.path(repo_path, 'DATA', 'RDS', 'QC'), opt$rds)
check_dir_exists(rds_path)

# Determine file extensions with suffix
if (is.null(opt$suffix)) {
  tsv_ext <- 'tsv'
  rds_ext <- 'rds'
} else {
  tsv_ext <- paste0(opt$suffix, '.tsv')
  rds_ext <- paste0(opt$suffix, '.rds')
}


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
# LFC matrix                                               #
############################################################

check_file_exists(opt$fc)

message(paste("Reading LFC matrix:", opt$fc))

# Read in LFC matrix
lfc_matrix <- read.delim(opt$fc, sep = "\t", header = T, check.names = F)

# Set annotation column names
annotation_colnames <- colnames(lfc_matrix)[1:opt$annotations]

############################################################
# Prepare LFC matrix                                       #
############################################################

message('Preparing LFC matrix...')

# Gather LFC matrix
lfc.narrow <- lfc_matrix %>%
  gather(sample_label, LFC, -all_of(annotation_colnames)) 

# Annotate with sample information
lfc.narrow <- lfc.narrow %>%
  left_join(sample_mapping %>% select(sample_label, cell_line_label, replicate, cancer_type), by = 'sample_label') %>%
  relocate('LFC', .after = 'cancer_type')

############################################################
# Scale LFC matrix                                         #
############################################################

# Get median of all replicates in a cell line by library type (e.g. safe|safe)
message('Calculating replicate medians per cell line by library type...')
lfc.narrow.medians <- lfc.narrow %>%
  group_by(cell_line_label, sgrna_group) %>%
  mutate(median_lfc_cell_line_lib_type = median(LFC)) %>%
  ungroup()

# Calculate safe and essential medians
message('Calculating safe and essential medians...')
sgrna_group_medians <- lfc.narrow.medians %>%
  filter(sgrna_group %in% c('Essential', 'Safe-targeting control')) %>%
  group_by(cell_line_label, sgrna_group) %>%
  summarise(median = median(LFC), .groups = 'keep') %>%
  spread(sgrna_group, median) %>%
  rename( safe_safe_median = `Safe-targeting control`, essential_median = `Essential`)

# Add a new column with the median of the safe_safes for that cell line
message('Combining median data frames...')
lfc.narrow.medians <- lfc.narrow.medians %>%
  left_join(sgrna_group_medians, by = 'cell_line_label')

# Calculated scaled LFC
# There will be a new essential median after minusing safe_safe median, so need to account for that by adding safe_safe_median to denominator
message('Scaling LFCs...')
lfc.narrow.scaled <- lfc.narrow.medians %>%
  mutate(scaled_LFC = ((LFC - safe_safe_median) / (safe_safe_median - essential_median)))

# Spread scaled LFCs
lfc.scaled.wide <- lfc.narrow.scaled %>%
  select(all_of(annotation_colnames), sample_label, scaled_LFC) %>%
  spread(sample_label, scaled_LFC)

# Check that the new matrix is the correct length
if (nrow(lfc_matrix) != nrow(lfc.scaled.wide)) {
  message(paste('Number of rows pre-scaling:', nrow(lfc_matrix)))
  message(paste('Number of rows post-scaling:', nrow(lfc.scaled.wide)))
  stop("Incorrect number of rows after scaling.")
}

# Write scaled LFCs to TSV
scaled_lfc.path <- file.path(opt$lfc_out, paste0('lfc_matrix.scaled.', tsv_ext))
write.table(lfc.scaled.wide, scaled_lfc.path, sep = "\t", quote = F, row.names = F)
message(paste("Scaled LFC TSV written to:", scaled_lfc.path))

############################################################
# Plot scaled LFC distribution                             #
############################################################

message("Building LFC distribution plot...")

# Combine scaled and unscaled LFC
lfc.narrow.all <- lfc.narrow.scaled %>% 
  select(sample_label, cell_line_label, cancer_type, sgrna_group, 'LFC' = scaled_LFC) %>%
  mutate('type' = 'scaled')
lfc.narrow.all <- rbind(lfc.narrow.all, lfc.narrow %>% select(sample_label, cell_line_label, cancer_type, sgrna_group, LFC) %>% mutate('type' = 'unscaled'))

# Set essentiality palette
sgrna_group_pal <- pal_d3(alpha = 0.3)(2)
    names(sgrna_group_pal) <- c('Non-essential', 'Essential')

# Plot LFC distributions
lfc_dist_path <- file.path(opt$plots_out, 'scaled_lfc_distributions.samples_removed.png')
scaled_lfc_boxplot <- plot_scaled_lfc_boxplot(lfc.narrow.all, sgrna_group_pal)
ggsave(filename = lfc_dist_path, plot = scaled_lfc_boxplot, device = 'png', dpi = 300 , width = 4000, height = 3000, units = 'px')

message(paste("LFC distribution plot saved to:", lfc_dist_path))

############################################################
# Plot scaled LFC distribution by guide type (violin)      #
############################################################

message("Building LFC distribution by guide source (violin)...")

# Plot LFC distributions by guide source
violin_guide_dist_path <- file.path(opt$plots_out, 'scaled_lfc_guide_source_violin.samples_removed.png')
scaled_lfc_source_violin <- plot_scaled_lfc_violin_by_guide_source(lfc.narrow.scaled)
ggsave(filename = violin_guide_dist_path, plot = scaled_lfc_source_violin, device = 'png', dpi = 300 , width = 3000, height = 4000, units = 'px')

message(paste("LFC guide source violin plot saved to:", violin_guide_dist_path))

############################################################
# Plot scaled LFC distribution by guide type (barplot)     #
############################################################

message("Building LFC distribution by guide source (barplot)...")

barplot_guide_dist_path <- file.path(opt$plots_out, 'scaled_lfc_guide_source_barplot.samples_removed.png')
scaled_lfc_barplot_dist <- plot_scaled_lfc_with_guide_type_distribution(lfc.scaled.wide, annotation_colnames)
ggsave(filename = barplot_guide_dist_path, plot = scaled_lfc_barplot_dist, device = 'png', dpi = 300 , width = 4000, height = 3000, units = 'px')

message(paste("LFC guide source violin plot saved to:", barplot_guide_dist_path))

############################################################
# RDS outputs                                              #
############################################################

# Write LFC medians to RDS
lfc_medians.rds.path <- file.path(rds_path, paste0('lfc_medians.unscaled.', rds_ext))
saveRDS(lfc.narrow.medians, file = lfc_medians.rds.path)
message(paste("Unscaled LFC median RDS written to:", lfc_medians.rds.path))

# Write scaled LFCs to RDS
scaled_lfc.rds.path <- file.path(rds_path, paste0('lfc_matrix.scaled.', rds_ext))
saveRDS(lfc.narrow.scaled, file = scaled_lfc.rds.path)
message(paste("Scaled LFC RDS written to:", scaled_lfc.rds.path))

message('Done.')