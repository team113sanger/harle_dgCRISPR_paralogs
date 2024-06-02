suppressPackageStartupMessages(suppressWarnings(library(optparse)))
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))

############################################################
# OPTIONS                                                  #
############################################################

option_list = list(
  make_option(c("-d", "--dir"), type = "character",
              help = "full path to repository", metavar = "character"),
  make_option(c("-m", "--mapping"), type = "character",
              help = "full path to sample mapping", metavar = "character"),
  make_option(c("-f", "--fc"), type = "character",
              help = "full path to lfc matrix", metavar = "character"),
  make_option(c("-c", "--counts"), type = "character",
              help = "full path to count matrix", metavar = "character"),
  make_option(c("--annotations"), type = "integer",
              help = "number of annotation columns", metavar = "integer"),
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

############################################################
# LFC matrix                                               #
############################################################

check_file_exists(opt$fc)

message(paste("Reading LFC matrix:", opt$fc))

# Read in LFC matrix
lfc_matrix <- read.delim(opt$fc, sep = "\t", header = T, check.names = F)

# message total number of guides
fnum_guides <- nrow(lfc_matrix)
message(paste("Guides in LFC matrix:", fnum_guides))

# message total number of gene pairs
fnum_gene_pairs <- lfc_matrix %>% pull(sorted_gene_pair) %>% unique() %>% length()
message(paste("Gene pairs in LFC matrix:", fnum_gene_pairs))

# Set annotation column names
lfc_annotation_colnames <- colnames(lfc_matrix)[1:opt$annotations]

############################################################
# Count matrix                                             #
############################################################

check_file_exists(opt$counts)

message(paste("Reading count matrix:", opt$counts))

# Read in count matrix
count_matrix <- read.delim(opt$counts, sep = "\t", header = T, check.names = F)

# message total number of guides
cnum_guides <- nrow(count_matrix)
message(paste("Guides in count matrix:", cnum_guides))

# message total number of gene pairs
cnum_gene_pairs <- count_matrix %>% pull(sorted_gene_pair) %>% unique() %>% length()
message(paste("Gene pairs in count matrix:", cnum_gene_pairs))

############################################################
# Annotation columns                                       #
############################################################

# Set annotation column names
count_annotation_colnames <- colnames(count_matrix)[1:opt$annotations]

# Check same number of annotations in both matrices
if(length(setdiff(count_annotation_colnames, lfc_annotation_colnames)) > 0) {
  stop("Different number annotation columns in count and LFC matrices.")
}

############################################################
# Samples per cell line                                    #
############################################################

message("Collating cell line and sample information...")

# Get list of samples in count matrix
sample_mapping.filt <- sample_mapping %>%
  filter(sample_label %in% colnames(count_matrix))

# Get replicates per cell line (stripped cell line name)
cell_line_samples <- list()
for (cl in unique(sample_mapping.filt$stripped_cell_line_name)) {
  cell_line_samples[[cl]] <- sample_mapping.filt %>%
    filter(stripped_cell_line_name == cl) %>%
    pull(sample_label) %>%
    unique()
}

############################################################
# Get single guide ids                                     #
############################################################

message("Extracting single guide ids...")

# Set list for guide ids
single_guide_ids <- list()
datasets <- c('combined', 'A', 'B')

# Get guide ids for combined (geneA and geneB) single gene-targeting guides and controls
single_guide_ids[['combined']] <- count_matrix %>%
  filter(guide_type %in% c("safe_targeting|safe_targeting", "gene|safe_targeting", "safe_targeting|gene")) %>%
  pull(id)
message(paste('Number of guides (AB):', length(single_guide_ids[['combined']])))

# Get guide ids for single geneA-targeting guides and controls
single_guide_ids[['A']] <- count_matrix %>%
  filter(guide_type %in% c("safe_targeting|safe_targeting", "gene|safe_targeting")) %>%
  pull(id)
message(paste('Number of guides (A):', length(single_guide_ids[['A']])))

# Get guide ids for single geneB-targeting guides and controls
single_guide_ids[['B']] <- count_matrix %>%
  filter(guide_type %in% c("safe_targeting|safe_targeting", "safe_targeting|gene")) %>%
  pull(id)
message(paste('Number of guides (B):', length(single_guide_ids[['B']])))

############################################################
# Convert count matrix for MAGeCK RRA                      #
############################################################

message("Converting singles count matrix for use with MAGeCK...")

# Create a list for MAGeCK-formatted count matrices
single_guides_mageck <- list()

for (d in datasets) {
  message(paste("Formatting MAGeCK count matrix for:", d))

  # Format library
  single_guides_mageck[[d]] <- count_matrix %>%
    filter(id %in% single_guide_ids[[d]]) %>%
    mutate('singles_target_gene' = ifelse(guide_type == "safe_targeting|safe_targeting", 'safe_targeting', singles_target_gene)) %>%
    rename(sgRNA = id) %>%
    rename(gene = singles_target_gene) %>%
    select(-any_of(count_annotation_colnames)) %>%
    unique()

  # Summarise library
  message(paste('Number of guides', d, ':', nrow(single_guides_mageck[[d]])))
  message(paste('Number of safe-targeting control guides', d, ':', nrow(single_guides_mageck[[d]] %>% filter(gene == 'safe_targeting'))))

  # Write library
  message(paste("Writing MAGeCK count matrix for:", d))
  singles_mageck.path <- file.path(opt$out, paste0('MAGeCK.singles_library.', d, '.tsv'))
  write.table(single_guides_mageck[[d]], singles_mageck.path, sep = "\t", quote = F, row.names = F)
  message(paste("Singles library (", d,  ") written to:", singles_mageck.path))
}

############################################################
# MAGeCK LSF commands                                      #
############################################################

# Prepare bsub commands
for (d in datasets) {
  message(paste("Preparing MAGeCK LSF commands for:", d))
  mageck_lsf_cmds <- list()
  for (cl in names(cell_line_samples)) {
    # Set output directory
    result_dir <- file.path(opt$dir, 'DATA', 'single_guide', '01_MAGeCK', cl, d)
    if (!dir.exists(result_dir)) { dir.create(result_dir, recursive = T) }

    # Set LFS parameters
    lfs_jobname <- paste('mageck', cl, d, sep = '_')
    lfs_out <- file.path(repo_path, 'LOGS', 'MAGeCK', paste0(lfs_jobname, '.o'))
    lfs_err <- file.path(repo_path, 'LOGS', 'MAGeCK', paste0(lfs_jobname, '.e'))
    lsf_mem <- " -R \"select[mem>4000] rusage[mem=4000] span[hosts=1]\" -M 4000"
    lfs_queue <- 'normal'

    # Set MAGeCK count matrix
    mageck_counts <- file.path(opt$out, paste0('MAGeCK.singles_library.', d, '.tsv'))

    # Set MAGeCK command
    mageck_cmd <- paste0("mageck test --norm-method 'none' --remove-zero 'none'",
                         ' -k ', mageck_counts,
                         " -t '", paste(cell_line_samples[[cl]], collapse = ','), "'",
                         " -c 'control_mean'",
                         " -n '", paste0(result_dir, '/MAGeCK'), "'")

    # Set LSF bsub command
    mageck_lsf_cmds[[cl]] <- paste('bsub -J', lfs_jobname,
                                    lsf_mem,
                                    '-o', lfs_out,
                                    '-e', lfs_err,
                                    '-q', lfs_queue,
                                    mageck_cmd)
  }

  # Write commands to file
  message(paste("Writing MAGeCK LSF commands for:", d))
  mageck_cmds_path <- file.path(opt$dir, 'SCRIPTS', 'single_guide', paste0('MAGeCK_bsub_commands_', d, '.sh'))
  write_lines(mageck_lsf_cmds, file(mageck_cmds_path), append = F)
  message(paste("MAGeCK LSF commands written to:", mageck_cmds_path))
}

############################################################
# Convert LFC matrix for BAGEL2                            #
############################################################

message("Converting singles lfc matrix for use with BAGEL2...")

# Create a list for BAGEL2-formatted LFC matrices
single_guides_bagel <- list()

for (d in datasets) {
  message(paste("Formatting BAGEL2 LFC matrix for:", d))

  # Format library
  single_guides_bagel[[d]] <- lfc_matrix %>%
    filter(id %in% single_guide_ids[[d]]) %>%
    mutate('singles_target_gene' = ifelse(guide_type == "safe_targeting|safe_targeting", 'safe_targeting', singles_target_gene)) %>%
    rename(GENE_CLONE = id) %>%
    rename(GENE = singles_target_gene) %>%
    select(-any_of(lfc_annotation_colnames)) %>%
    unique()

  # Summarise library
  message(paste('Number of guides', d, ':', nrow(single_guides_bagel[[d]])))
  message(paste('Number of safe-targeting control guides', d, ':', nrow(single_guides_bagel[[d]] %>% filter(GENE == 'safe_targeting'))))

  # Write library
  message(paste("Writing BAGEL2 count matrix for:", d))
  singles_bagel.path <- file.path(opt$out, paste0('BAGEL2.singles_library.', d, '.tsv'))
  write.table(single_guides_bagel[[d]], singles_bagel.path, sep = "\t", quote = F, row.names = F)
  message(paste("Singles library (", d,  ") written to:", singles_bagel.path))
}

############################################################
# BAGEL2 LSF commands                                      #
############################################################

# Prepare bsub commands
for (d in datasets) {
  message(paste("Preparing BAGEL2 LSF commands for:", d))
  bagel_lsf_cmds <- list()
  for (cl in names(cell_line_samples)) {
    # Set output directory
    result_dir <- file.path(opt$dir, 'DATA', 'single_guide', '02_BAGEL2', cl, d)
    if (!dir.exists(result_dir)) { dir.create(result_dir, recursive = T) }

    # Set LFS parameters
    lfs_jobname <- paste('bagel', cl, d, sep = '_')
    lfs_out <- file.path(repo_path, 'LOGS', 'BAGEL2', paste0(lfs_jobname, '.o'))
    lfs_err <- file.path(repo_path, 'LOGS', 'BAGEL2', paste0(lfs_jobname, '.e'))
    lsf_mem <- " -R \"select[mem>4000] rusage[mem=4000] span[hosts=1]\" -M 4000"
    lfs_queue <- 'normal'

    # Set BAGEL2 LFC matrix
    bagel_lfc <- file.path(opt$out, paste0('BAGEL2.singles_library.', d, '.tsv'))

    # Set BAGEL2 gene command
    bagel2_gene_out <- file.path(result_dir, 'BAGEL2.gene.bf')
    bagel2_gene_cmd <- paste0("BAGEL.py bf",
                               ' -i ', bagel_lfc,
                               " -c '", paste(cell_line_samples[[cl]], collapse = ','), "'",
                               " -e ", opt$ess,
                               " -n ", opt$noness,
                               " -o ", bagel2_gene_out)

    # Set LSF bsub command (gene)
    bagel_lsf_cmds[[paste0(cl, '_gene')]] <- paste('bsub -J', lfs_jobname,
                                                    lsf_mem,
                                                    '-o', lfs_out,
                                                    '-e', lfs_err,
                                                    '-q', lfs_queue,
                                                    bagel2_gene_cmd)

    # Set BAGEL2 sgrna command
    bagel2_sgrna_out <- file.path(result_dir, 'BAGEL2.sgrna.bf')
    bagel2_sgrna_cmd <- paste0("BAGEL.py bf",
                               ' -i ', bagel_lfc,
                               " -c '", paste(cell_line_samples[[cl]], collapse = ','), "'",
                               " -e ", opt$ess,
                               " -n ", opt$noness,
                               " -r ",
                               " -o ", bagel2_sgrna_out)

    # Set LSF bsub command (sgrna)
    bagel_lsf_cmds[[paste0(cl, '_sgrna')]] <- paste('bsub -J', lfs_jobname,
                                                     lsf_mem,
                                                     '-o', lfs_out,
                                                     '-e', lfs_err,
                                                     '-q', lfs_queue,
                                                     bagel2_sgrna_cmd)
          }

  # Write commands to file
  message(paste("Writing BAGEL2 LSF commands for:", d))
  bagel_cmds_path <- file.path(opt$dir, 'SCRIPTS', 'single_guide', paste0('BAGEL2_bsub_commands_', d, '.sh'))
  write_lines(bagel_lsf_cmds, file(bagel_cmds_path), append = F)
  message(paste("BAGEL2 LSF commands written to:", bagel_cmds_path))
}

message('Done.')
