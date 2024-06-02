suppressPackageStartupMessages(suppressWarnings(library(optparse)))
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))
suppressPackageStartupMessages(suppressWarnings(library(scales)))
suppressPackageStartupMessages(suppressWarnings(library(ggpubr)))

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
# Download DepMap relative copy number                     #
############################################################

message('Downloading DepMap relative copy number...')

# Read in DepMap 22Q2 gene copy number per cell line
# Gene level copy number data, log2 transformed with a pseudo count of 1. 
tmpfile <- tempfile()
options(timeout = 500)
download.file(url = 'https://ndownloader.figshare.com/files/34989937', destfile = tmpfile, quiet = T)

message('Reading DepMap relative copy number...')

depmap_cn <- read.delim(tmpfile, sep = ",", header = T)
colnames(depmap_cn)[1] <- 'depMapID'

############################################################
# Prepare DepMap CN                                        #
############################################################

message('Filtering DepMap copy number...')

# Filter to only include screened cell lines 
depmap_cn <- depmap_cn %>% 
  filter(depMapID %in% sample_mapping$depMapID)

# Gather and split gene information
depmap_cn <- depmap_cn %>% 
  gather(gene_label, relative_copy_number_log2ps1, -depMapID) %>%
  separate(gene_label, sep = "\\.\\.", into = c('gene', 'entrez_id')) %>%
  mutate(entrez_id = str_replace_all(entrez_id, "\\.", ""))

############################################################
# Expand CN by library metadata and cell line              #
############################################################

message('Expanding library by cell line...')

depmap_cn.ann <- expand.grid(expanded_library$sorted_gene_pair, sample_mapping$depMapID) %>%
  unique() %>%
  select('sorted_gene_pair' = 'Var1', 'depMapID' = 'Var2' ) %>%
  filter(!is.na(depMapID))

message('Adding additional library metadata...')

depmap_cn.ann <- depmap_cn.ann %>%
  left_join(expanded_library %>% select(sorted_gene_pair, sorted_gene_pair, targetA, targetB), by = 'sorted_gene_pair', relationship = 'many-to-many') %>%
  unique()
  
message('Adding copy number for targetA...')

depmap_cn.ann <- depmap_cn.ann %>%
  left_join(depmap_cn %>% 
              select(-entrez_id) %>%
              rename_at(vars(-gene, -depMapID), ~ paste('targetA', ., sep = '__')), by = c('targetA' = 'gene', 'depMapID')) %>%
  unique()

message('Adding copy number for targetB...')

depmap_cn.ann <- depmap_cn.ann %>%
  left_join(depmap_cn %>% 
              select(-entrez_id) %>%
              rename_at(vars(-gene, -depMapID), ~ paste('targetB', ., sep = '__')), by = c('targetB' = 'gene', 'depMapID')) %>%
  unique()

message('Add cell line metadata to DepMap copy number...')

depmap_cn.ann <- depmap_cn.ann %>%
  left_join(sample_mapping %>% select(depMapID, cell_line_label, cancer_type), by = 'depMapID', relationship = 'many-to-many') %>%
  unique() %>%
  select(depMapID, cell_line_label, cancer_type, 
         sorted_gene_pair, sorted_gene_pair, 
         targetA, targetB, 
         targetA__relative_copy_number_log2ps1, 
         targetB__relative_copy_number_log2ps1)

message('Plot DepMap copy number...')

depmap_cn.ann.narrow <- rbind(depmap_cn.ann %>% select(cell_line_label, 'gene' = targetA, 'cn' = targetA__relative_copy_number_log2ps1), 
                              depmap_cn.ann %>% select(cell_line_label, 'gene' = targetB, 'cn' = targetB__relative_copy_number_log2ps1)) 
depmap_cn.ann.narrow <- depmap_cn.ann.narrow %>% unique()

depmap_cn_plot <- 
  ggplot(depmap_cn.ann.narrow, aes(x = gene, y = cn, label = paste(cell_line_label, gene, sep = " : "))) +
    geom_point(size = 0.8) +
    geom_text(data = subset(depmap_cn.ann.narrow, cn > 2.5), nudge_y = 0.1) +
    geom_hline(yintercept = 1, color = 'red') +
    scale_y_continuous(breaks = pretty_breaks(10)) +
    theme_pubr(base_size = 14) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    labs(x = '', y = 'DepMap relative copy number (log2 +1)')
depmap_cn_plot_path = file.path(opt$dir, 'DATA', 'postprocessing', 'DepMap_relative_copy_number.png')
ggsave(filename = depmap_cn_plot_path, plot = depmap_cn_plot, device = 'png', dpi = 300 , width = 4000, height = 3000, units = 'px')

############################################################
# Write outputs to file                                    #
############################################################

message('Writing to file...')

# Save as a TSV
cn.filepath <- file.path(opt$dir, 'DATA', 'postprocessing', 'intermediate_tables', 'depmap_relative_copy_number.tsv')
write.table(depmap_cn.ann, file = cn.filepath, sep = "\t", row.names = F, quote = F)
message(paste0('DepMap relative copy number TSV written to:', cn.filepath))

# Save as an Rdata object
cn.rds.filepath <- file.path(opt$dir, 'DATA', 'RDS', 'postprocessing', 'depmap_relative_copy_number.rds')
saveRDS(depmap_cn.ann, file = cn.rds.filepath)
message(paste0('DepMap relative copy number RDS written to:', cn.rds.filepath))

message('Done.')