suppressPackageStartupMessages(suppressWarnings(library(optparse)))
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))

############################################################
# OPTIONS                                                  #
############################################################

option_list = list(
  make_option(c("-d", "--dir"), type = "character",
              help = "repository directory path", metavar = "character"),
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
# Expanded library                                         #
############################################################

check_file_exists(opt$annotation)

message(paste("Reading expanded library:", opt$annotation))

# Read in library
expanded_library <- read.delim(opt$annotation, sep = "\t", header = T)

############################################################
# Initialise list                                          #
############################################################

# Initialise empty list
common_essentiality_gene_lists <- list()

############################################################
# Download BAGEL common essentials                         #
############################################################

message('Downloading BAGEL essentials...')

# Read in BAGEL2 common essentials (c6af217)
tmpfile <- tempfile()
download.file(url = 'https://raw.githubusercontent.com/hart-lab/bagel/master/CEGv2.txt', destfile = tmpfile, quiet = T)
common_essentiality_gene_lists[['BAGEL2 common essential']][['label']] <- 'is_common_essential_bagel2'
common_essentiality_gene_lists[['BAGEL2 common essential']][['gene']] <- read.delim(tmpfile, sep = "\t", header = T)$GENE
common_essentiality_gene_lists[['BAGEL2 common essential']][['n_genes']] <- length(common_essentiality_gene_lists[['BAGEL2 common essential']][['gene']])

############################################################
# Download BAGEL common non-essentials                     #
############################################################

message('Downloading BAGEL non-essentials...')

# Read in BAGEL2 common non-essentials (c6af217)
tmpfile <- tempfile()
download.file(url = 'https://raw.githubusercontent.com/hart-lab/bagel/master/NEGv1.txt', destfile = tmpfile, quiet = T)
common_essentiality_gene_lists[['BAGEL2 common non-essential']][['label']] <- 'is_common_nonessential_bagel2'
common_essentiality_gene_lists[['BAGEL2 common non-essential']][['gene']] <- read.delim(tmpfile, sep = "\t", header = T)$GENE
common_essentiality_gene_lists[['BAGEL2 common non-essential']][['n_genes']] <- length(common_essentiality_gene_lists[['BAGEL2 common non-essential']][['gene']])

############################################################
# Download DepMap common essentials                        #
############################################################

message('Downloading DepMap common essentials (ceres)...')

# Read in DepMap 22Q2 common essentials (ceres, pre_chronos)
tmpfile <- tempfile()
download.file(url = 'https://ndownloader.figshare.com/files/34990024', destfile = tmpfile, quiet = T)
common_essentiality_gene_lists[['DepMap (CERES) common essential']][['label']] <- 'is_common_essential_depmap_ceres'
common_essentiality_gene_lists[['DepMap (CERES) common essential']][['gene']] <- read.delim(tmpfile, sep = " ", col.names = c('gene', 'entrez_id'), header = F, skip = 1) %>% mutate(entrez_id = str_replace_all(entrez_id, "[()]", "")) %>% pull(gene)
common_essentiality_gene_lists[['DepMap (CERES) common essential']][['n_genes']] <- length(common_essentiality_gene_lists[['DepMap (CERES) common essential']][['gene']])

message('Downloading DepMap common essentials (chronos)...')

# Read in DepMap 22Q2 common essentials (post_chronos)
tmpfile <- tempfile()
download.file(url = 'https://ndownloader.figshare.com/files/34990027', destfile = tmpfile, quiet = T)
common_essentiality_gene_lists[['DepMap (Chronos) common essential']][['label']] <- 'is_common_essential_depmap_chronos'
common_essentiality_gene_lists[['DepMap (Chronos) common essential']][['gene']] <- read.delim(tmpfile, sep = " ", col.names = c('gene', 'entrez_id'), header = F, skip = 1) %>% mutate(entrez_id = str_replace_all(entrez_id, "[()]", "")) %>% pull(gene)
common_essentiality_gene_lists[['DepMap (Chronos) common essential']][['n_genes']] <- length(common_essentiality_gene_lists[['DepMap (Chronos) common essential']][['gene']])

############################################################
# Download DepMap common non-essentials                    #
############################################################

message('Downloading DepMap common non-essentials (ceres)...')

# Read in DepMap 22Q2 common non-essentials (ceres, pre_chronos)
tmpfile <- tempfile()
download.file(url = 'https://ndownloader.figshare.com/files/34990051', destfile = tmpfile, quiet = T)
common_essentiality_gene_lists[['DepMap (CERES) common non-essential']][['label']] <- 'is_common_nonessential_depmap_ceres'
common_essentiality_gene_lists[['DepMap (CERES) common non-essential']][['gene']] <- read.delim(tmpfile, sep = " ", col.names = c('gene', 'entrez_id'), header = F, skip = 1) %>% mutate(entrez_id = str_replace_all(entrez_id, "[()]", "")) %>% pull(gene)
common_essentiality_gene_lists[['DepMap (CERES) common non-essential']][['n_genes']] <- length(common_essentiality_gene_lists[['DepMap (CERES) common non-essential']][['gene']])

############################################################
# Map essentials to library                                #
############################################################

message ('Mapping essentiality to library...')

# Get required fields from library 
common_essentiality <- expanded_library %>% 
  select(sorted_gene_pair, targetA, targetB) %>%
  unique()

# Loop over gene lists and add common gene list flags to library-----------
for (gln in names(common_essentiality_gene_lists)){ # Loop over common gene lists
  l <- common_essentiality_gene_lists[[gln]][['label']] # Get column header
  g <- common_essentiality_gene_lists[[gln]][['gene']] # Get gene lists
  
  # Set flags for each gene list per gene pair (1 = present, 0 = absent)
  a <- paste('targetA', l, sep = '__')
  b <- paste('targetB', l, sep = '__')
  common_essentiality <- common_essentiality %>%
    mutate(!!a := ifelse(targetA %in% g, 1, 0), 
           !!b := ifelse(targetB %in% g, 1, 0))
}

############################################################
# Write outputs to file                                    #
############################################################

message('Writing to file...')

# Save as a TSV
common_essentiality.filepath <- file.path(opt$dir, 'DATA', 'postprocessing', 'intermediate_tables', 'depmap_and_bagel_common_genes.tsv')
write.table(common_essentiality, file = common_essentiality.filepath, sep = "\t", row.names = F, quote = F)
message(paste0('Common essential and non-essential TSV written to:', common_essentiality.filepath))

# Save as an Rdata object
common_essentiality.rds.filepath <- file.path(opt$dir, 'DATA', 'RDS', 'postprocessing', 'depmap_and_bagel_common_genes.rds')
saveRDS(common_essentiality, file = common_essentiality.rds.filepath)
message(paste0('Common essential and non-essential RDS written to:', common_essentiality.rds.filepath))

# Save as an Rdata object
common_essentiality_gene_lists.rds.filepath <- file.path(opt$dir, 'DATA', 'RDS', 'postprocessing', 'common_essentiality_gene_lists.rds')
saveRDS(common_essentiality_gene_lists, file = common_essentiality_gene_lists.rds.filepath)
message(paste0('Common essential and non-essential gene lists RDS written to:', common_essentiality_gene_lists.rds.filepath))

message('Done.')