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
# Download CMP copy number                                 #
############################################################

message('Downloading CMP WES copy number...')

# Read in Cell Model Passport gene copy number per cell line (20220623)
tmpfile <- tempfile()
download.file(url = "https://cog.sanger.ac.uk/cmp/download/WES_pureCN_CNV_genes_20220623.zip", destfile = tmpfile, quiet = T, mode = 'wb')
conn <- unz(tmpfile, 'WES_pureCN_CNV_genes_total_copy_number_20220623.csv')
cmp_cn <- read.csv(conn, skip = 1)
colnames(cmp_cn)[1] <- 'symbol'
cmp_cn <- cmp_cn[-c(1:2),]
unlink(tmpfile)

############################################################
# Download CMP metadata.                                   #
############################################################

message('Downloading CMP model metadata...')

# Read in model metadata
tmpfile <- tempfile()
download.file(url = 'https://cog.sanger.ac.uk/cmp/download/model_list_20220205.csv', destfile = tmpfile, quiet = T)
cmp_meta <- read.delim(tmpfile, sep = ",", header = T)

############################################################
# Filter metadata for screened lines                       #
############################################################

message('Filtering model metadata...')

# Remove lines not included in our screen from metadata
screened_models <- cmp_meta %>%
  filter(BROAD_ID %in% sample_mapping$depMapID) %>%
  select(model_id, 'depMapID' = 'BROAD_ID')

############################################################
# Add to expanded library.                                 #
############################################################

message('Get gene pairs from library...')

# Reduce guide level library down to just gene pairs
expanded_library <- expanded_library %>%
  select(sorted_gene_pair, targetA, targetB) %>%
  unique()

############################################################
# Filtering CN for screened lines and genes                #
############################################################

message('Removing unwanted lines and genes from CN...')

# Remove lines and genes not included in our screen from CN
cmp_cn.filt <- cmp_cn %>% 
  gather(model_id, total_copy_number, -symbol) %>% # models are column names
  mutate(total_copy_number = as.numeric(total_copy_number)) %>%
  filter(model_id %in% screened_models$model_id) %>% # Filter for screened lines
  filter(symbol %in% c(expanded_library$targetA, expanded_library$targetB)) %>% # Filter for screened genes
  left_join(screened_models, by = 'model_id')

############################################################
# Expand CN by library metadata and cell line              #
############################################################

message('Expanding library by cell line...')

cmp_cn.ann <- expand.grid(expanded_library$sorted_gene_pair, sample_mapping$depMapID) %>%
  unique() %>%
  select('sorted_gene_pair' = 'Var1', 'depMapID' = 'Var2' ) %>%
  filter(!is.na(depMapID))

message('Adding additional library metadata...')

cmp_cn.ann <- cmp_cn.ann %>%
  left_join(expanded_library %>% select(sorted_gene_pair, targetA, targetB), by = 'sorted_gene_pair', relationship = 'many-to-many') %>%
  unique()

message('Adding screened models...')

cmp_cn.ann <- cmp_cn.ann %>%
  left_join(screened_models, by = 'depMapID') %>%
  unique()

message('Adding copy number for targetA...')

cmp_cn.ann <- cmp_cn.ann %>%
  left_join(cmp_cn.filt %>% 
              select(model_id, depMapID, 'gene' = symbol, total_copy_number) %>%
              rename_at(vars(-gene, -depMapID, -model_id), ~ paste('targetA', ., sep = '__')), by = c('targetA' = 'gene', 'depMapID', 'model_id')) %>%
  unique()

message('Adding copy number for targetB...')

cmp_cn.ann <- cmp_cn.ann %>%
  left_join(cmp_cn.filt %>% 
              select(model_id, depMapID, 'gene' = symbol, total_copy_number) %>%
              rename_at(vars(-gene, -depMapID, -model_id), ~ paste('targetB', ., sep = '__')), by = c('targetB' = 'gene', 'depMapID', 'model_id')) %>%
  unique()

message('Adding sample mapping...')

cmp_cn.ann <- cmp_cn.ann %>%
  left_join(sample_mapping %>% select(cell_line_label, cancer_type, depMapID), by = 'depMapID', relationship = 'many-to-many') %>%
  unique()
              
message('Postprocessong...')

cmp_cn.ann <- cmp_cn.ann %>%
  rename('CMP_model_id' = 'model_id') %>%
  select(depMapID, CMP_model_id, cell_line_label, cancer_type, 
         sorted_gene_pair, 
         targetA, targetB, 
         targetA__total_copy_number, 
         targetB__total_copy_number)

message('Plot DepMap copy number...')

cmp_cn.ann.narrow <- rbind(cmp_cn.ann %>% select(cell_line_label, 'gene' = targetA, 'cn' = targetA__total_copy_number), 
                           cmp_cn.ann %>% select(cell_line_label, 'gene' = targetB, 'cn' = targetB__total_copy_number)) 
cmp_cn.ann.narrow <- cmp_cn.ann.narrow %>% unique()

cmp_cn_plot <- 
  ggplot(cmp_cn.ann.narrow, aes(x = gene, y = cn, label = paste(cell_line_label, gene, sep = " : "))) +
    geom_point(size = 0.5) +
    geom_text(data = subset(cmp_cn.ann.narrow, cn > 15), nudge_y = 1) +
    scale_y_continuous(breaks = pretty_breaks(10)) +
    theme_pubr(base_size = 14) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    labs(x = '', y = 'Cell Model Passport total copy number')
cmp_cn_plot_path = file.path(opt$dir, 'DATA', 'postprocessing', 'CMP_total_copy_number.png')
ggsave(filename = cmp_cn_plot_path, plot = cmp_cn_plot, device = 'png', dpi = 300 , width = 4000, height = 3000, units = 'px')


############################################################
# Write outputs to file                                    #
############################################################

message('Writing to file...')

# Save as a TSV
cn.filepath <- file.path(opt$dir, 'DATA', 'postprocessing', 'intermediate_tables', 'cmp_total_copy_number.tsv')
write.table(cmp_cn.ann, file = cn.filepath, sep = "\t", row.names = F, quote = F)
message(paste0('Cell model passport copy number TSV written to:', cn.filepath))

# Save as an Rdata object
cn.rds.filepath <- file.path(opt$dir, 'DATA', 'RDS', 'postprocessing', 'cmp_total_copy_number.rds')
saveRDS(cmp_cn.ann, file = cn.rds.filepath)
message(paste0('Cell model passport copy number RDS written to:', cn.rds.filepath))

message('Done.')