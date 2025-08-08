# This script reads expression data, filters for top hits, processes the data, and plots heatmaps of TPM values


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(stringr)


# Read in input files -----------------------------------------------------

# Set top level directory
top_dir <- getwd()

# Normalised expression data per patient from https://www.ebi.ac.uk/gxa/experiments/E-MTAB-5423/Results
by_patient <- read.table(file.path(top_dir, 'E-MTAB-5423-query-results.tpms.tsv'),sep = '\t',header = TRUE)

# Expression averages across tumour and tissue types https://www.ebi.ac.uk/gxa/experiments/E-MTAB-5200/Results
# Including normal tissue GTEx data for comparison 
by_tumour_type <- read.table(file.path(top_dir, 'E-MTAB-5200-query-results.tpms.tsv'),sep = '\t',header = TRUE)

output_plot_dir <- file.path(top_dir, 'MANUSCRIPT', 'PLOTS')

annotated_library <- read.table(file.path(top_dir, 'METADATA','libraries','paralog_library.tsv'),sep = '\t',header = TRUE)


# Process data ------------------------------------------------------------

#Identifying gene pairs in the library
paired_library.pairs <- annotated_library %>% 
  filter( guide_type == 'gene|gene') %>% 
  select( sorted_gene_pair ) %>% 
  unique() %>%
  separate(sorted_gene_pair, sep = "\\|", into = c( 'l_gene', 'r_gene' ), remove = F )

paired_library.pairs$sorted_gene_pair <- sub("\\|","_",paired_library.pairs$sorted_gene_pair)

#Replace NA with 0 
by_patient <- by_patient %>% replace(is.na(.), 0)

by_tumour_type <- by_tumour_type %>% replace(is.na(.), 0)

#Limit expression data to gene pairs in the screen 
by_tumour_type.subset <- by_tumour_type %>% 
  filter(Gene.Name %in% paired_library.pairs$l_gene | Gene.Name %in% paired_library.pairs$r_gene )

by_tumour_type.subset_long <- by_tumour_type.subset %>% pivot_longer(cols=!Gene.ID & !Gene.Name,names_to = "label",values_to = "TPM")

by_patient.subset <- by_patient %>% 
  filter(Gene.Name %in% paired_library.pairs$l_gene | Gene.Name %in% paired_library.pairs$r_gene )
missing_genes <- setdiff( c( paired_library.pairs$l_gene, paired_library.pairs$r_gene), by_patient.subset$Gene.Name )

# Labelling tumours vs normal
by_tumour_type.subset_long <- by_tumour_type.subset_long %>% 
  mutate(class=case_when(
    str_detect(label,"GTEx") ~ "normal - GTEX",
    str_detect(label,"adjacent") ~ "normal - tumour adjacent",
    TRUE ~ "tumour"
  ))

#Adding the tumour type
tumour_TPMs <- by_tumour_type.subset_long %>% 
  filter( class == 'tumour') %>% 
  separate(label, sep = "\\.\\.", into = c( 'tumour_type', 'tissue_type' ), remove = F )

#Adding tissue type and filter out normal - tumour adjacent
normal_TPMs <- by_tumour_type.subset_long %>% 
  filter( class == 'normal - GTEX') %>% 
  separate(label, sep = "\\.\\.\\.", into = c( 'normal', 'GTEx','tissue_type' ), remove = F ) %>%
  select(-c(normal,GTEx))

#Averaging across the expression values for different regions of the same tissue 
normal_TPM_averages <- normal_TPMs %>%
  group_by(Gene.ID,Gene.Name,tissue_type) %>% 
  summarise( average_normal_TPM = mean(TPM) )

# Comparing TPMs in the normal with TPMs in tumour
tumour_normal_TPMs <- tumour_TPMs %>%
  left_join(normal_TPM_averages,by=c('Gene.ID','Gene.Name','tissue_type'))

# Checking expression of each gene in tumour and normal - gene A tumour TPM gene A normal TPM 
gene_pair_tumour_normal_TPMs <- data.frame()
for ( i in 1:length( paired_library.pairs$sorted_gene_pair ) ) {
  gene_pair <- as.vector( paired_library.pairs$sorted_gene_pair[i] )
  geneA <- paired_library.pairs$l_gene[i]
  geneB <- paired_library.pairs$r_gene[i]
  
  if ( !geneA %in% missing_genes & !geneB %in% missing_genes ) {
    tmp.df <- tumour_normal_TPMs %>% 
      filter( Gene.Name == geneA | Gene.Name == geneB ) %>% 
      select( Gene.Name, tumour_type, tissue_type, TPM) %>% 
      spread( Gene.Name, TPM ) %>%
      mutate( 'sorted_pair_id' = gene_pair,
              'geneA.tumour.TPM' = get( geneA ),
              'geneB.tumour.TPM' = get( geneB ),
              'geneA'=geneA,
              'geneB'=geneB) %>%
      select(tumour_type,tissue_type, sorted_pair_id,geneA,geneB, geneA.tumour.TPM, geneB.tumour.TPM)
    if ( nrow( gene_pair_tumour_normal_TPMs ) == 0 ) {
      gene_pair_tumour_normal_TPMs <- tmp.df
    } else {
      gene_pair_tumour_normal_TPMs<- rbind( gene_pair_tumour_normal_TPMs, tmp.df )
    }
  }
}                         

gene_pair_normal_TPMs <- data.frame()
for ( i in 1:length( paired_library.pairs$sorted_gene_pair ) ) {
  gene_pair <- as.vector( paired_library.pairs$sorted_gene_pair[i] )
  geneA <- paired_library.pairs$l_gene[i]
  geneB <- paired_library.pairs$r_gene[i]
  
  if ( !geneA %in% missing_genes & !geneB %in% missing_genes ) { 
    tmp.df <- tumour_normal_TPMs %>% 
      filter( Gene.Name == geneA | Gene.Name == geneB ) %>% 
      select( Gene.Name, tumour_type, tissue_type, average_normal_TPM) %>% 
      spread( Gene.Name, average_normal_TPM ) %>% 
      mutate( 'sorted_pair_id' = gene_pair,
              'geneA.normal.TPM' = get( geneA ),
              'geneB.normal.TPM' = get( geneB ),
              'geneA'=geneA,
              'geneB'=geneB) %>%
      select(tumour_type,tissue_type, sorted_pair_id,geneA,geneB, geneA.normal.TPM, geneB.normal.TPM)
    if ( nrow( gene_pair_normal_TPMs ) == 0 ) {
      gene_pair_normal_TPMs <- tmp.df 
    } else {
      gene_pair_normal_TPMs<- rbind( gene_pair_normal_TPMs, tmp.df )
    }
  }
}         

gene_pair_tumour_normal_TPMs <- gene_pair_tumour_normal_TPMs %>%
  left_join(gene_pair_normal_TPMs,by=c('tumour_type','tissue_type','sorted_pair_id','geneA','geneB'))

expression_data <- gene_pair_tumour_normal_TPMs

# Select only top hits
top_hits <- expression_data |> 
  filter(sorted_pair_id %in% c('CNOT7_CNOT8', 'GDI1_GDI2', 'SAR1A_SAR1B', 'SEC23A_SEC23B', 
                               'ASF1A_ASF1B', 'PDS5A_PDS5B', 'SLC25A28_SLC25A37', 'CCNL1_CCNL2'))

# Separate into Gene A and Gene B and combine back
merged_data <- bind_rows(
  top_hits |> select(tumour_type, tissue_type, sorted_pair_id, gene = geneA, normal = geneA.normal.TPM, tumour = geneA.tumour.TPM),
  top_hits |> select(tumour_type, tissue_type, sorted_pair_id, gene = geneB, normal = geneB.normal.TPM, tumour = geneB.tumour.TPM)
)

# Adjust TPMs -------------------------------------------------------------

# Define function for TPM adjustment
adjust_TPM <- function(df, col_name) {
  df |> mutate(
    !!paste0(col_name, "_TPM") := case_when(
      !!sym(col_name) < 1 ~ '0-1',
      !!sym(col_name) >= 1 & !!sym(col_name) < 10 ~ '1-10',
      !!sym(col_name) >= 10 & !!sym(col_name) < 20 ~ '10-20',
      !!sym(col_name) >= 20 & !!sym(col_name) < 40 ~ '20-40',
      !!sym(col_name) >= 40 & !!sym(col_name) < 60 ~ '40-60',
      !!sym(col_name) >= 60 & !!sym(col_name) < 80 ~ '60-80',
      !!sym(col_name) >= 80 & !!sym(col_name) < 100 ~ '80-100',
      !!sym(col_name) >= 100 & !!sym(col_name) < 200 ~ '100-200',
      !!sym(col_name) >= 200 & !!sym(col_name) < 300 ~ '200-300',
      !!sym(col_name) >= 300 & !!sym(col_name) < 400 ~ '300-400',
      !!sym(col_name) >= 400 & !!sym(col_name) < 500 ~ '400-500',
      !!sym(col_name) >= 500 & !!sym(col_name) < 600 ~ '500-600',
      !!sym(col_name) >= 600 & !!sym(col_name) < 700 ~ '600-700'
    )
  )
}

# Adjust TPM for normal and tumour
adjusted_data <- merged_data |> adjust_TPM("normal") |> adjust_TPM("tumour")


# Plot heatmaps -----------------------------------------------------------

plot_TPM <- function(data, tpm_column, title, y_axis, colors) {
  ggplot(data, aes_string(y = y_axis, x = "gene", fill = tpm_column)) + 
    geom_tile(color = "gray30", aes(width = 0.9, height = 0.9), linewidth = 0.4) + 
    scale_fill_manual(
      breaks = c('0-1', '1-10', '10-20', '20-40', '40-60', '60-80', '80-100', '100-200', '200-300', '300-400', '400-500', '500-600', '600-700'),
      values = colors
    ) +
    coord_cartesian(clip = "off") + 
    theme_classic() + 
    labs(x = '', y = '', title = title) +
    theme(
      plot.margin = unit(c(1, 1, 1, 1), "lines"), 
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8, face = "italic"),
      axis.text.y = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 12),
      legend.position = "right",
      legend.justification = "left",
      legend.box.spacing = unit(0, "pt")
    ) + 
    guides(fill = guide_legend(title = "TPM")) +
    geom_vline(xintercept = seq(2.5, 14.5, by = 2), linetype = "dashed", colour = "gray30", linewidth = 0.5)
}

# Define color scales for normal and tumour plots
colour_palette <- c("#B71B1BFF", "#E3F2FDFF", "#BADEFAFF", "#90CAF8FF", "#64B4F6FF", "#41A5F4FF", "#2096F2FF", "#1E87E5FF", "#1976D2FF", "#1465BFFF", "#0C46A0FF", "#FFF49DFF", "#FABF2CFF", "#F8A725FF")

# Plot Normal data --------------------------------------------------------

normal_heatmap <- plot_TPM(adjusted_data, "normal_TPM", 'Normal Tissue - GTEx data', "tissue_type", colour_palette)


# Plot Tumour data --------------------------------------------------------

tumour_heatmap <- plot_TPM(adjusted_data, "tumour_TPM", 'Tumour Tissue - TCGA data', "tumour_type", colour_palette)

# Save plots --------------------------------------------------------------

ggsave(file.path(output_plot_dir, 'F7B__top_hits_GTEX_expression_tumour.png'), tumour_heatmap, dpi = 300, width = 8, height = 6)
ggsave(file.path(output_plot_dir, 'F7C__top_hits_GTEX_expression_normal.png'), normal_heatmap, dpi = 300, width = 8, height = 6)
