# This script compares dual guide library hit genes to those expressed in normal tissues in GTEX and not expressed in tumours in PCAWG

# Load libraries ----------------------------------------------------------


library(tidyverse)
library(VennDiagram)
library(RColorBrewer)


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

screen_results <- read.table(file.path(top_dir, 'DATA','postprocessing','combined_gene_level_results.binary.tsv'),sep = '\t',header = TRUE)


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


# Expression of screened genes -------------------------------------------------------


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


# Counting number of patients where gene is expressed or not expressed
by_patient.subset_long <- by_patient.subset %>% pivot_longer(cols=!Gene.ID & !Gene.Name,names_to = "patient",values_to = "TPM")

by_patient.subset_counts <- by_patient.subset_long %>% 
  group_by( Gene.ID,Gene.Name ) %>% 
  summarise( n=n(),
             not_expressed = sum( TPM < 1 ),
             pct_not_expressed = round( ( not_expressed / n ) * 100,  2 ) , 
             expressed = sum( TPM >= 1 ),
             pct_expressed = round( ( expressed / n ) * 100, 2 ) )

# Find frequencies of expressed/not expressed in particular tissue types. 
by_patient.subset_long <- by_patient.subset_long %>% 
  separate(patient, sep = "\\.\\.\\.", into = c( 'label', 'tumour' ), remove = F ) %>%
  separate(label, sep = "\\.\\.", into = c( 'tumour_type', 'tissue_type','patientID' ), remove = F ) %>%
  select(-c(label))
by_patient.subset_long_tumour_only <- by_patient.subset_long %>%
  filter(tumour != "normal")

# Count number of times gene is expressed or not expressed in patient tumours by type
by_patient.subset_counts_tumour_type <- by_patient.subset_long_tumour_only %>% 
  group_by( Gene.ID,Gene.Name,tumour_type,tissue_type) %>% 
  summarise( n=n(),
             not_expressed = sum( TPM < 1 ),
             pct_not_expressed = round( ( not_expressed / n ) * 100,  2 ) , 
             expressed = sum( TPM >= 1 ),
             pct_expressed = round( ( expressed / n ) * 100, 2 ) )


# Expression by tumour type ------------------------------------------------------------


# Gene pair TPMs including cancer type
gene_pair_TPMs_by_patient_w_tumour_types <- data.frame()
for ( i in 1:length( paired_library.pairs$sorted_gene_pair ) ) {
  gene_pair <- as.vector( paired_library.pairs$sorted_gene_pair[i] )
  geneA <- paired_library.pairs$l_gene[i]
  geneB <- paired_library.pairs$r_gene[i]
  
  if ( !geneA %in% missing_genes & !geneB %in% missing_genes ) {
    tmp.df <- by_patient.subset_long_tumour_only %>% 
      filter( Gene.Name == geneA | Gene.Name == geneB ) %>% 
      select( Gene.Name,patient,tumour_type,tissue_type,tumour,patientID,TPM ) %>% 
      spread( Gene.Name, TPM ) %>% 
      mutate( 'sorted_pair_id' = gene_pair,
              'geneA.tumour.TPM' = !!sym(geneA),
              'geneB.tumour.TPM' = !!sym(geneB),
              'geneA'=geneA,
              'geneB'=geneB,
              class = case_when(  !!sym(geneA) < 1 & !!sym(geneB) >= 1  ~ "not_expressed:expressed",
                                  !!sym(geneA) >= 1 & !!sym(geneB) < 1 ~ "expressed:not_expressed",
                                  !!sym(geneA) >= 1 & !!sym(geneB) >= 1 ~ "expressed:expressed",
                                  !!sym(geneA) < 1 & !!sym(geneB) < 1 ~ "not_expressed:not_expressed",
                                  TRUE ~ as.character( 'unclassified' ) ) ) %>%
      select( patient,patientID,sorted_pair_id,geneA,geneB,tumour_type,tissue_type,tumour,geneA.tumour.TPM, geneB.tumour.TPM,class)
    
    if ( nrow( gene_pair_TPMs_by_patient_w_tumour_types ) == 0 ) {
      gene_pair_TPMs_by_patient_w_tumour_types <- tmp.df
    } else {
      gene_pair_TPMs_by_patient_w_tumour_types <- rbind( gene_pair_TPMs_by_patient_w_tumour_types, tmp.df )
    }
  }
}


# Finding frequency of expression of gene pairs in different tumour types
gene_pair_TPMs_by_tumour_type.summary <- gene_pair_TPMs_by_patient_w_tumour_types %>% 
  group_by( sorted_pair_id, tumour_type,class) %>% 
  summarise( n_patients = n() ) %>% 
  spread( class, n_patients ) %>%
  replace(is.na(.), 0) %>%
  ungroup() %>%
  mutate( total_patients = rowSums( .[3:ncol(.)] ),
          pct_geneA_hits = round( ( get( 'expressed:not_expressed') / total_patients ) * 100 , 2 ),
          pct_geneB_hits = round( ( get( 'not_expressed:expressed') / total_patients ) * 100 , 2 ),
          pct_geneA_and_geneB_hits = round( pct_geneA_hits + pct_geneB_hits, 2 ) ) 

# Frequencies of expressed and not expressed in particular tumour types
gene_pair_TPMs_by_tumour_type.summary <- gene_pair_TPMs_by_tumour_type.summary %>%
  left_join(gene_pair_tumour_normal_TPMs,by=c('tumour_type','sorted_pair_id')) %>%
  select(sorted_pair_id,tumour_type,tissue_type,geneA,geneB,`expressed:expressed`,`expressed:not_expressed`,`not_expressed:expressed`,`not_expressed:not_expressed`,total_patients,pct_geneA_hits,pct_geneB_hits,pct_geneA_and_geneB_hits,geneA.tumour.TPM,geneB.tumour.TPM,geneA.normal.TPM,geneB.normal.TPM)

# Define ubiquitously expressed as both genes TPM > 1 in all tissue types measured in GTEX 
# Filter instances of one or more genes not expressed
pcawg_gtex_no_expression_normal <- gene_pair_TPMs_by_tumour_type.summary %>%
  filter(geneA.normal.TPM<1|geneB.normal.TPM<1)

# List of genes that have tissues where at least one is not expressed
lost_expression_normal <- unique(pcawg_gtex_no_expression_normal$sorted_pair_id)

# List of gene pairs that are ubiquitously expressed in normal tissue
ubiquitously_expressed_normal <- gene_pair_TPMs_by_tumour_type.summary %>%
  filter(!sorted_pair_id %in% lost_expression_normal)
ubiquitously_expressed_normal<- unique(ubiquitously_expressed_normal$sorted_pair_id)

# Pct_gene_A_and_gene_B_hits: > 0 
# Percentage of patients of particular tumour type that have loss of expression of one member of the pair, but not the other
pcawg_gtex_loss_expression_tumour <- gene_pair_TPMs_by_tumour_type.summary %>%
  filter(pct_geneA_and_geneB_hits>0)

# Lost expression in tumour: loss of expression of one member of the pair in a tumour type 
lost_expression_tumour <- unique(pcawg_gtex_loss_expression_tumour$sorted_pair_id)


# Plot venn diagram  ------------------------------------------------------------

# Finding overlap of:
# - screen hit gene pairs in at least one cell line
# - gene pairs with loss of expression of one gene in a tumour
# - and both genes ubiquitously expressed in normal tissue

# Identify genes that are hits in at least one cell line
hit_matrix <- screen_results |> select(sorted_gene_pair, bassik__Lung_NSCLC, bassik__Melanoma, bassik__Pancreas, bassik__total)
hit_matrix_binary <- hit_matrix |> mutate("hits_binary" = case_when(bassik__total >0 ~ 1, .default = 0))
hit_matrix_pairs <- hit_matrix_binary %>%
  mutate(sorted_gene_pair = str_replace(sorted_gene_pair,"\\|","_"))
screen_hits <- hit_matrix_pairs %>%
  filter(hits_binary==1)

myCol <- brewer.pal(3, "Pastel2")

pdf(file.path(output_plot_dir, 'F7D__screen_hits_GTEX_expression_tumour_normal.pdf'), width = 4.8, height = 4.8)

grid.newpage()
venn.plot <- venn.diagram(
  x = list(
    lost_expression_tumour,
    ubiquitously_expressed_normal,
    screen_hits$sorted_gene_pair
  ),
  category.names = c("Lost expression\n in tumour", 
                     "Ubiquitously \nexpressed normal", 
                     "Screen hits"),
  filename = NULL,  
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  fontfamily = "sans",
  cat.fontfamily = "sans",
  cat.just = rep(list(c(0.7, 0.7)), 3),
  cat.pos = c(-20, 27, 135),
  cat.dist = c(0.05, 0.07, 0.06)
)
grid.draw(venn.plot) 
dev.off()  


