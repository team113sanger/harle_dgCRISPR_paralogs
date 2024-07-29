# Purpose: visualisation of PanCancer/PCAWG clinical and CNV data for top genes

# Load libraries ----------------------------------------------------------
library(tidyverse)

# Read in input files -----------------------------------------------------

# Set top level directory
top_dir <- getwd()

# Read directly from Figshare
clinical <- read.delim('https://figshare.com/ndownloader/files/47993605', header = TRUE)
cnv <- read.delim('https://figshare.com/ndownloader/files/47993560', header = TRUE)

# Change this path for output directory!
output_plot_dir <- file.path(top_dir, 'MANUSCRIPT', 'PLOTS')

# Merge data --------------------------------------------------------------
cnv_clinical <- merge(clinical, cnv, by.x = 'Sample.ID', by.y = 'SAMPLE_ID')

# Pivot data --------------------------------------------------------------
cnv_pivot <- cnv_clinical |> 
  select(Sample.ID, Cancer.Type, CNOT7:CCNL2) |>
  pivot_longer(cols = CNOT7:CCNL2, names_to = 'Gene', values_to = 'CNV')

# Count patients by cancer type -------------------------------------------
cancer_type_count <- cnv_clinical |> 
  count(Cancer.Type) |> 
  rename(Number_of_patients = n)

# Count CNV occurrences by cancer type and gene ---------------------------
genes_count <- cnv_pivot |> 
  group_by(Cancer.Type, Gene) |> 
  count(CNV) |> 
  rename(CNV_gene = n) |> 
  ungroup()

# Merge counts ------------------------------------------------------------
merged_counts <- merge(cancer_type_count, genes_count)

# Calculate percentages and categorise CNV values -------------------------
percentage_cnv <- merged_counts |> 
  mutate(Percentage = (CNV_gene / Number_of_patients) * 100) |>
  mutate(Copy_Number = case_when(CNV == -2 ~ '-2',
                                 CNV == -1 ~ '-1',
                                 CNV == 0 ~ '0',
                                 CNV == 1 ~ '1',
                                 CNV == 2 ~ '2')) |> 
  ungroup() |> 
  drop_na()

# Create plot -------------------------------------------------------------
cnv_plot <- ggplot(percentage_cnv, aes(x = Cancer.Type, y = Percentage, fill = Copy_Number)) + 
  geom_bar(position = "stack", stat = "identity", width = 0.7) +  
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14), 
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14), 
        plot.title = element_text(size = 14),
        panel.background = element_rect(fill = "white")) + 
  labs(title = '', y = 'Percentage of Patients', x = '') +
  guides(fill = guide_legend(title = "CNV")) +
  scale_fill_manual(breaks = c('2','1','0','-1','-2'),
                    values = c("#B71B1BFF","#FFCCD2FF","#C7E5C9FF", "#90CAF8FF", "#0C46A0FF")) + 
  facet_wrap("Gene", ncol = 2) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "lines")) +
  theme(strip.text = element_text(face = "italic"))

# Save plot ---------------------------------------------------------------
ggsave(file.path(output_plot_dir, 'F7A__CNV_for_top_hits.png'), cnv_plot, dpi = 300, width = 14, height = 20)
