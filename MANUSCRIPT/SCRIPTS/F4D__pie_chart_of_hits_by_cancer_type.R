## Pie chart of hits by cancer type

# Prepare data ------------------------------------------------------------

# Source data preparation script
source(file.path('MANUSCRIPT', 'SCRIPTS', 'prepare_data.R'))

# Load required libraries
library(tidyverse)
library(ggVennDiagram)

# Select relevant columns from the data
binary_results <- binary_results |>
  select(sorted_gene_pair, Lung, Melanoma, Pancreatic)

# Create lists of gene pairs for each condition
conditions <- list(
  Lung = binary_results %>% filter(Lung >= 1) %>% pull(sorted_gene_pair),
  Melanoma = binary_results %>% filter(Melanoma >= 1) %>% pull(sorted_gene_pair),
  Pancreas = binary_results %>% filter(Pancreatic >= 1) %>% pull(sorted_gene_pair)
)

# Create lists of gene pairs for each condition and prepare data for ggVennDiagram
venn_data <- map(c("Lung", "Melanoma", "Pancreatic"), ~ binary_results %>%
                   filter(get(.x) >= 1) %>%
                   pull(sorted_gene_pair) %>%
                   as.character()) %>%
  set_names(c("Lung", "Melanoma", "Pancreas"))

# Draw the Venn diagram
# Note names are coloured and not black
venn.plot <- ggVennDiagram(venn_data, 
                           label_alpha = 0.5, 
                           label = 'count',
                           set_color = c("#0072b2", "#e69f00", "#009e73"))+
  scale_fill_gradient(low = "#caddee", high = "#305294") +
  theme(legend.position = "none") 

# Save the Venn diagram 
ggsave(file.path(output_plot_dir, 'F4D__pie_chart_of_hits_by_cancer_type.png'), venn.plot, dpi = 300, width = 12, height = 10)
