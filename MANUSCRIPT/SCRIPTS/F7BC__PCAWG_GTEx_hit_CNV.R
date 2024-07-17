# This script reads expression data, filters for top hits, processes the data, and plots heatmaps of TPM values


# Load libraries ----------------------------------------------------------

library(tidyverse)

# Read in input files -----------------------------------------------------

# Set top level directory
top_dir <- getwd()

# Change this path for scripts to run!
expression_data <- read_delim("PCAWG_GTEx_tumour_normal_gene_pair_TPMs.tsv", delim = "\t", trim_ws = TRUE)

# Change this path for output directory!
output_plot_dir <- file.path(top_dir, 'MANUSCRIPT', 'PLOTS')

# Process data ------------------------------------------------------------

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