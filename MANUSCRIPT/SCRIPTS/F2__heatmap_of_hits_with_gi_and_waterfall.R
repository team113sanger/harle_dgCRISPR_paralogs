## Heatmap and waterfall plot - plot both landscape and portrait for presentation v publication

# Prepare data ------------------------------------------------------------

# Source data preparation script
source(file.path('MANUSCRIPT', 'SCRIPTS', 'prepare_data.R'))


# Add alpha to colours ----------------------------------------------------

# Function to add alpha to vector of colours (without needing a new library dependency)
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}


# Get hit gene pairs which are hits in 3 or more lines --------------------

# Extract (74) gene pairs which are hits in 3 or more cell lines 
hits <- binary_results |> 
  filter(bassik__total >= 3) |> 
  pull(sorted_gene_pair)


# Get cell lines from full results ----------------------------------------

# Extract non-redundant list of cell line labels
cell_lines <- full_results |>
  select(cell_line_label, cancer_type) |>
  unique()


# Prepare data frame of ordered hits for heatmap --------------------------

# Filter full results to get only hit gene pairs being displayed
# Keep only the columns we need
# Bolt on the total number of cell lines in which pair was a hit binary results 
heatmap_df <- full_results |> 
  filter(sorted_gene_pair %in% hits) |>
  select(sorted_gene_pair, cell_line_label, mean_norm_gi, is_bassik_hit) |> 
  left_join(binary_results |> select(sorted_gene_pair, bassik__total), by = 'sorted_gene_pair')

# Use heatmap data frame to determine the order of the gene pairs in all plots
ordered_gene_pairs <- unique(heatmap_df$sorted_gene_pair[order(heatmap_df$bassik__total)])
heatmap_df$sorted_gene_pair <- factor(heatmap_df$sorted_gene_pair, levels = ordered_gene_pairs)

# Use heatmap data frame to determine the order of the cell lines in all plots
ordered_cell_lines <- unique(full_results$cell_line_label[order(full_results$cancer_type, full_results$cell_line_label)]) 
heatmap_df$cell_line_label <- factor(heatmap_df$cell_line_label, levels = ordered_cell_lines)


# Get total number of hit gene pairs per cell line ------------------------


### Total number of hit gene pairs per cell line ----------------------------

# From the unfiltered data set, determine the number of gene pairs that are hits for each cell line
hit_summary <- get_total_cell_line_hits(full_results)

# Reorder cell line labels to match heatmap
hit_summary$cell_line_label <- factor(hit_summary$cell_line_label, levels = ordered_cell_lines)

# Add total number of hits per cell line to heatmap data frame
heatmap_df <- heatmap_df |> 
  left_join(hit_summary, by = 'cell_line_label')


# Summarise hits per cancer type (cell lines) -----------------------------

# Calculate number of cell lines (per type) in which gene pair is a hit for stacked bar plot
cell_line_summary <- binary_results |>
  filter(sorted_gene_pair %in% hits) |>
  select(sorted_gene_pair, Lung, Melanoma, Pancreatic) |>
  gather(cancer_type, n_cell_lines, -sorted_gene_pair)

# Set order of gene pairs
cell_line_summary$sorted_gene_pair <- factor(cell_line_summary$sorted_gene_pair, levels = ordered_gene_pairs)


# Set up palettes for heatmap ---------------------------------------------

# Colour palette (from -ve to +ve)
heatmap_colors <- colorRampPalette(c('blue', 'white'))(200)

# Add alpha to fade palette 
alpha_cancer_type_palette <- add.alpha(cancer_type_palette, 0.5)
names(alpha_cancer_type_palette) <- c(levels(cell_line_summary$cancer_type))

# Plot landscape version of heatmap and waterfall plot --------------------

# Build lanscape heatmap
landscape_gi_heatmap <- 
  ggplot(heatmap_df, aes(y = cell_line_label, x = sorted_gene_pair, fill = mean_norm_gi)) + 
  geom_tile(data = subset(heatmap_df, is_bassik_hit == 1), # Plot only the hits
            color = "gray30", 
            aes(width = 0.9, height = 0.9), 
            linewidth = 0.4) + 
  scale_fill_gradientn(colors = heatmap_colors, limits = c(-4,0)) + # Set gradient
  coord_cartesian(clip = "off") + # Don't cut off annotation on right-hand side
  theme_classic() + # Set main theme
  labs(x = '', y = '') +
  theme(plot.margin = unit(c(0.001, 2, 1, 1), "lines"), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8, face = "italic"),
        axis.text.y = element_text(size = 8),
        legend.text= element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.position = "top",
        legend.justification = "left",
        legend.box.spacing = unit(0, "pt"))+ # The spacing between the plotting area and the legend box (unit) 
  guides(fill = guide_legend(title = "GI Score"))+ # Rotate x-axis labels (cell lines)
  geom_hline(yintercept = 10.5, color = 'gray30') + # Add vertical line to separate cancer types
  geom_hline(yintercept = 18.5, color = 'gray30') # Add vertical line to separate cancer types

# Bar plot of number of hits per cell line, grouped by cancer type (landscape)
landscape_cell_line_barplot <- 
  ggplot(cell_line_summary, aes(fill = cancer_type, y = n_cell_lines, x = sorted_gene_pair)) + 
  geom_bar(position = "stack", stat = "identity", colour="black") + 
  scale_fill_manual(values = alpha_cancer_type_palette) + # Color by cancer type
  labs(fill = "Cancer Type", x = '', y = 'Number of cell lines') +
  theme_classic() +
  theme(axis.text.x = element_blank(), # Remove gene pairs as labels on y-axis
        axis.ticks = element_blank(), # Remove ticks on both axes
        legend.position="top",
        legend.justification="left",
        legend.box.spacing = unit(0, "pt"),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.text= element_text(size=8),
        legend.title = element_text(size=8)) +
  geom_hline(yintercept = c(5,10,15,20,25), linetype='dashed', col = 'gray', alpha=0.5)+
  theme(plot.margin = unit(c(1, 2, 0.001, 1), "lines"))

# Combine plots into single figure
landscape_combined_figure <- egg::ggarrange(
  landscape_cell_line_barplot, 
  landscape_gi_heatmap, 
  nrow = 2,
  ncol=1,
  widths = c(),
  heights = c(3, 6)) 

# Save landscape plot
#ggsave(file.path(output_plot_dir, 'heatmap_without_total_hits_in_italics_landscape.png'), landscape_combined_figure , dpi = 300, width = 210, height = 200, units = "mm")               


# Plot portrait version of heatmap and waterfall plot ---------------------
## NOTE: this is Figure 2 in the manuscript

portrait_gi_heatmap <- 
  ggplot(heatmap_df, aes(x = cell_line_label, y = sorted_gene_pair, fill = mean_norm_gi)) + 
  geom_tile(data = subset(heatmap_df, is_bassik_hit == 1), # Plot only the hits
            color = "gray30", 
            aes(width = 0.9, height = 0.9), 
            linewidth = 0.4) + 
  scale_fill_gradientn(colors = heatmap_colors, limits = c(-4,0)) + # Set gradient
  coord_cartesian(clip = "off") + # Don't cut off annotation on right-hand side
  theme_classic(base_size = 8) + # Set main theme
  labs(x = '', y = '') +
  theme(plot.margin = unit(c(1, 0, 1, 1), "lines"), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(face = 'italic'),
        legend.text= element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.position="top",
        legend.justification="left",
        legend.box.spacing = unit(0, "pt"))+ # The spacing between the plotting area and the legend box (unit) 
  guides(fill = guide_legend(title = "GI Score"))+ # Rotate x-axis labels (cell lines)
  geom_vline(xintercept = 10.5, color = 'gray30') + # Add vertical line to separate cancer types
  geom_vline(xintercept = 18.5, color = 'gray30') # Add vertical line to separate cancer types

# Bar plot of number of hits per cell line, grouped by cancer type (portrait)
portrait_cell_line_barplot <- 
  ggplot(cell_line_summary, aes(fill = cancer_type, x = n_cell_lines, y = sorted_gene_pair)) + 
  geom_bar(position = "stack", stat = "identity", colour="black") + 
  scale_fill_manual(values = alpha_cancer_type_palette) + # Color by cancer type
  labs(fill = "Cancer Type", y = '', x = 'Number of cell lines') +
  theme_classic(base_size=8) +
  theme(axis.text.y = element_blank(), # Remove gene pairs as labels on y-axis
        axis.ticks = element_blank(), # Remove ticks on both axes
        legend.position="top",
        legend.justification="left",
        legend.box.spacing = unit(0, "pt")) +
  geom_vline(xintercept = c(5,10,15,20,25), linetype='dashed', col = 'gray', alpha=0.5)+
  theme(plot.margin = unit(c(1, 1, 1, 0), "lines"))

# Combine plots into single figure
portrait_combined_figure <- egg::ggarrange(
  portrait_gi_heatmap,
  portrait_cell_line_barplot, 
  nrow = 1,
  ncol=2,
  widths = c(2,1.5),
  heights = c()) 

# Save combined plot
ggsave(file.path(output_plot_dir, 'F2__heatmap_without_total_hits_portrait.png'), portrait_combined_figure , dpi = 300, width = 210, height = 220, units = "mm")               

