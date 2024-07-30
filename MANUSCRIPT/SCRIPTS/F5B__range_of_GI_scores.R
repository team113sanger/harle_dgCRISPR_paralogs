## The range of the GI scores over each cell line (to highlight that some cell lines have really low GI scores and others cell lines have no low GI scores).

# Prepare data ------------------------------------------------------------

# Source data preparation script
source(file.path('MANUSCRIPT', 'SCRIPTS', 'prepare_data.R'))


# Range of GI scores within lines -----------------------------------------

# Filter for only hits
is_hit_only <- hits |> 
  filter(is_bassik_hit == 1)

# Set order of cell lines
ordered_cell_lines <- unique(is_hit_only$cell_line_label[order(is_hit_only$cancer_type, is_hit_only$cell_line_label)]) 
is_hit_only$cell_line_label <- factor(is_hit_only$cell_line_label, levels = ordered_cell_lines)

# Violin plot 
gi_range_of_hits <-
  ggplot(is_hit_only, aes(x = cell_line_label, y = mean_norm_gi, colour = cancer_type)) + 
    geom_violin() + 
    geom_boxplot(width = 0.1) + 
    scale_colour_manual(values = cancer_type_palette) +
    labs(x = "", y = 'GI Score of Hits', colour = 'Cancer Type') +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), #To place text in centre if combining with second figure add hjust=0.5, margin=margin(l= 20)
          legend.position = "top",
          legend.justification = "left",
          legend.box.spacing = unit(0.1, "pt"),
          legend.key.size = unit(0.5, "lines"))  +
    theme(plot.margin = unit(c(1, 1, 1, 1), "lines"))

# Save GI scores of different cell lines 
ggsave(file.path(output_plot_dir, 'F5B__GI_range_of_hits_per_cell_line.png'), gi_range_of_hits, dpi = 300, width = 120, height =100, units = "mm") 
