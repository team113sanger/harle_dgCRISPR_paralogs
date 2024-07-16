## Plotting the GI score of each hit in descending order of how many lines the pair was a hit (coloured by cancer type)

# Prepare data ------------------------------------------------------------

# Source data preparation script
source(file.path('MANUSCRIPT', 'SCRIPTS', 'prepare_data.R'))


# Draw GI dot plot with all hits ------------------------------------------

# Build dot plot
all_gi_scores <- ggplot(gi, aes(y = reorder(sorted_gene_pair, Total_hits), x = mean_norm_gi, colour = cancer_type)) + 
  geom_point(alpha = 0.5) + 
  scale_colour_manual(values = cancer_type_palette) +
  scale_y_discrete(position = 'right')+ 
  labs(x = "GI Score of Hits", y ='', colour = 'Cancer Type') +
  theme_bw(base_size = 6) +
  theme(axis.text.y.right = element_text(face = "italic", size = 6), #To place text in centre if combining with second figure add hjust=0.5, margin=margin(l= 20)
        legend.position = "top",
        legend.justification = "left",
        legend.box.spacing = unit(0.1, "pt"),
        legend.key.size = unit(0.5, "lines"))  +
  theme(plot.margin = unit(c(1, 1, 1, 1), "lines"))

# Save GI plot of all hits
ggsave(file.path(output_plot_dir, 'F4A__all_gi_scores.png'), all_gi_scores, dpi = 300, width = 100, height = 280, units = "mm") 