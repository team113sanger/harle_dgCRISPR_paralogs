## Ranking of top 50% of hits to show if the most common hits were also the strongest hits within each line

# Prepare data ------------------------------------------------------------

# Source data preparation script
source(file.path('MANUSCRIPT', 'SCRIPTS', 'prepare_data.R'))


# Calculate ranking and determine if top 24 hits are also top ranked  -----

# Get hits for GI
hits_for_gi <- get_hits_for_gi(full_results)

# Add ranking
top_ranked <- hits_for_gi |> 
  filter(is_bassik_hit == 1) |> 
  group_by(cell_line_label) |> 
  mutate(Rank = rank(mean_norm_gi, ties.method = "first")) |> 
  ungroup() 

# Filter for only the strongest hits (more than 50% of cell lines)
strongest_hits_ranking <- top_ranked |> 
  filter(total_hits >= 14)

# Plot strongest hits
strongest_hits_ranking_plot <- ggplot(strongest_hits_ranking, aes(x = reorder(sorted_gene_pair, -total_hits), y = Rank, colour = cancer_type)) + 
  geom_point() + 
  scale_colour_manual(values = cancer_type_palette) +
  labs(x = "", y = 'Rank of Gene Pair Hit Per Cell Line', colour = 'Cancer Type') +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "italic"), #To place text in centre if combining with second figure add hjust=0.5, margin=margin(l= 20)
        legend.position = "top",
        legend.justification = "left",
        legend.box.spacing = unit(0.1, "pt"),
        legend.key.size = unit(0.5, "lines"))  +
  theme(plot.margin = unit(c(1, 1, 1, 1), "lines"))

# Save strongest hits ranking plot
ggsave(file.path(output_plot_dir, 'F5C__cell_line_ranking_of_hits_in_more_than_50pct_of_lines.png'), strongest_hits_ranking_plot, dpi = 300, width = 120, height =100, units = "mm") 
