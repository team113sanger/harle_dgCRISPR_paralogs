## Plotting A549 top hits for imaging GI scores

# Prepare data ------------------------------------------------------------

# Source data preparation script
source(file.path('MANUSCRIPT', 'SCRIPTS', 'prepare_data.R'))


# Filter and label data ---------------------------------------------------

# Filter hits to get only those in A549 
# Label whether the pair is a hit or not in A549
imaging_pairs_in_A549 <- hits |> 
  filter(cell_line_label == 'A549' & Total_hits >= 23) |>
  mutate(Hit = case_when((is_bassik_hit == 1) ~ 'Hit', 
                         (is_bassik_hit == 0) ~ 'Not a Hit'))

p <- ggplot(imaging_pairs_in_A549, 
            aes(x = reorder(sorted_gene_pair, -Total_hits), 
                y = mean_norm_gi, 
                colour = Hit)) + 
  geom_point(alpha = 0.8) + 
  scale_colour_manual(values = c("#0072B2","#D55E00")) +
  labs(x = "", y = 'GI Score', colour = 'Hit Status') +
  theme_bw(base_size = 8) +
  ylim(-2, 0)+
  geom_hline(yintercept = -0.5, colour = 'red', linetype = 'dotted') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "italic"),
        legend.position = "top",
        legend.justification = "left",
        legend.box.spacing = unit(0.1, "pt"),
        legend.key.size = unit(0.5, "lines"),
        plot.margin = unit(c(1, 1, 1, 1), "lines"))

#Save strongest hits ranking plot
ggsave(file.path(output_plot_dir, 'F6A__A549_top_hits_for_imaging_GI_scores.png'), p, dpi = 300, width = 120, height = 100, units = "mm") 
