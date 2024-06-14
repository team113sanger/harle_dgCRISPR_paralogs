## Plotting  the median GI of a gene pair against the number of times the gene pair was a hit, adding a line to show the correlation.

# Prepare data ------------------------------------------------------------

# Source data preparation script
source(file.path('MANUSCRIPT', 'SCRIPTS', 'prepare_data.R'))


# GI vs number of times a hit ---------------------------------------------

# Get hits for GI
hits_for_gi <- get_hits_for_gi(full_results)

# Add ab_line 
# (lm = linear model) - this would be for mean of individual values
gi.lm <- lm(mean_norm_gi ~ total_hits, hits_for_gi)

#Unique median values only for plot
hits_for_gi_unique <- hits_for_gi  |> 
  select(total_hits, median_gi) |> 
  unique()

# Add ab_line 
# (lm = linear model) -  this would be for median of combined values
gi.lm <- lm(median_gi ~ total_hits, hits_for_gi_unique)

# Plot GI v Hits 
median_gi_plot<- ggplot(hits_for_gi_unique, aes(x = total_hits, y = median_gi)) + 
  geom_point(colour = '#0072B2', alpha = 0.5) + 
  labs(y = "Median GI Score of Hits", x = 'Total Number of Cell Lines Gene Pair is a Hit') +
  theme_classic(base_size = 10)  +
  theme(plot.margin = unit(c(1, 1, 1, 1), "lines")) +
  geom_abline(slope = coef(gi.lm)[["total_hits"]], 
              intercept = coef(gi.lm)[["(Intercept)"]])

#Save Median GI (calculated from mean norm GI per hit) to total hits
ggsave(file.path(output_plot_dir, 'F5A__median_GI_of_hits_per_pair_v_total_hits.png'), median_gi_plot, dpi = 300, width = 100, height = 80, units = "mm") 
