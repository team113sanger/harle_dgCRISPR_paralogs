## Pie chart of all genes pairs split as: never a hit, context dependent i.e. hit in less than 50% of lines, Strong hit i.e. a hit in more than 50% of lines 

# Prepare data ------------------------------------------------------------

# Source data preparation script
source(file.path('MANUSCRIPT', 'SCRIPTS', 'prepare_data.R'))


# Percentage of library that was a hit ever -------------------------------

# Pie chart including 0
all_pairs <- binary_results |> 
  mutate(hit_frequency = case_when(
    (`bassik__total` %in% 0)  ~'Never a Hit', 
    (`bassik__total` %in% 1:14)  ~'Context Dependent Hit\n < 50% of lines',
    (`bassik__total` %in% 15:26)  ~'Strong Hit > 50% of lines')) |> 
  select(sorted_gene_pair, hit_frequency)

# Get frequency
all_pairs_freq <- all_pairs |> 
  group_by(hit_frequency) |> 
  summarise(n=n())

# Reorder classifications
order_of_classifications <- c('Never a Hit','Context Dependent Hit\n < 50% of lines','Strong Hit > 50% of lines')
ordered_freq <- all_pairs_freq |> 
  arrange(factor(hit_frequency, levels = order_of_classifications))

# Pie chart of all gene pairs
pie_of_all <- ggplot(ordered_freq , aes(x = "", y = n , fill= reorder(hit_frequency, -n))) +
  geom_bar(stat = "identity", width = 0.5, color = "#999999", alpha = 0.6) +
  coord_polar("y", start = 0) +
  theme_void() + 
  scale_fill_manual(values = c("#E69F00","#F0E442", "#009E73")) +
  labs(title = '', y = '', x = '')+ 
  labs(fill = 'Number of Gene Pairs') +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 10)+
  theme(plot.margin = unit(c(1, 1, 1, 1), "lines"), 
        legend.title = element_text(hjust = 0.5, size = 18),
        legend.text = element_text(size = 18),
        legend.box.spacing = unit(5, "pt"),
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(1, 'cm'))

# Save plot
ggsave(file.path(output_plot_dir, 'F4B__pie_chart_frequency_of_all_pairs.png'), pie_of_all, dpi = 300, width = 12, height = 10)
