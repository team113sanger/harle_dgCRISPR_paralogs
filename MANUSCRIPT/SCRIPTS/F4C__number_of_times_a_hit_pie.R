## Plotting pie chart of number of times a hit split by: 1-5, 6-10, 11-15, 16-20 and 21-25

# Prepare data ------------------------------------------------------------

# Source data preparation script
source(file.path('MANUSCRIPT', 'SCRIPTS', 'prepare_data.R'))


# How many times is each gene pair a hit ----------------------------------

# Summarise hits
number_of_hits <- hits |> select(sorted_gene_pair, Total_hits) |> unique() |> filter(Total_hits >=1)
frequency_of_hits <- tabulate(number_of_hits$Total_hits) 
frequency_of_hits <- as.data.frame(frequency_of_hits) |>
  rownames_to_column(var = "number_of_times_a_hit")

# Group frequencies into bins
grouped_frequencies <- frequency_of_hits |>
  mutate(hit_frequency = case_when(
    (number_of_times_a_hit %in% 1:5)   ~'1-5', 
    (number_of_times_a_hit %in% 6:10)  ~'6-10',
    (number_of_times_a_hit %in% 11:15) ~'11-15',
    (number_of_times_a_hit %in% 16:20) ~'16-20',
    (number_of_times_a_hit %in% 21:25) ~'21-25'))

# Get frequency of each bin
summed_grouped <- grouped_frequencies |> 
  group_by(hit_frequency) |> 
  mutate(frequency = sum(frequency_of_hits)) |>
  select(hit_frequency, frequency) |> 
  unique()

# Generate pie chart
pie_of_hits <- ggplot(summed_grouped , aes(x = "", y = frequency, fill = fct_inorder(hit_frequency))) +
  geom_bar(stat = "identity", width = 0.5, color = "#999999", alpha = 0.6) +
  coord_polar("y", start = 0) +
  theme_void() + 
  scale_fill_manual(values = c("#E69F00","#F0E442","#009E73","#56B4E9", "#0072B2")) +
  labs(title = '', y = '', x = '')+ 
  labs(fill = 'Number of Cell Lines\nGene Pair is called a hit') +
  theme(plot.margin = unit(c(1, 1, 1, 1), "lines"), 
        legend.title = element_text(hjust = 0.5, size = 18),
        legend.text = element_text(size = 18),
        legend.box.spacing = unit(10, "pt"),
        legend.key.height = unit(2, 'cm'),
        legend.key.width = unit(2, 'cm'))

# Save file
ggsave(file.path(output_plot_dir, 'F4C__pie_chart_hit_frequency_altered_colours.png'), pie_of_hits , dpi = 300, width = 12, height = 10)