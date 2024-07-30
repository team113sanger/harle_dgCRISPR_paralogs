# 7th_May_24 - plotting the mouse phenotyping data for Slc25a28. 
# The data was downloaded from https://www.mousephenotype.org on 7th May 2024 (Release 21.0)

# Load libraries
library(tidyverse)
library(curl)

# Set output plot directory
top_dir <- getwd()
output_plot_dir <- file.path(top_dir, 'MANUSCRIPT', 'PLOTS')

# Download the dataset
url <- "https://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/release-21.0/results/statistical-results-ALL.csv.gz"
temp <- tempfile()
curl_download(url, temp)
Release_21 <- read.csv(gzfile(temp), stringsAsFactors = FALSE)
unlink(temp)

# Extract sata for gene(s) of interest
# Calculate log10(p_value)
Slc25a28 <- Release_21 |> 
  filter(marker_symbol == "Slc25a28") |>
  mutate(p_value_log10 = log10(p_value), .after = p_value)

# Plot the data
p <- Slc25a28 |>
  ggplot(aes(x = fct_inorder(procedure_name), 
             y = p_value_log10, 
             colour = parameter_name, 
             label = parameter_name)) +
  geom_point(show.legend = FALSE) +
  scale_y_reverse() +
  theme_classic() + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 70, vjust = 1, hjust = 1))

# Save plot
ggsave(file.path(output_plot_dir, 'EF6D__Slc25a28_IMPC.png'), p, dpi = 300, width = 300, height = 200, units = "mm") 
