## Prepare data for manuscript
## This file is sourced by other files in this directory

# Load libraries ----------------------------------------------------------
library(tidyverse)
library(grDevices)
library(ggpubr)
library(egg)
library(readxl)
library(ggVennDiagram)
library(corrr)


# Set colour palette ------------------------------------------------------

# Colour blind friendly R palette with alpha
palette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Set palette, colored by cancer type and add alpha to fade palette
cancer_type_palette <- c('#0072B2', '#E69F00', '#009E73') 


# Set top level directory path --------------------------------------------

# Change this path for scripts to run!
top_dir <- getwd()


# Read in input files -----------------------------------------------------

# Combined gene results
full_results <- read.delim(file.path(top_dir, "combinatorial_crispr_screen_analysis", "DATA", "postprocessing", "combined_gene_level_results.tsv"), header = T, sep = "\t")

# Binary results table
binary_results <- read.delim(file.path(top_dir, "combinatorial_crispr_screen_analysis", "DATA", "postprocessing", "combined_gene_level_results.binary.tsv"), header = T)


# Set output path ---------------------------------------------------------

# Where to write plots to 
output_plot_dir <- file.path(top_dir, 'MANUSCRIPT', 'PLOTS')


# Rename cancer types for figures -----------------------------------------

# Update cancer types for binary results
binary_results <- binary_results |>
  rename('Lung' = 'bassik__Lung_NSCLC',
         'Melanoma' = 'bassik__Melanoma',
         'Pancreatic' = 'bassik__Pancreas')

# Update cancer types in full results
full_results <- full_results |> 
  mutate(across('cancer_type', str_replace, 'Lung NSCLC', 'Lung')) |>
  mutate(across('cancer_type', str_replace, 'Melanoma', 'Melanoma')) |>
  mutate(across('cancer_type', str_replace, 'Pancreas', 'Pancreatic'))


# Summarise hits and GIs --------------------------------------------------

# Add total hits column to the full results (with mean gis)
hits <- full_results |> 
  select(sorted_gene_pair, cell_line_label, cancer_type, mean_norm_gi , is_bassik_hit) |>
  group_by(sorted_gene_pair) |> 
  mutate('Total_hits'= sum(is_bassik_hit)) |> 
  ungroup() 

# Calculate the mean and median GI
gi <- hits |> 
  filter(is_bassik_hit == 1) |> 
  group_by(sorted_gene_pair) |> 
  mutate(mean_GI = mean(mean_norm_gi)) |>
  mutate(median_GI = median(mean_norm_gi))


# Helper functions --------------------------------------------------------
get_total_cell_line_hits <- function(data) {
  # Get only those results which are a hit in bassik analysis
  totals <- data |>
    filter(is_bassik_hit == 1) |> 
    group_by(cell_line_label, cancer_type) |>
    count(name = 'total__hits')
  return(totals)
}

get_hits_for_gi <- function(data) {
  hits_for_gi <- data |> 
    filter(is_bassik_hit == 1) |>
    group_by(sorted_gene_pair) |> 
    mutate('total_hits' = sum(is_bassik_hit)) |> 
    mutate('median_gi' = median(mean_norm_gi)) |> 
    ungroup()
  return(hits_for_gi)
}
