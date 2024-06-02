# Author: Victoria Harle

# Purpose:
# Determine and remove outliers
# Robust Z-transform data with respect to controls

#####################################
#       CHANGE PATH TO DATA         #
#####################################

# Top level directory - can be changed 
plate_directory <- getwd()

#####################################
#   CHANGE BELOW WITH CAUTION       #
#####################################

# Load dependencies
library(tidyverse)
library(ggrepel)
library(scales) 
library(ggpubr)

# Source helper functions
source(file.path(plate_directory, 'SCRIPTS', 'helper.R'))

# Read in raw data list from RDS (see processing_raw_data.R)
data.list <- readRDS(file = file.path(plate_directory, 'DATA', 'raw_processed_plate_list_with_medians.rds'))

# Collate nested list into single data frame
# Remove blank wells
data.df <- convert_nested_list_to_df(data.list)

# Remove wells which have less than 10 analysed fields
data.df <- data.df |> filter(`Number of Analyzed Fields` >= 10)

# Remove wells which have less than 1000 objects
data.df <- data.df |> filter(`Cells Final - Number of Objects` >= 1000)

# Prepare raw control data for PCA (objects and features)
# Use all data otherwise we get straight lines as many columns have little variance or many values close to zero
# This can be tweaked later, but gives a good overview of the data set
raw_pca_data <- data.df |>
  filter(!is.na(Replicate) & Group_Target == 'Control' & Target != 'Parental') |>
  unite(id, Plate, Replicate, Well, Target, Group_Target, sep = '__', remove = F) |>
  select(id, contains('Median per Well'), contains('Object')) |>
  drop_na() |>
  column_to_rownames('id') 

# Build raw control PCA object
# The bits you mainly care about (the principal components are in raw_pca_obj$x)
raw_pca_obj <- prcomp(raw_pca_data)

# Build raw control PCA plots
# To help work out the axis limits run summary(raw_pca_obj$x[,1:3]) and round up the min and max a bit
raw_control_pca_plots <- build_pca_plots(raw_pca_obj, 
                                         pc1_lims = c(-21000, 11000),
                                         pc2_lims = c(-2000, 5000), 
                                         pc3_lims = c(-8000, 2000))

# Save the raw PCA plot of controls for all feature and object columns
# Note: only processing done so far is removing wells with low numbers of objects or analysed fields
ggsave(file.path(plate_directory, 'PLOTS', 'raw_pca.control.all_features_and_objects.png'), raw_control_pca_plots, dpi = 300, width = 12, height = 10)

# Get total number of objects per well for each replicate
# Get the sum of all objects across all plates (full data set)
# Calculate scaling factor as proportion contribution of each replicate to full data set
object_summary <- data.df |> 
  select(Plate:Replicate, `Cells Final - Number of Objects`) |>
  group_by(Replicate) |>
  summarise('total_objects_per_well' = sum(`Cells Final - Number of Objects`), .groups = 'keep') |>
  ungroup() |>
  mutate('total_objects' = sum(total_objects_per_well)) |>
  mutate('scaling_factor' = total_objects_per_well / total_objects)

# Divide each well value by the corresponding scaling factor
# Note: there is no within-plate scaling to controls here
scaled_objects.narrow <- data.df |>
  select(Plate:Replicate, `Non-Proliferative - Number of Objects`:`Cells Final - Number of Objects`) |>
  pivot_longer(cols = `Non-Proliferative - Number of Objects`:`Cells Final - Number of Objects`, values_to = 'raw') |>
  left_join(object_summary |> select(Replicate, scaling_factor), by = c('Replicate'), multiple = "all") |>
  mutate('scaled' = raw / scaling_factor)

# Spread data set so object names are columns again
scaled_objects.wide <- scaled_objects.narrow |>
  pivot_wider(names_from = name, values_from = c(raw, scaled)) 

# Prepare raw control data for PCA (objects)
raw_control_pca_data_objects <- scaled_objects.wide |>
  filter(!is.na(Replicate) & Group_Target == 'Control' & Target != 'Parental') |>
  unite(id, Plate, Replicate, Well, Target, Group_Target, sep = '__', remove = F) |>
  select(id, matches("raw.*Object"), -contains('Cells')) |>
  drop_na() |>
  column_to_rownames('id') 

# Build raw control PCA object
# The bits you mainly care about (the principal components are in raw_control_objects_pca_obj$x)
raw_control_objects_pca_obj <- prcomp(raw_control_pca_data_objects)

# Build raw control PCA plots
# To help work out the axis limits run summary(raw_control_objects_pca_obj$x[,1:3]) and round up the min and max a bit
raw_control_object_pca_plots <- build_pca_plots(raw_control_objects_pca_obj, 
                                                   pc1_lims = c(-5000, 11000),
                                                   pc2_lims = c(-2000, 6000), 
                                                   pc3_lims = c(-3000, 6000))

# Save the raw PCA plot of controls for feature columns only
# Note: only processing done so far is removing wells with low numbers of objects or analysed fields
ggsave(file.path(plate_directory, 'PLOTS', 'raw_pca.control.objects_only.png'), raw_control_object_pca_plots, dpi = 300, width = 12, height = 10)

# Prepare scaled control data for PCA (objects)
scaled_control_pca_data_objects <- scaled_objects.wide |>
  filter(!is.na(Replicate) & Group_Target == 'Control' & Target != 'Parental') |>
  unite(id, Plate, Replicate, Well, Target, Group_Target, sep = '__', remove = F) |>
  select(id, matches("scaled.*Object"), -contains('Cells')) |>
  drop_na() |>
  column_to_rownames('id') 

# Build scaled control PCA object
# The bits you mainly care about (the principal components are in scaled_control_objects_pca_obj$x)
scaled_control_objects_pca_obj <- prcomp(scaled_control_pca_data_objects)

# Build scaled control PCA plots
# To help work out the axis limits run summary(scaled_control_objects_pca_obj$x[,1:3]) and round up the min and max a bit
scaled_control_object_pca_plots <- build_pca_plots(scaled_control_objects_pca_obj, 
                                                   pc1_lims = c(-18000, 30000),
                                                   pc2_lims = c(-7000, 25000), 
                                                   pc3_lims = c(-9000, 21000))

# Save the scaled PCA plot of controls for object columns only
ggsave(file.path(plate_directory, 'PLOTS', 'scaled_pca.control.objects_only.png'), scaled_control_object_pca_plots, dpi = 300, width = 12, height = 10)

# Prepare data to compare scaled vs raw object counts as boxplots
# Remove unwanted columns
# Filter to get only object counts and exclude unclassified
scaled_vs_raw_obj.narrow <- scaled_objects.narrow |>
  select(-scaling_factor) |>
  mutate(name = gsub(' - Number of Objects', '', name)) |>
  filter(name != 'Cells Final' & Target != 'Parental' & name != 'Unclassified') |>
  pivot_longer(cols = c(raw, scaled), names_to = 'dataset')

# Loop over the classifications (e.g. Apoptotic, Proliferative)
# Build boxplot where Target is on the x-axis (e.g. ASF1A) and count is on the y-axis
# Facet rows are scaled or raw and columns are the Group_Target (e.g. ASF1A_ASF1B)
for (n in unique(scaled_vs_raw_obj.narrow$name)) {
  # Build barplot for cell classification
  scaled_vs_norm_obj_plot <- 
    ggplot(scaled_vs_raw_obj.narrow |> filter(name == n), 
           aes(x = Target, y = value, fill = Group_Target)) + 
      geom_boxplot() +
      facet_grid(dataset ~ Group_Target, space = 'free_x', scales = 'free') + 
      scale_y_continuous(breaks = pretty_breaks(6)) +
      theme_pubr() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  # Save barplot for cell classification
  ggsave(file.path(plate_directory, 'PLOTS', paste0('scaled_vs_raw_barplot.', n, '.objects_only.png')), 
         scaled_vs_norm_obj_plot, dpi = 300, width = 12, height = 10)
}

# Calculate median per feature for all controls except Parental (i.e. Control_1, Control_2 and Control_1|Control_2)
# Selects the controls
# Selects the annotation columns (Plate:Replicate) and the feature columns (i.e. Median per Well)
# Moves features from being column names to rows
# Groups by Replicate and Feature
# Calculates the median of the controls per Replicate for each Feature
# Calculates the absolute deviation (i.e. takes the positive value when subtracting median of controls from each raw value per well)
# Calculates the median of the control absolute deviations per Replicate per Feature (MAD)
# Gets only the columns we want and makes sure rows are unique
control_median_per_feature <- data.df |>
  filter(Target != 'Parental' & Group_Target == 'Control') |>
  select(Plate:Replicate, contains('Median per Well')) |>
  pivot_longer(cols = `Cells Final - Cell 488 Axial Length Ratio - Median per Well`:`Cells Final - Total Spot Area - Median per Well`, names_to = 'Feature', values_to = 'raw') |>
  group_by(Replicate, Feature) |>
  mutate('control_median' = median(raw),
         'control_abs_deviation' = abs(raw - control_median),
         'control_mad' = median(control_abs_deviation)) |>
  select(Plate, Replicate, Feature, control_median, control_mad) |>
  unique()

# Selects the annotation columns (Plate:Replicate) and the feature columns (i.e. Median per Well)
# Moves features from being column names to rows
# Adds the control statistics (e.g. median and MAD per Feature per Replicate)
# Robust Z-score for features = (1.4826 * (x - median(control))) / mad(control)
# Where control MAD is 0, set the Robust_Z to 0 (instead of -Inf or Inf as dividing by 0)
robust_z <- data.df |>
  select(Plate:Replicate, contains('Median per Well')) |>
  pivot_longer(cols = `Cells Final - Cell 488 Axial Length Ratio - Median per Well`:`Cells Final - Total Spot Area - Median per Well`, names_to = 'Feature', values_to = 'raw') |>
  left_join(control_median_per_feature, by = c('Plate', 'Replicate', 'Feature'), multiple = "all") |>
  mutate('Robust_Z' = (1.4826 * (raw - control_median)) / control_mad) |>
  mutate('Robust_Z' = ifelse(Robust_Z == Inf | Robust_Z == -Inf | is.na(Robust_Z), 0 , Robust_Z))

# Remove unwanted columns
# Make features column names rather than rows
robust_z.wide <- robust_z |>
  select(-control_median, -control_mad) |>
  pivot_wider(names_from = Feature, values_from = c(raw, Robust_Z)) 
  
# Prepare raw control data for PCA (Feature Median per Well)
raw_control_pca_data_features <- data.df |>
  filter(!is.na(Replicate) & Group_Target == 'Control' & Target != 'Parental') |>
  unite(id, Plate, Replicate, Well, Target, Group_Target, sep = '__', remove = F) |>
  select(id, contains("Median per Well")) |>
  drop_na() |>
  column_to_rownames('id') 

# Build raw control PCA object (features only)
# The bits you mainly care about (the principal components are in raw_control_features_pca_obj$x)
raw_control_features_pca_obj <- prcomp(raw_control_pca_data_features)

# Build raw control PCA plots (features only)
# To help work out the axis limits run summary(raw_control_features_pca_obj$x[,1:3]) and round up the min and max a bit
# Note: we see some straight lines in scatter plot as many features have either low variance or are all close to zero
raw_control_features_pca_plots <- build_pca_plots(raw_control_features_pca_obj, 
                                                  pc1_lims = c(-500, 400),
                                                  pc2_lims = c(-25, 5), 
                                                  pc3_lims = c(-5, 5))

# Save the raw PCA plot of controls for feature columns only
ggsave(file.path(plate_directory, 'PLOTS', 'raw_pca.control.features_only.png'), raw_control_features_pca_plots, dpi = 300, width = 12, height = 10)

# Prepare scaled control data for PCA (Robust Z of Feature Median per Well)
scaled_control_pca_data_features <- robust_z.wide |>
  filter(!is.na(Replicate) & Group_Target == 'Control' & Target != 'Parental') |>
  unite(id, Plate, Replicate, Well, Target, Group_Target, sep = '__', remove = F) |>
  select(id, matches("Robust_Z")) |>
  drop_na() |>
  column_to_rownames('id') 

# Build scaled control PCA object (features only)
# The bits you mainly care about (the principal components are in scaled_control_features_pca_obj$x)
scaled_control_features_pca_obj <- prcomp(scaled_control_pca_data_features)

# Build scaled control PCA plots (features only)
# To help work out the axis limits run summary(scaled_control_features_pca_obj$x[,1:3]) and round up the min and max a bit
scaled_control_features_pca_plots <- build_pca_plots(scaled_control_features_pca_obj, 
                                                   pc1_lims = c(-50, 200),
                                                   pc2_lims = c(-200, 100), 
                                                   pc3_lims = c(-50, 50),
                                                   label_outliers_pc1 = TRUE, pc1_label_lims = c(30, -40),
                                                   label_outliers_pc2 = TRUE, pc2_label_lims = c(50, -60),
                                                   label_outliers_pc3 = TRUE, pc3_label_lims = c(50, -50))

# Save the scaled PCA plot of controls for feature columns only
ggsave(file.path(plate_directory, 'PLOTS', 'scaled_pca.control.features_only_all_excluded.png'), scaled_control_features_pca_plots, dpi = 300, width = 12, height = 10)
