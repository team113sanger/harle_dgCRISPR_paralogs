# Script description


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(readxl)

# Set paths ---------------------------------------------------------------

# Set top level directory
top_dir <- getwd()

# Source helper script
source(file.path(top_dir, 'SCRIPTS', 'helper.R'))

# Set path to input data
plate_directory <- file.path(top_dir, 'DATA')

# Set path to output data
output_dir <- file.path(top_dir, 'RESULTS')


# Prepare reusable plate labels -------------------------------------------

all_plate_labels <- prepare_plate_labels(plate_directory)

# Save all plate labels as a TSV
write.table(all_plate_labels[['all_plate_labels']],
            file.path(output_dir, 'all_plate_labels.tsv'),
            sep = "\t", row.names = F, quote = F)

# Save plate label list as RDS
saveRDS(all_plate_labels[['plate_label_list']], file.path(output_dir, 'plate_label_list.rds'))


# Raw plate data ----------------------------------------------------------

# Get list of raw plate data files
raw_plate_data_files <- list.files(file.path(plate_directory, 'raw_data'),
                                   pattern = "Objects_Population - Cells Final.txt", 
                                   full.names = T, ignore.case = T, recursive = T)

# Read raw plate data into list (raw_plate_list[[plate name]][[replicate]])
raw_plate_list <- read_raw_plates(raw_plate_data_files, all_plate_labels[['plate_label_list']])


# Calculate summary statistics across all groups --------------------------

# Calculate cell class (e.g. Proliferative, Enlarged..) statistics per well for each plate
cell_class_stats <- calculate_cell_class_stats_by_plate(raw_plate_list)

# Calculate number of analysed fields per well for each plate (i.e. number of unique Field values in each well)
num_analysed_fields <- calculate_num_analysed_fields_by_plate(raw_plate_list)

# Calculate mean, median and standard deviation per well
# Note: despite parallelisation this step can take several minutes to complete
characteristic_stats <- calculate_characteristic_stats_by_plate(raw_plate_list)

# Collate summary statistics into a single data frame per plate
raw_processed_plate_list <- 
  prepare_summary_stat_df(all_plate_labels[['plate_label_list']], cell_class_stats, num_analysed_fields, characteristic_stats)


# Save data objects -------------------------------------------------------

# Write raw plate data list to RDS
print('Writing raw plate list to RDS.')
saveRDS(raw_plate_list, file.path(output_dir, 'raw_plate_list.rds'))

# Write raw cell classification stats to RDS
print('Writing raw cell classification statistics list to RDS.')
saveRDS(cell_class_stats, file.path(output_dir, 'cell_class_statistics_list.rds')) 

# Write raw processed plate data list including medians to RDS
print('Writing raw processed plate list with medians to RDS.')
saveRDS(raw_processed_plate_list, file.path(output_dir, 'raw_processed_plate_list_with_medians.rds'))P: