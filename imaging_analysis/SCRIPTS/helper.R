# Purpose:
# Functions for ingesting and processing imaging data

# Function to read in raw plate data --------------------------------------

read_raw_plates <- function(files = NULL, plate_label_list = NULL) {
  
  # Check we have files in the vector
  if (0 == length(files) | is.null(files)) { 
    stop(print('Raw plate files cannot be read: no raw plate files were in the list.')) 
  }
  
  # Check we have labels in list
  if (is.null(plate_label_list) | 0 == length(plate_label_list)) { 
    stop(print('Raw plate files cannot be read: no plate labels were in the list.')) 
  }
   
  # Create empty list for raw plate data sets
  list_of_plates <- list()
  
  # Loop over raw plate data files
  for (i in files) {
    # Logging
    print(paste("Reading in raw plate file:", i))
    
    # Get plate name from file name
    tmp_plate <- paste0('Plate_', gsub(".*Plate (.*)__.*", "\\1", i))
    
    # Get replicate from file name
    tmp_rep <- ifelse(!grepl('Set', i), 'N1', paste0('N', gsub(".*Screen Set (.*) Plate.*", "\\1", i)))
    
    # Read in plate data to data frame
    tmp_data <- suppressMessages(read_delim(i, delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 9))
    
    # Remove columns which are all NaN, NA or 0
    columns_to_remove <- names(which(colSums(is.na(tmp_data)) == nrow(tmp_data)))
    print(paste("Columns to remove (contain all NA or NaN):", paste(columns_to_remove, collapse = ', ')))
    tmp_data <- tmp_data |>
      select(-columns_to_remove)
    
    # Add plate label
    tmp_data <- tmp_data |> 
      mutate('Plate' = tmp_plate, 'Replicate' = tmp_rep, .before = 'Row')
    
    # Merge Row and Column to get Well
    tmp_data <- tmp_data |> unite(Well, Row, Column, sep = ',', remove = F) 
    
    # Add plate labels
    tmp_data_with_labels <- plate_label_list[[tmp_plate]] |> 
      full_join(tmp_data, by = c('Plate', 'Well', 'Row', 'Column'), multiple = "all")
    
    # Remove wells which are blank (Target == BLANK)
    tmp_data_with_labels <- tmp_data_with_labels |>
      filter(Target != 'BLANK')
    
    # Add to list
    list_of_plates[[tmp_plate]][[tmp_rep]] <- tmp_data_with_labels
    
    # Clean up
    rm(list = c(ls(pattern = 'tmp')))
  }
  
  # Clean up 
  rm(i)
  
  # Check list has populated
  if (is.null(list_of_plates) | 0 == length(list_of_plates)) { 
    stop(print('Raw plate files cannot be read: list is empty or null.')) 
  }
  
  # Return populated plate list
  return(list_of_plates)
}


# Function to read in old processed plate data ----------------------------

read_old_processed_plates <- function(files = NULL, plate_label_list = NULL) {
  
  # Check we have files in the vector
  if (0 == length(files) | is.null(files)) { 
    stop(print('Old processed plate files cannot be read: no raw plate files were in the list.')) 
  }
  
  # Check we have labels in list
  if (is.null(plate_label_list) | 0 == length(plate_label_list)) { 
    stop(print('Old processed plate files cannot be read: no plate labels were in the list.')) 
  }
  
  # Create empty list for raw plate data sets
  list_of_plates <- list()
  
  # Loop over plate data files
  for (i in files) {
    
    # Logging
    print(paste("Reading in old processed plate file:", i))
    
    # Get plate name from filename
    tmp_plate <- gsub("_N.*.txt", '', basename(i))
    
    # Get replicate from filename
    tmp_rep <- gsub("Plate.*_(.*).txt", '\\1', basename(i))
    
    # Read in plate data to data frame
    tmp_data <- suppressMessages(read_delim(i, delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 8))
    
    # Remove columns which are all NaN, NA or 0
    columns_to_remove <- names(which(colSums(is.na(tmp_data)) == nrow(tmp_data)))
    print(paste("Columns to remove (contain all NA or NaN):", paste(columns_to_remove, collapse = ', ')))
    tmp_data <- tmp_data |>
      select(-c(columns_to_remove))
    
    # Add plate label
    tmp_data <- tmp_data |> 
      mutate('Plate' = tmp_plate, 'Replicate' = tmp_rep, .before = 'Row')
    
    # Merge Row and Column to get Well
    tmp_data <- tmp_data |> unite(Well, Row, Column, sep = ',', remove = F) 
    
    # Remove unwanted columns
    tmp_data <- tmp_data |> 
      select(-Timepoint, -`Time [s]`)
    
    # Add plate labels
    tmp_data_with_labels <- plate_label_list[[tmp_plate]] |> 
      full_join(tmp_data, by = c('Plate', 'Well', 'Row', 'Column'))
    
    # Remove wells which are blank (Target == BLANK)
    tmp_data_with_labels <- tmp_data_with_labels |>
      filter(Target != 'BLANK')
    
    # Add to list
    list_of_plates[[tmp_plate]][[tmp_rep]] <- tmp_data_with_labels
    
    # Clean up
    rm(list = c(ls(pattern = 'tmp')))
  }
  
  # Clean up 
  rm(i)
  
  # Check list has populated
  if (is.null(list_of_plates) | 0 == length(list_of_plates)) { 
    stop(print('Old processed plate files cannot be read: list is empty or null.')) 
  }
  
  # Return populated plate list
  return(list_of_plates)
}


# Calculate cell class statistics per plate -------------------------------

calculate_cell_class_stats_by_plate <- function(data = NULL) {
  # Check list is populated
  if (is.null(data) | 0 == length(data)) { 
    stop(print('Cannot calculate cell class statistics: input list is empty or null.')) 
  }
  
  # Set up empty list
  list_of_cell_class_stats <- list() 
  
  # Loop over plates
  for (tmp_plate in names(data)) {
   # if(tmp_plate != 'Plate_4') {next;}
    # Loop over replicates
    for (tmp_rep in names(data[[tmp_plate]])) {
      # Logging
      print(paste("Calculating cell class stats:", tmp_plate, tmp_rep))
      
      # Pull data for a single plate
      tmp_data <- data[[tmp_plate]][[tmp_rep]]
      
      # Get the total number of objects
      tmp_total_objects <- tmp_data |> 
        group_by(Plate, Replicate, Well, Column, Row, `Cells Final - Class`) |> 
        count(name = 'Number of Objects')
      
      # Add the total number of objects per class
      tmp_total_objects_per_class <- tmp_total_objects |>
        filter(`Cells Final - Class` != 'NA') |>
        mutate('Cells Final - Class' = case_when(`Cells Final - Class` == 'A' ~ 'Non-Proliferative',
                                                 `Cells Final - Class` == 'B' ~ 'Proliferative',
                                                 `Cells Final - Class` == 'C' ~ 'Apoptotic',
                                                 `Cells Final - Class` == 'D' ~ 'Enlarged',
                                                 `Cells Final - Class` == 'UnClassified' ~ 'Unclassified'))
        
      # Spread the table so cell classes become column names
      tmp_total_objects_per_class <- tmp_total_objects_per_class |>
        pivot_wider(names_from = `Cells Final - Class`, values_from = `Number of Objects`) 
      
      # Calculate percentage of cells by class
      tmp_total_objects_per_class <- tmp_total_objects_per_class |>
        rowwise() |>
        mutate('Cells Final - Number of Objects' = sum(c_across(c(`Non-Proliferative`:`Unclassified`)), na.rm = T)) |>
        mutate('% Non-Proliferative' = ((`Non-Proliferative` / `Cells Final - Number of Objects`) * 100),
               '% Proliferative' = ((`Proliferative` / `Cells Final - Number of Objects`) * 100),
               '% Apoptotic' = ((`Apoptotic` / `Cells Final - Number of Objects`) * 100),
               '% Enlarged' = ((`Enlarged` / `Cells Final - Number of Objects`) * 100),
               '% Unclassified' = ((`Unclassified` / `Cells Final - Number of Objects`) * 100))
      
      # Rename columns to match existing summarised data columns
      tmp_total_objects_per_class <- tmp_total_objects_per_class |>
        rename('Non-Proliferative - Number of Objects' = 'Non-Proliferative',
               'Proliferative - Number of Objects' = 'Proliferative',
               'Apoptotic - Number of Objects' = 'Apoptotic',
               'Enlarged - Number of Objects' = 'Enlarged',
               'Unclassified - Number of Objects' = 'Unclassified')
      
      # Add to the list 
      list_of_cell_class_stats[[tmp_plate]][[tmp_rep]] <- tmp_total_objects_per_class
    
      # Clean up
      rm(list = c('tmp_data', 'tmp_total_objects_per_class', 'tmp_total_objects'))
    }
  }
  
  # Check list is populated
  if (is.null(list_of_cell_class_stats) | 0 == length(list_of_cell_class_stats)) { 
    stop(print('Cannot calculate cell class statistics: output list is empty or null.')) 
  }
  
  # Clean up
  rm(list = c('tmp_plate', 'tmp_rep'))
  
  # Return cell class statistics for all plates
  return(list_of_cell_class_stats)
}


# Calculate number of analysed fields per plate ---------------------------

calculate_num_analysed_fields_by_plate <- function(data = NULL) {
  # Check list is populated
  if (is.null(data) | 0 == length(data)) { 
    stop(print('Cannot calculate number of analysed fields: input list is empty or null.')) 
  }
  
  # Set up empty list
  list_of_num_fields_stats <- list() 
  
  # Loop over plates
  for (tmp_plate in names(data)) {
    # Loop over replicates
    for (tmp_rep in names(data[[tmp_plate]])) {
      # Logging
      print(paste("Calculating number of analysed fields:", tmp_plate, tmp_rep))
      
      # Pull data for a single plate
      tmp_data <- data[[tmp_plate]][[tmp_rep]]
      
      tmp_num_analysed_fields <- tmp_data |>
        group_by(Plate, Replicate, Well, Column, Row) |> 
        summarise('Number of Analyzed Fields' = n_distinct(`Field`), .groups = 'keep')
      
      # Add to the list 
      list_of_num_fields_stats[[tmp_plate]][[tmp_rep]] <- tmp_num_analysed_fields
      
      # Clean up
      rm(list = c('tmp_data', 'tmp_num_analysed_fields'))
    }
  }
  
  # Check list is populated
  if (is.null(list_of_num_fields_stats) | 0 == length(list_of_num_fields_stats)) { 
    stop(print('Cannot calculate number of analysed fields: output list is empty or null.')) 
  }
  
  # Clean up
  rm(list = c('tmp_plate', 'tmp_rep'))
  
  # Return cell class statistics for all plates
  return(list_of_num_fields_stats)
}


# Calculate mean, median and standard deviation per plate -----------------

calculate_characteristic_stats_by_plate <- function(data = NULL) {
  # Check list is populated
  if (is.null(data) | 0 == length(data)) { 
    stop(print('Cannot calculate characteristic stats: input list is empty or null.')) 
  }
  
  # Set up empty list
  list_of_characteristic_stats <- list() 
  
  # Loop over plates
  for (tmp_plate in names(data)) {
    # Loop over replicates
    for (tmp_rep in names(data[[tmp_plate]])) {
      # Logging
      print(paste("Calculating characteristic stats:", tmp_plate, tmp_rep))
      
      # Pull data for a single plate
      tmp_data <- data[[tmp_plate]][[tmp_rep]]
      
      # Select only the columns required and collapse the columns
      tmp_characteristic_stats <- tmp_data |> 
        select(Plate, Replicate, Well, Column, Row, 
               `Cells Final - Cell Area [µm²]`:`Cells Final - Number of Spots per Area of Cytoplasm`, 
               `Cells Final - Regression A-B`, `Cells Final - Non-Proliferative`:`Cells Final - Enlarged`) |>
        pivot_longer(cols = `Cells Final - Cell Area [µm²]`:`Cells Final - Enlarged`)
      
      print(paste('Num rows before removing empty data:', nrow(tmp_characteristic_stats)))
      
      # Remove values which have NA
      tmp_characteristic_stats <-  tmp_characteristic_stats |>
        filter(!is.na(value) & value != 'NaN' & value != '' & !is.nan(value))
      
      print(paste('Num rows after removing empty data:', nrow(tmp_characteristic_stats)))
      
      # Calculate characteristic stats (e.g. mean, median and standard deviation) per well
      tmp_characteristic_stats <- tmp_characteristic_stats |>
        group_by(Plate, Replicate, Well, Column, Row, name) |>
        summarise('Mean per Well' = mean(value, na.rm = TRUE),
                  'Median per Well' = median(value, na.rm = TRUE),
                  'StdDev per Well' = sd(value, na.rm = TRUE), 
                  .groups = 'keep')
      
      # Spread the column stats back out wide again
      tmp_characteristic_stats <- tmp_characteristic_stats |>
        pivot_wider(names_from = name, 
                    names_glue = "{name} - {.value}",
                    values_from = c(`Mean per Well`, `Median per Well`, `StdDev per Well`))
      
      # Add to the list 
      list_of_characteristic_stats[[tmp_plate]][[tmp_rep]] <- tmp_characteristic_stats
      
      # Clean up
      rm(list = c('tmp_data', 'tmp_characteristic_stats'))
    }
  }
  
  # Clean up
  rm(list = c('tmp_plate', 'tmp_rep'))
  
  # Return cell class statistics for all plates
  return(list_of_characteristic_stats)
}

# Prepare summary statistic data frame ------------------------------------

prepare_summary_stat_df <- function(plate_label_list = NULL, cell_class = NULL, analysed_fields = NULL, characteristic_stats = NULL) {
  
  # Set up empty list
  list_of_processed_plate_stats <- list() 
  
  # Loop over plates
  for (tmp_plate in names(cell_class)) {
    # Loop over replicates
    for (tmp_rep in names(cell_class[[tmp_plate]])) {
      # Logging
      print(paste("Collating plate stats:", tmp_plate, tmp_rep))
      
      # Pull data for a single plate
      tmp_plate_labels         <- plate_label_list[[tmp_plate]]
      tmp_cell_class           <- cell_class[[tmp_plate]][[tmp_rep]]
      tmp_analysed_fields      <- analysed_fields[[tmp_plate]][[tmp_rep]]
      tmp_characteristic_stats <- characteristic_stats[[tmp_plate]][[tmp_rep]]
      
      # Collate stats for a single plate
      tmp_stats <- tmp_plate_labels |>
        right_join(tmp_characteristic_stats, by = c('Plate', 'Well', 'Row', 'Column')) |>
        full_join(tmp_cell_class, by = c('Plate', 'Replicate', 'Well', 'Column', 'Row')) |>
        full_join(tmp_analysed_fields, by = c('Plate', 'Replicate', 'Well', 'Column', 'Row'))
      
      # Add to the list 
      list_of_processed_plate_stats[[tmp_plate]][[tmp_rep]] <- tmp_stats
      
      # Clean up
      rm(list = c('tmp_stats', 'tmp_plate_labels', 'tmp_cell_class', 'tmp_analysed_fields', 'tmp_characteristic_stats'))
    }
  }
  
  # Clean up
  rm(list = c('tmp_plate', 'tmp_rep'))
  
  # Return cell class statistics for all plates
  return(list_of_processed_plate_stats)
}


# Compare old and new processed data sets ---------------------------------

compare_old_and_new_data <- function(old_data = NULL, new_data = NULL) {
  
  # Set up empty data frame for mismatching data
  tmp_diffs <- data.frame()
  
  # Loop over plates
  for (tmp_plate in names(new_data)) {
    # Loop over replicates
    for (tmp_rep in names(new_data[[tmp_plate]])) {

      # Logging
      print(paste("Comparing old and new processed plate data:", tmp_plate, tmp_rep))
      
      # Pull data for a single plate
      tmp_old_data <- old_data[[tmp_plate]][[tmp_rep]]
      tmp_new_data <- new_data[[tmp_plate]][[tmp_rep]]
      
      # Check all columns exist in both data sets
      # Remove all white space as not the same between column names in same data set
      tmp_diff_colnames <- setdiff(gsub(' ', '', colnames(tmp_old_data)), gsub(' ', '', colnames(tmp_new_data)))
      if (length(tmp_diff_colnames) > 0){
        print(paste("Old and new column names differ:", paste(tmp_diff_colnames, collapse = ', ')))
      } else {
        print("Old and new column names are the same.")
      }
      
      # Pivot old data to long format (Features: column names)
      # Remove empty wells
      tmp_old_data <- tmp_old_data |>
        filter(Target != 'BLANK') |>
        pivot_longer(cols = !c(`Plate`:`Replicate`), names_to = 'Features', values_to = 'old') |>
          mutate('Features' = gsub(' ', '', Features))
      
      # Pivot new data to long format (Features: column names)
      tmp_new_data <- tmp_new_data |>
        filter(Target != 'BLANK') |>
        pivot_longer(cols = !c(`Plate`:`Replicate`), names_to = 'Features', values_to = 'new') |>
        mutate('Features' = gsub(' ', '', Features))
      
      # Join old and new data where column names are shared
      tmp_all_data <- tmp_old_data |> 
        full_join(tmp_new_data, by = c('Plate', 'Position', 'Well', 'Row', 'Column', 'Target', 'Group_Target', 'siRNA_target_A', 'siRNA_target_B', 'Replicate', 'Features'))
      
      # Compare old and new values for number of objects
      # is_equal = 1 (same), is_equal = 0 (differ)
      tmp_num_obj <- tmp_all_data |>
        filter(grepl('NumberofObjects', Features)) |>
        mutate('is_equal' = ifelse(old == new, 1, 0))
      tmp_num_obj_diff <- tmp_num_obj |> filter(is_equal == 0) |> mutate('Reason' = 'number of object values do not match')
      
      # Print out if number of objects differ
      if (nrow(tmp_num_obj_diff) > 0) {
        print(paste("Total number of objects with matching values differing (mismatch):", nrow(tmp_num_obj_diff)))
        print(paste("Names of number of object columns with matching values differing (mismatch):", paste(tmp_num_obj$Features, collapse = ', ')))
        if (nrow(tmp_diffs) == 0) {
          tmp_diffs <- tmp_num_obj_diff
        } else {
          tmp_diffs <- rbind(tmp_diffs, tmp_num_obj_diff)
        }
      } else {
        print("Number of object column values are the same.")
      }
      
      # Compare columns containing mean values to 6 decimal places
      tmp_mean_vals <- tmp_all_data |>
        filter(grepl('Mean', Features) & !grepl('Median', Features)) |>
        mutate('is_equal' = ifelse(round(old, 6) == round(new, 6), 1, 0))
      
      # Collect mean values which differ (looks like small rounding differences?)
      tmp_mean_diffs <- tmp_mean_vals |> filter(is_equal == 0) |> mutate('Reason' = 'mean values do not match')

      # Print out if mean values differ when rounded to 6dp
      if (nrow(tmp_mean_diffs) > 0) {
        print(paste("Total mean columns with matching values differing (mismatch):", nrow(tmp_mean_diffs), 'of', nrow(tmp_all_data)))
        if (nrow(tmp_diffs) == 0) {
          tmp_diffs <- tmp_mean_diffs
        } else {
          tmp_diffs <- rbind(tmp_diffs, tmp_mean_diffs)
        }
      } else {
        print("All mean column values found are the same.")
      }

      # Collect mean values where we can't check the difference
      tmp_mean_missing <- tmp_mean_vals |> filter(is.na(is_equal)) |> mutate('Reason' = 'mean values missing in old or new data')
      
      # Print out if mean values differ when rounded to 6dp
      if (nrow(tmp_mean_missing) > 0) {
        print(paste("Total mean columns with matching values differing (missing):", nrow(tmp_mean_missing), 'of', nrow(tmp_all_data)))
        if (nrow(tmp_diffs) == 0) {
          tmp_diffs <- tmp_mean_missing
        } else {
          tmp_diffs <- rbind(tmp_diffs, tmp_mean_missing)
        }
      } else {
        print("No mean column values are missing.")
      }
      
      # Clean up
      rm(list = c('tmp_diff_colnames', 'tmp_old_data', 'tmp_new_data', 'tmp_all_data',
                  'tmp_num_obj', 'tmp_num_obj_diff', 'tmp_mean_vals', 'tmp_mean_diffs', 'tmp_mean_missing'))
    }
  }
  
  # Return data that differs
  return(tmp_diffs)
}


# Build PCA combined plot -------------------------------------------------

build_pca_plots <- function(pca_obj = NULL, 
                            pc1_lims = NULL, pc2_lims = NULL, pc3_lims = NULL, 
                            label_outliers_pc1 = FALSE, label_outliers_pc2 = FALSE, label_outliers_pc3 = FALSE,
                            pc1_label_lims = NULL, pc2_label_lims = NULL, pc3_label_lims = NULL) {
  
  # Create PC data frame from PCA object
  pca_df <- 
    pca_obj$x |>
    as.data.frame() |>
    rownames_to_column('id') |>
    separate(id, into = c('Plate', 'Replicate', 'Well', 'Target', 'Group_Target'), sep = '__', remove = F)
  
  # Get proportion of explained variance per PC
  var_explained_df <- data.frame(PC = paste0("PC", 1:ncol(pca_obj$x)),
                                 var_explained = pca_obj$sdev ^ 2 / sum(pca_obj$sdev ^ 2))
  
  # Scatter plot of PC1 and PC2 (shape = control, color = plate)
  p1_p2 <- 
    ggplot(pca_df, aes(x = PC1, y = PC2, color = Replicate, shape = Target, label = id)) + 
    geom_point(size = 2) + 
    scale_x_continuous(limits = pc1_lims, breaks = pretty_breaks(10)) +
    scale_y_continuous(limits = pc2_lims, breaks = pretty_breaks(10)) +
    theme_pubr()
  
  # Scatter plot of PC1 and PC3 (shape = control, color = plate)
  p1_p3 <- 
    ggplot(pca_df, aes(x = PC1, y = PC3, color = Replicate, shape = Target, label = id)) + 
    geom_point(size = 2) + 
    scale_x_continuous(limits = pc1_lims, breaks = pretty_breaks(10)) +
    scale_y_continuous(limits = pc3_lims, breaks = pretty_breaks(10)) +
    theme_pubr() 
  
  # Scatter plot of PC2 and PC3 (shape = control, color = plate)
  p2_p3 <- 
    ggplot(pca_df, aes(x = PC2, y = PC3, color = Replicate, shape = Target, label = id)) + 
    geom_point(size = 2) + 
    scale_x_continuous(limits = pc2_lims, breaks = pretty_breaks(10)) +
    scale_y_continuous(limits = pc3_lims, breaks = pretty_breaks(10)) +
    theme_pubr()
  
  # Add labels for PC1 outliers to scatter plots
  if (label_outliers_pc1) {
    p1_p2 <- p1_p2 + geom_text_repel(data = subset(pca_df, PC1 > pc1_label_lims[1] | PC1 < pc1_label_lims[2]), 
                                      max.overlaps = 10, size = 3)
    p1_p3 <- p1_p3 + geom_text_repel(data = subset(pca_df, PC1 > pc1_label_lims[1] | PC1 < pc1_label_lims[2]), 
                                      max.overlaps = 10, size = 3) 
  }
  
  # Add labels for PC2 outliers to scatter plots
  if (label_outliers_pc2) {
    p1_p2 <- p1_p2 + geom_text_repel(data = subset(pca_df, PC2 > pc2_label_lims[1] | PC2 < pc2_label_lims[2]), 
                                      max.overlaps = 10, size = 3) 
    p2_p3 <- p2_p3 + geom_text_repel(data = subset(pca_df, PC2 > pc2_label_lims[1] | PC2 < pc2_label_lims[2]), 
                                      max.overlaps = 10, size = 3) 
  }
  
  # Add labels for PC3 outliers to scatter plots
  if (label_outliers_pc3) {
    p1_p3 <- p1_p3 + geom_text_repel(data = subset(pca_df, PC3 > pc3_label_lims[1] | PC3 < pc3_label_lims[2]), 
                                      max.overlaps = 10, size = 3) 
    p2_p3 <- p2_p3 + geom_text_repel(data = subset(pca_df, PC3 > pc3_label_lims[1] | PC3 < pc3_label_lims[2]), 
                                      max.overlaps = 10, size = 3) 
  }
  
  # Scree barplot of proportion of variance explained by each PC
  scree <-
    ggplot(var_explained_df[1:5,], aes(x = PC, y = var_explained)) +
    geom_col() +
    scale_y_continuous(limits = c(0, 1), breaks = pretty_breaks(5)) +
    xlab('Principal component (PC)') +
    ylab('Proportion of variance explained by PC') +
    theme_pubr() 
  
  p <- 
    ggarrange(p1_p2, p1_p3, p2_p3, scree, ncol = 2, nrow = 2, common.legend = TRUE, legend = "top")
  
  return(p)
}


# Convert nested list to data frame --------------------------------------- 

convert_nested_list_to_df <- function(nested_list = NULL) {
  tmp_df <- data.frame()
  for (tmp_plate in names(nested_list)) {
    for (tmp_rep in names(nested_list[[tmp_plate]])) {
      if (0 == nrow(tmp_df)) {
        tmp_df <- nested_list[[tmp_plate]][[tmp_rep]]
      } else {
        tmp_df <- bind_rows(tmp_df, nested_list[[tmp_plate]][[tmp_rep]])
      }
    }
  }
  return(tmp_df)
}


