# Purpose:
# Functions for ingesting and processing imaging data


# Prepare reusable plate labels -------------------------------------------

prepare_plate_labels <- function(path = NULL) {
  # Set labels --------------------------------------------------------------
  
  # Set key value pairs for names and their labels 
  target_names <- c('Control_1' = 'Control', 'Control_2' = 'Control', 'Control_1|Control_2' = 'Control', 'Parental' = 'Control',
                    'ASF1A' = 'ASF1A_ASF1B', 'ASF1B' = 'ASF1A_ASF1B', 'ASF1A|ASF1B' = 'ASF1A_ASF1B', 
                    'CNOT7' = 'CNOT7_CNOT8', 'CNOT8' = 'CNOT7_CNOT8', 'CNOT7|CNOT8' = 'CNOT7_CNOT8',
                    'CCNL1' = 'CCNL1_CCNL2', 'CCNL2' = 'CCNL1_CCNL2', 'CCNL1|CCNL2' = 'CCNL1_CCNL2',
                    'SLC25A37' = 'SLC25A37_SLC25A28', 'SLC25A28' = 'SLC25A37_SLC25A28', 'SLC25A37|SLC25A28' = 'SLC25A37_SLC25A28',
                    'GDI1' = 'GDI1_GDI2', 'GDI2' = 'GDI1_GDI2', 'GDI1|GDI2' = 'GDI1_GDI2',
                    'PDS5A' = 'PDS5A_PDS5B', 'PDS5B' = 'PDS5A_PDS5B', 'PDS5A|PDS5B' = 'PDS5A_PDS5B',
                    'SAR1A' = 'SAR1A_SAR1B', 'SAR1B' = 'SAR1A_SAR1B', 'SAR1A|SAR1B' = 'SAR1A_SAR1B',
                    'SEC23A' = 'SEC23A_SEC23B', 'SEC23B' = 'SEC23A_SEC23B', 'SEC23A|SEC23B' = 'SEC23A_SEC23B',
                    'EAF1' = 'EAF1_EAF2', 'EAF2' = 'EAF1_EAF2', 'EAF1|EAF2' = 'EAF1_EAF2',
                    'INTS6' = 'INTS6_INTS6L', 'INTS6L' = 'INTS6_INTS6L', 'INTS6|INTS6L' = 'INTS6_INTS6L',
                    'TTC7A' = 'TTC7A_TTC7B', 'TTC7B' = 'TTC7A_TTC7B', 'TTC7A|TTC7B' = 'TTC7A_TTC7B')
  
  
  # Read in raw plate data --------------------------------------------------
  
  # Get list of plate label files
  plate_label_files <- list.files(path, pattern = "Plate_.*_target.*", full.names = T, ignore.case = T, recursive = T)
  
  # Create list for plate labels
  plate_label_list <- list()
  
  # Loop over plate data files
  for (i in plate_label_files) {
    # Get plate name from filename
    tmp_plate <- gsub("(.*)_target.*", "\\1", basename(i))
    print(tmp_plate)
    # Read in plate labels
    tmp_data <- suppressMessages(readxl::read_excel(path = i, sheet = 1))
    # Add plate to data frame and rename Well_co to Well
    tmp_data <- tmp_data |> 
      mutate('Plate' = tmp_plate, .before = 'Position') |>
      rename('Well' = 'Well_co')
    # Set group target using target_names and Target columns in data frame containing all plates
    tmp_data <- tmp_data |> 
      rowwise() |>
      mutate('Group_Target' = list(target_names[grep(paste0("^", str_escape(Target), "$"), names(target_names))]), .after = 'Target')
    # Unlist group targets
    tmp_data$Group_Target <- as.character(tmp_data$Group_Target)
    # Set group targets to BLANK when target is BLANK
    tmp_data <- tmp_data |> mutate(Group_Target = ifelse(Target == 'BLANK', 'BLANK', Group_Target))
    # Add to list
    plate_label_list[[tmp_plate]] <- tmp_data
    # Clean up
    rm(list = c(ls(pattern = 'tmp')))
  }
  
  # Clean up 
  rm(i)
  
  # Combine all plate labels into a single data frame
  all_plate_labels <- data.table::rbindlist(list(plate_label_list[[1]],
                                                 plate_label_list[[2]],
                                                 plate_label_list[[3]], 
                                                 plate_label_list[[4]]))
  return(list('all_plate_labels' = all_plate_labels, 
              'plate_label_list' = plate_label_list))
}


# Function to read in raw plate data --------------------------------------

read_raw_plates <- function(files = NULL, plate_label_list = NULL) {
  
  # Check we have files in the vector
  if (0 == length(files) | is.null(files)) { 
    stop('Raw plate files cannot be read: no raw plate files were in the list.') 
  }
  
  # Check we have labels in list
  if (is.null(plate_label_list) | 0 == length(plate_label_list)) { 
    stop('Raw plate files cannot be read: no plate labels were in the list.') 
  }
  
  # Create empty list for raw plate data sets
  list_of_plates <- list()
  
  # Loop over raw plate data files
  for (i in files) {
    # Logging
    print(paste("Reading in raw plate file:", i))
    
    # Get plate name from file name
    tmp_plate <- paste0('Plate_', gsub(".*Plate (.*)__.*", "\\1", i))
    print(paste("Temporary Plate:", tmp_plate))
    
    # Get replicate from file name
    tmp_rep <- ifelse(!grepl('Set', i), 'N1', paste0('N', gsub(".*Screen Set (.*) Plate.*", "\\1", i)))
    print(paste("Temporary Replicate:", tmp_rep))
    
    # Read in plate data to data frame
    tmp_data <- suppressMessages(read_delim(i, delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 9))
    print("Temporary Data (first 6 rows):")
    print(head(tmp_data))
    
    # Remove columns which are all NaN, NA or 0
    columns_to_remove <- names(which(colSums(is.na(tmp_data)) == nrow(tmp_data)))
    print(paste("Columns to remove (contain all NA or NaN):", paste(columns_to_remove, collapse = ', ')))
    tmp_data <- tmp_data |> select(-all_of(columns_to_remove))
    print("Temporary Data after removing columns (first 6 rows):")
    print(head(tmp_data))
    
    # Add plate label
    tmp_data <- tmp_data |> 
      mutate(Plate = tmp_plate, Replicate = tmp_rep, .before = 'Row')
    print("Temporary Data after adding labels (first 6 rows):")
    print(head(tmp_data))
    
    # Merge Row and Column to get Well
    tmp_data <- tmp_data |> mutate(Well = paste(Row, Column, sep = ','))
    print("Temporary Data after creating Well (first 6 rows):")
    print(head(tmp_data))
    
    # Check if plate labels exist
    if (!is.null(plate_label_list[[tmp_plate]])) {
      # Add plate labels
      tmp_data_with_labels <- plate_label_list[[tmp_plate]] |> 
        full_join(tmp_data, by = c('Plate', 'Well', 'Row', 'Column'), multiple = "all")
      print("Temporary Data with labels (first 6 rows):")
      print(head(tmp_data_with_labels))
      
      # Remove wells which are blank (Target == BLANK)
      tmp_data_with_labels <- tmp_data_with_labels |> filter(Target != 'BLANK')
      print("Temporary Data after filtering blanks (first 6 rows):")
      print(head(tmp_data_with_labels))
      
      # Add to list
      list_of_plates[[tmp_plate]][[tmp_rep]] <- tmp_data_with_labels
    } else {
      warning(paste("No labels found for plate:", tmp_plate))
    }
    
    # Clean up
    rm(list = ls(pattern = 'tmp'))
  }
  
  # Clean up 
  rm(i)
  
  # Check list has populated
  if (is.null(list_of_plates) | 0 == length(list_of_plates)) { 
    stop('Raw plate files cannot be read: list is empty or null.') 
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

