# Check file exists
check_file_exists <- function(filepath) {
  if (!file.exists(filepath)) {
    stop(sprintf("File does not exist: %s", filepath))
  }
}

# Check directory exists
check_dir_exists <- function(dirpath) {
  if (!dir.exists(dirpath)) {
    stop(sprintf("Directory does not exist: %s", dirpath))
  }
}

# Create directory
create_directory <- function(dirpath) {
  if (!dir.exists(dirpath)) {
    dir.create(dirpath, recursive = T) # Create directory if it doesn't exist
    check_dir_exists(dirpath)
  }
  return(TRUE)
}

################################################################################
#* --                                                                      -- *#
#* --                     convert_variable_to_integer()                    -- *#
#* --                                                                      -- *#
################################################################################

#' Converts variable value to integers
#'
#' @description
#' Converts variable value to integer (where necessary)
#'
#' @param x input value
#'
#' @export convert_variable_to_integer
convert_variable_to_integer <-
  function(x = NULL) {
    y <- x
    # Validate input
    if (is.null(y) || length(y) == 0)
      stop("Cannot convert variable to integer, x is null.")
    # Only assign variable if it's not null
    if (!is.null(y)) {
    # Skip if the variable is already an integer
      if (!is.integer(y)) {
        y <- tryCatch({
          strtoi(x)
        }, error = function(e) {
          # Stop if there is an error
          stop(paste("Cannot make variable an integer:", x))
        })
      }
      if (is.na(y))
        stop(paste("Could not convert value to integer:", x))
    }
    return(y)
  }

###############################################################################
#* --                                                                     -- *#
#* --                           check_dataframe()                         -- *#
#* --                                                                     -- *#
###############################################################################

#' Check dataframe is not empty and indices are valid
#'
#' @description Check dataframe is not empty and indices exist.
#' @export check_dataframe
#' @param data a dataframe.
#' @param indices vector of indices to check.
#' @param check_na whether to check for NA values.
#' @param check_nan whether to check for NaN values.
check_dataframe <-
  function(data = NULL, indices = NULL, check_na = FALSE, check_nan = FALSE) {
    if (is.null(data))
      stop("Dataframe is null.")
    if (nrow(data) == 0) # Assumes ncol also by default
      stop("Dataframe has no rows.")
    # Check indices don't exceed dataframe columns
    if (!is.null(indices)) {
      for (i in 1:length(indices)) {
        indices[i] <- tryCatch({
          convert_variable_to_integer(indices[i])
        }, error = function(e) {
          # Stop if there is an error
          stop(paste("Cannot convert index to integer:", indices[i]))
        })
        if (indices[i] > ncol(data))
          stop(paste0("Index exceeds dataframe limits: ", indices[i], ", ", ncol(data)))
      }
    }
    # Check whether data frame contains NaN values
    if (check_nan) {
      if (sum(apply(data, 2, function(x) any(is.nan(x)))) != 0)
        stop("Dataframe contains NaN values.")
    }
    # If required, check whether data frame contains NA values
    if (check_na) {
      if (sum(apply(data, 2, function(x) any(is.na(x)))) != 0)
        stop("Dataframe contains NA values.")
    }
    return(TRUE)
  }

###############################################################################
#* --                                                                     -- *#
#* --                 check_is_numeric_and_is_integer()                   -- *#
#* --                                                                     -- *#
###############################################################################

#' Check value is numeric and is an integer.
#'
#' @description Check value is numeric and is an integer.
#' @export check_is_numeric_and_is_integer
#' @param x a value.
#' @returns logical.
check_is_numeric_and_is_integer <-
  function(x) {
    if (is.null(x)) {
      return(FALSE)
    } else if (is.na(x)) {
      return(FALSE)
    } else if (!is.numeric(x)) {
      return(FALSE)
    } else if (x %% 1 != 0) {
      return(FALSE)
    } else {
      return(TRUE)
    }
  }

###############################################################################
#* --                                                                     -- *#
#* --                       process_column_indices()                      -- *#
#* --                                                                     -- *#
###############################################################################

#' Process column indices into vector.
#'
#' @description
#' Convert a character string of column indices into a numeric vector.
#'
#' @examples
#' process_column_indices('2,3,4,5')
#' process_column_indices('2:5')
#' process_column_indices('2-5')
#' process_column_indices('2:4,5')
#'
#' @param columns character string of column indices.
#'
#' @return a numeric vector of column indices.
#' @import dplyr
#' @importFrom stringr str_split
#' @export process_column_indices
process_column_indices <-
  function(columns = NULL) {
    # Error if there are no values or input is null
    if (length(columns) < 1 || is.null(columns))
      stop("Cannot process columns (<1 or null).")
    # Split character into a vector
    columns <- columns %>% paste(collapse = ",") %>% str_split(',', simplify = T)
    # Loop through values and expand ranges
    processed_columns <- vector()
    for (i in 1:length(columns)) {
      # If the input is a range, expand the range
      if (grepl("[:-]", columns[i])) {
        # Split the start and end values
        numeric_column_range <- suppressWarnings(
          as.numeric(
            str_split(columns[i], "[\\:\\-]", n = 2, simplify = T)))
        # Check the range values are integers
        if (!check_is_numeric_and_is_integer(numeric_column_range[1]) ||
            !check_is_numeric_and_is_integer(numeric_column_range[2]))
          stop(paste("Column indices contain a non-integer value in range:", columns[i]))
        # Add expanded range into the column vector
        processed_columns <- c(processed_columns,
                               c(numeric_column_range[1]:numeric_column_range[2]))
      } else {
        # If not a range, check it's an integer and return the value
        numeric_column_index <- suppressWarnings(as.numeric(columns[i]))
        # Check value is an integer
        if (!check_is_numeric_and_is_integer(numeric_column_index)) {
          stop(paste("Column indices contain a non-integer value:", columns[i]))
        }
        # Add value into column vector
        processed_columns <- c(processed_columns, numeric_column_index)
      }
    }
    return(processed_columns)
  }

###############################################################################
#* --                                                                     -- *#
#* --                    get_guides_failing_filter()                      -- *#
#* --                                                                     -- *#
###############################################################################
#' Identify guides failing count filter.
#'
#' @description
#' Identify guides which fail count filter.
#'
#' @param count_matrix count matrix.
#' @param id_column index of column containing unique sgRNA identifiers.
#' @param count_column indices of columns containing counts.
#' @param filter_indices indices of columns on which to apply filter.
#' @param filter_method filter method from one of: all, any, mean or median (Default: all).
#' @param min_reads minimum number of reads for filter (Default: 30).
#'
#' @import dplyr
#' @importFrom matrixStats rowMedians
#' @return dataframe
#' @export get_guides_failing_filter
get_guides_failing_filter <-
  function(count_matrix = NULL,
           id_column = 1,
           count_column = NULL,
           filter_indices = NULL,
           filter_method = 'all',
           min_reads = 30) {
    # Check input data is not null
    if (is.null(count_matrix))
      stop("Cannot get guides to filter from count matrix, count matrix is null.")
    if (is.null(count_column))
      stop("Cannot get guides to filter from count matrix, count_column is null.")
    if (is.null(min_reads))
      stop("Cannot get guides to filter from count matrix, min_reads is null.")
    if (min_reads < 0)
      stop(paste("Cannot get guides to filter from count matrix, min_reads is < 0:", min_reads))
    if (is.null(filter_indices))
      stop("Cannot get guides to filter from count matrix, filter_indices is null.")
    if (is.null(filter_method))
      stop("Cannot get guides to filter from count matrix, filter_method is null.")
    if (!filter_method %in% c('all', 'any', 'mean', 'median'))
      stop(paste("Cannot get guides to filter from count matrix, filter_method is not valid (all, any, mean or median):",
                 filter_method))
    # Try to make each column an integer if it isn't already
      for (i in c('id_column', 'min_reads')) {
        if (!is.null(get(i)))
          assign(i, convert_variable_to_integer(get(i)))
      }
    # Check count matrix
    check_dataframe(count_matrix)
    # Expand count and filter indices
    count_column <- process_column_indices(count_column)
    filter_indices <- process_column_indices(filter_indices)
    # Check that filter indices are in count columns
    if (length(setdiff(filter_indices, count_column)) != 0)
      stop("Cannot get guides to filter from count matrix, filter indices not in count columns.")
    # Prepare count matrix
    count_matrix <- count_matrix[,c(id_column, filter_indices)]
    # Set up empty vector for filtered reads
    id_column_name <- colnames(count_matrix)[1]
    if (filter_method == 'all') {
      filtered_guides <- data.frame(count_matrix[,1],
        apply(count_matrix[2:ncol(count_matrix)], 2, function(x) x < min_reads))
      filtered_guides <- filtered_guides[,1]
      filtered_guides <- count_matrix %>%
                         filter_at(vars(-!!id_column_name), all_vars(. < min_reads)) %>%
                          pull(!!id_column_name)
    } else if (filter_method == 'any') {
      filtered_guides <- count_matrix %>%
                          filter_at(vars(-!!id_column_name), any_vars(. < min_reads)) %>%
                          pull(!!id_column_name)
    } else if (filter_method == 'mean') {
      mean_cols <- colnames(count_matrix)[-1]
      filtered_guides <- count_matrix %>%
                          mutate(filter_mean = rowMeans(as.matrix(.[mean_cols]), na.rm = FALSE)) %>%
                          filter(filter_mean < min_reads ) %>%
                          pull(!!id_column_name)
    } else if (filter_method == 'median') {
      median_cols <- colnames(count_matrix)[-1]
      filtered_guides <- count_matrix %>%
                          mutate(filter_median = rowMedians(as.matrix(.[median_cols]), na.rm = FALSE)) %>%
                          filter(filter_median < min_reads) %>%
                          pull(!!id_column_name)
    } else {
      filtered_guides <- vector()
    }
    return(filtered_guides)
  }

###############################################################################
#* --                                                                     -- *#
#* --                            calculate_lfc()                          -- *#
#* --                                                                     -- *#
###############################################################################

#' Calculate log fold changes.
#'
#' @description Calculate log fold changes.
#'
#' @details
#' Takes the per-guide mean of control samples, adds a pseduocount to
#' all counts (including control mean) and calculates the log2 fold change.
#' Requires column indices be defined:
#' \itemize{
#' \item `id_column` - column containing guide (sgRNA) identifiers (Default = 1).
#' \item `gene_column` - column containing gene symbols/identifiers (Default = 2).
#' \item `control_indices` - indices of columns containing control sample counts.
#' \item `treatment_indices` - indices of columns containing treatment sample counts.
#' }
#'
#' @param data sample count matrix.
#' @param id_column the index of column containing unique sgRNA identifiers.
#' @param gene_column the index of column containing gene symbols.
#' @param control_indices vector indices of columns containing control sample counts.
#' @param treatment_indices vector indices of columns containing treatment sample counts.
#' @param pseudocount a pseudocount (Default: 0.5).
#'
#' @return a data frame containing sample log fold changes.
#' @export calculate_lfc
calculate_lfc <-
  function(data = NULL,
           id_column = 1,
           gene_column = 2,
           control_indices = NULL,
           treatment_indices = NULL,
           pseudocount = 0.5) {
    # Validate input
    if (is.null(data))
      stop("Cannot calculate LFCs, data is null.")
    if (is.null(control_indices))
      stop("Cannot calculate LFCs, control_indices is null.")
    if (is.null(treatment_indices))
      stop("Cannot calculate LFCs, treatment_indices is null.")
    if (is.null(pseudocount))
      stop("Cannot calculate LFCs, pseudocount is null.")
    if (!is.numeric(pseudocount))
      stop("Cannot calculate LFCs, pseudocount is not numeric.")
    check_is_numeric_and_is_integer(id_column)
    check_is_numeric_and_is_integer(gene_column)
    # Check data
    check_dataframe(data, check_na = F, check_nan = T)
    # Process control and treatment indices
    control_indices <- process_column_indices(control_indices)
    treatment_indices <- process_column_indices(treatment_indices)
    count_indices <- c(control_indices, treatment_indices)
    # Check that control and treatment indices don't overlap
    if (sum(duplicated(count_indices)) > 0)
      stop(paste("Cannot calculate LFC, duplicated indices:",
                 paste(control_indices, sep = ','),
                 paste(treatment_indices, sep = ',')))
    # Get only the data we need
    data <- data[,c(id_column, gene_column, count_indices)]
    # Reformat control indices and get number of samples
    ncontrol <- length(process_column_indices(control_indices))
    control_indices <- c(3:(3 + (ncontrol - 1)))
    ntreatment <- length(process_column_indices(treatment_indices))
    treatment_indices <- c((2 + ncontrol + 1):(2 + ncontrol + ntreatment))
    nsamples <- ncontrol + ntreatment
    # Get mean of control indices
    data <- data.frame(data[1:2],
                       'control_means' = apply(data[3:(4 + (ncontrol - 2))], 1, function(x) mean(x)),
                       data[treatment_indices], check.names = FALSE)
    # Add pseudocount to counts
    data <- add_pseudocount(data,
                            pseudocount = pseudocount,
                            indices = c(3:(3 + ntreatment)))
    # Calculate log fold changes
    lfc <- data.frame(data[1:2],
                      log2(data[,4:(3 + ntreatment)] / data[,3]), check.names = FALSE)
    # Preserve column names when only one treatment sample
    if (length(treatment_indices) == 1) {
      colnames(lfc)[3] <- colnames(data)[4]
    }
    return(lfc)
  }

###############################################################################
#* --                                                                     -- *#
#* --                         add_pseudocount()                           -- *#
#* --                                                                     -- *#
###############################################################################

#' Add pseudocount
#'
#' @description Add pseudocount to selected dataframe columns
#' @export add_pseudocount
#' @param data a data frame.
#' @param pseudocount a pseudocount (Default: 5).
#' @param indices column indices.
#' @param ... parameters for `check_dataframe`.
#'
#' @import dplyr
#' @return a data frame.
#' @export add_pseudocount
add_pseudocount <-
  function(data = NULL,
           pseudocount = 5,
           indices = NULL,
           ...) {
    # Validate inputs
    if (is.null(data))
      stop("Cannot add pseudocount, data is null.")
    if (is.null(pseudocount))
      stop("Cannot add pseudocount, pseudocount is null.")
    if (!is.numeric(pseudocount))
      stop("Cannot add pseudocount, pseudocount is not numeric.")
    # Check dataframe
    check_dataframe(data)
    if (is.null(indices)) {
      # If no indices given, apply to full data frame
      # Add pseudocount to all data frame columns
      warning("No indices given, adding pseudocount to all indices")
      data <- data %>% mutate(across(all_of(indices), ~ . + pseudocount))
    } else {
      # Process column indices
      indices <- process_column_indices(indices)
      # Add pseudocount to selected data frame columns
      data <- tryCatch({
        data %>% mutate(across(all_of(indices), ~ . + pseudocount))
      }, error = function(e) {
        # Stop if there is an error
        stop("Cannot add pseudocount to dataframe.")
      })
    }
    # Check data frame
    check_dataframe(data, ...)
    # Return data frame
    return(data)
  }

###############################################################################
#* --                                                                     -- *#
#* --                         get_guide_matrix()                          -- *#
#* --                                                                     -- *#
###############################################################################

get_guide_matrix <- function(lb = NULL, control_guide_file = NULL) {
  message("Prepare empty dual guide matrix...")
  
  # Check control file exists and read in control guide names (one per line)
  check_file_exists(control_guide_file)
  control_guides <- scan(file = control_guide_file, what = 'character')
  
  # Empty guide matrix
  guide_matrix <- data.frame(id = character(), sorted_gene_pair = character(), g1 = character(), g2 = character())
  
  # Get only gene|gene guides
  message("Getting dual guides...")
  message(paste("Total number of guides:", nrow(lb)))
  
  # Remove safe targeting controls
  df <- lb %>%
    filter(guide_type != 'safe_targeting|safe_targeting')
  message(paste("Number of guides without controls:", nrow(df)))
  
  # Get dual guide genes
  duals <- df %>% filter(guide_type == 'gene|gene')
  message(paste("Number of dual guides:", nrow(duals)))
  
  # Get single guide genes (A)
  singlesA <- df %>% filter(guide_type == 'gene|safe_targeting')
  message(paste("Number of single guides (A):", nrow(singlesA)))
  
  # Get single guide genes (B)
  singlesB <- df %>% filter(guide_type == 'safe_targeting|gene')
  message(paste("Number of single guides (B):", nrow(singlesB)))
  
  # Loop over gene|gene guides and get corresponding singles
  for (i in 1:nrow(duals)) {
    # Get dual information
    gp <- duals$sorted_gene_pair[i]
    g12 <- duals$id[i]
    
    # Get left guide id
    g1 <- singlesA %>% 
      filter(sorted_gene_pair == gp) %>% # Make sure we're looking at the correct gene pair
      separate(sgrna_ids, into = c('l_guide', 'r_guide'), sep = "\\|", remove = F) %>% # Split the guide ids
      filter(l_guide == str_split(g12, "\\|", simplify = T)[1]) %>% # Make sure we match the correct left guide
      select(id) %>% unlist() %>%  as.vector()
    
    # Get right guide id
    g2 <- singlesB %>% 
      filter(sorted_gene_pair == gp) %>% # Make sure we're looking at the correct gene pair
      separate(sgrna_ids, into = c('l_guide', 'r_guide'), sep = "\\|", remove = F) %>% # Split the guide ids
      filter(r_guide == str_split(g12, "\\|", simplify = T)[2]) %>% # Make sure we match the correct left guide
      select(id) %>% unlist() %>%  as.vector()
 
    # Add guide pair to data frame
    guide_matrix.tmp <- data.frame(id = g12, sorted_gene_pair = gp, g1 = g1, g2 = g2)
    guide_matrix <- rbind(guide_matrix, guide_matrix.tmp)
  }

  message('Returning guide matrix...')
  return(guide_matrix)
}