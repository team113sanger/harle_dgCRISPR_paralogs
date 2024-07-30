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
#* --                     read_sample_metadata()                           -- *#
#* --                                                                      -- *#
################################################################################

read_sample_metadata <- function(path = NULL) {
  # Check path exists
  check_file_exists(path)

  # Read in sample metadata
  sample_metadata <- read.delim(path, sep = "\t", header = T)
  # Sort cancer type
  sample_metadata$cancer_type <- factor(sample_metadata$cancer_type, levels = sort(unique(sample_metadata$cancer_type)))

  # Add a plot group to allow for separation of controls
  sample_metadata <- sample_metadata %>%
    rowwise() %>%
    mutate('qc_group_with_controls' = ifelse(grepl('control', sample_label), 'control', qc_group))
  sample_metadata$qc_group_with_controls <- factor(sample_metadata$qc_group_with_controls,
                                                   levels = unique(c('control', sample_metadata %>% filter(qc_group_with_controls != 'control') %>% pull(qc_group_with_controls))))
  sample_metadata$qc_group <- factor(sample_metadata$qc_group,
                                     levels = unique(sample_metadata[order(sample_metadata$qc_group_with_controls, sample_metadata$cancer_type), 'qc_group']))
  # Sort sample_metadata
  sample_metadata <- sample_metadata[order(sample_metadata$qc_group_with_controls, sample_metadata$cancer_type),]

  return(sample_metadata)
}

################################################################################
#* --                                                                      -- *#
#* --                     collate_mageck_gene_results()                    -- *#
#* --                                                                      -- *#
################################################################################


collate_mageck_gene_results <- function(path = NULL, sample_mapping = NULL) {
  # Find MAGeCK gene results files
  message('Finding MAGeCK gene results files...')
  mageck_result_files <- list.files(path = path, pattern = 'MAGeCK.gene_summary.txt', full.names = F, recursive = T)
  message(paste('Number of MAGeCK gene results:', length(mageck_result_files)))

  # Create empty data frame
  mageck_gene_results <- data.frame()

  # Loop over MAGeCK gene results files
  message('Reading and processing MAGeCK gene results files...')
  for (result_file in mageck_result_files) {
    # Get cell line label
    cl <- str_split(result_file, pattern = "/", simplify = T)[1]
    cl <- sample_mapping %>% filter(stripped_cell_line_name == cl) %>% pull(cell_line_label) %>% unique()
    # Get data set
    d <- str_split(result_file, pattern = "/", simplify = T)[2]
    # Read in results
    mageck_tmp <- read.delim(file.path(path, result_file), header = T, check.names = T)
    # Add data set and cell line label to result
    mageck_tmp <- mageck_tmp %>%
      mutate('dataset' = d, 'cell_line_label' = cl, .before = 'id') %>%
      rename('gene' = 'id')
    # Add binary for significantly depleted and enriched genes
    mageck_tmp <- mageck_tmp %>%
      mutate('is_depleted_mageck' = ifelse(neg.fdr < 0.05, 1, 0),
             'is_enriched_mageck' = ifelse(pos.fdr < 0.05, 1, 0))
    # Add to existing dataset results
    if (nrow(mageck_gene_results) == 0) {
      mageck_gene_results <- mageck_tmp
    } else {
      mageck_gene_results <- bind_rows(mageck_gene_results, mageck_tmp)
    }
  }
  message('MAGeCK gene results collated...')
  return(mageck_gene_results)
}



################################################################################
#* --                                                                      -- *#
#* --                      collate_bagel_gene_results()                    -- *#
#* --                                                                      -- *#
################################################################################

collate_bagel_gene_results <- function(path = NULL, sample_mapping = NULL) {
  # Find BAGEL2 guide results files
  message('Finding BAGEL2 guide results files...')
  BAGEL2_result_files <- list.files(path = path, pattern = 'BAGEL2.sgrna.bf', full.names = F, recursive = T)
  message(paste('Number of BAGEL2 guide results:', length(BAGEL2_result_files)))

  # Create empty data frame
  BAGEL2_gene_results <- data.frame()

  # Loop over BAGEL2 guide results files
  message('Reading and processing BAGEL2 guide results files...')
  for (result_file in BAGEL2_result_files) {
    # Get cell line label
    cl <- str_split(result_file, pattern = "/", simplify = T)[1]
    cl <- sample_mapping %>% filter(stripped_cell_line_name == cl) %>% pull(cell_line_label) %>% unique()
    # Get data set
    d <- str_split(result_file, pattern = "/", simplify = T)[2]
    # Read in results
    BAGEL2_tmp <- read.delim(file.path(path, result_file), header = T, check.names = T)
    # Remove calculate mean BFs per cell line
    BAGEL2_tmp <- BAGEL2_tmp %>%
      group_by(GENE) %>%
      summarise('BF' = mean(BF)) %>%
      mutate('dataset' = d, 'cell_line_label' = cl, .before = 'GENE')
    # Add to existing dataset results
    if (nrow(BAGEL2_gene_results) == 0) {
      BAGEL2_gene_results <- BAGEL2_tmp
    } else {
      BAGEL2_gene_results <- bind_rows(BAGEL2_gene_results, BAGEL2_tmp)
    }
  }
  message('BAGEL2 gene results collated...')
  return(BAGEL2_gene_results)
}

################################################################################
#* --                                                                      -- *#
#* --                             roc_metrics()                            -- *#
#* --                                                                      -- *#
################################################################################

roc_metrics <- function(data = NULL) {
  datasets <- c('A', 'B', 'combined')
  roc_list <- list()
  cell_lines <- unique(data$cell_line_label)

  for (d in datasets) {
    for (cl in cell_lines) {
      # Prepare data
      data.filt <- data %>%
        filter(cell_line_label == cl & dataset == d) %>%
        select('gene' = 'GENE', values = 'BF', classification)

      # Get all classified genes (essential and non-essential)
      all_genes <- data.filt %>%
        filter(classification != 'unknown') %>% pull(gene) %>% unique()

      # Get essential genes
      essential_genes <- data.filt %>%
        filter(classification == 'Essential') %>% pull(gene) %>% unique()

      # Prepare the response dataframe for pROC
      # Essential = 1, Non-essential = 0
      essentiality <- rep(0, length(all_genes))
      names (essentiality) <- all_genes
      essentiality[essential_genes] <- 1

      # Get predicted values for each observation
      essentiality_data <- list()
      essentiality_data[['essentiality']] <- essentiality
      essentiality_data[['predictor']] <- data.filt %>% filter(gene %in% all_genes) %>% pull(values)
      essentiality_data[['min']] <- min(essentiality_data[['predictor']], na.rm = TRUE)
      essentiality_data[['modified_predictor']] <-  essentiality_data[['predictor']] - essentiality_data[['min']]

      # Preparing ROC
      roc_obj <- suppressMessages(roc(essentiality_data[['essentiality']], essentiality_data[['modified_predictor']]))
      roc_coords <- pROC::coords(roc_obj, ret = c( 'all' ), transpose = T)
      roc_coords['threshold',] <- roc_coords['threshold',] + essentiality_data[['min']]

      # Add id
      id <- min(which(roc_coords['ppv',] > (1 - 0.05)))
      if( id == "Inf" ){
        id <- min(which(roc_coords['ppv',] >= (1 - 0.05)))
      }
      if( id == "Inf" ){
        id <- max(which(round(roc_coords['ppv',]) >= (1 - 0.05)))
      }

      # Get best precision
      best_prec <- c(roc_coords['threshold',id],
                     roc_coords['specificity',id],
                     roc_coords['sensitivity',id],
                     roc_coords['ppv',id])

      # Get thresholds
      bestPrecisionTh <- rbind(best_prec)
      colnames(bestPrecisionTh) <- c('thresholds', 'specificity', 'sensitivity', 'ppv')

      # Add data set and cell line label
      bestPrecisionTh <- bestPrecisionTh %>%
        as_tibble() %>%
        mutate('dataset' = d, 'cell_line_label' = cl, 'AUC' = roc_obj$auc[1], .before = 'thresholds')

      # ROC list
      roc_coords.df <- data.frame('dataset' = d,
                                  'cell_line_label' = cl,
                                  'FPR' = (1 - roc_obj$specificities),
                                  'TPR' = roc_obj$sensitivities)
      if ('roc_coords' %in% names(roc_list)) {
        roc_list[['roc_coords']] <- rbind(roc_list[['roc_coords']], roc_coords.df)
        roc_list[['thresholds']] <- rbind(roc_list[['thresholds']], bestPrecisionTh)
      } else {
        roc_list[['roc_coords']] <- roc_coords.df
        roc_list[['thresholds']] <- bestPrecisionTh
      }
    }
  }
  return(roc_list)
}

################################################################################
#* --                                                                      -- *#
#* --                             plot_roc()                               -- *#
#* --                                                                      -- *#
################################################################################

# Stacked plot of read pair mapping proportions
plot_roc <- function(data = NULL, pal = NULL) {
  p <-
    ggplot(data, aes(x = FPR, y = TPR, group = cell_line_label, color = cancer_type)) +
    geom_line(linewidth = 0.5) +
    scale_color_manual(values = pal, name = 'Cancer type') +
    facet_grid(dataset ~ cancer_type) +
    xlab('FPR (1 - Specificity)') +
    ylab('TPR (Sensitivity)') +
    theme_pubr() +
    theme(text = element_text(size = 12),
          panel.border = element_blank(),
          strip.background = element_blank())
  return(p)
}

################################################################################
#* --                                                                      -- *#
#* --                         plot_n_genes_binary_depleted()               -- *#
#* --                                                                      -- *#
################################################################################

plot_n_genes_binary_depleted <- function(data = NULL, pal = NULL) {
  data.summary <- data %>%
    select(cell_line_label, cancer_type, targetA__is_single_depleted, targetB__is_single_depleted) %>%
    gather(classification, pass, -cell_line_label, -cancer_type) %>%
    group_by(cell_line_label, cancer_type, classification) %>%
    summarise(.groups = 'keep', n = sum(pass == 1))
  
  p <- 
    ggplot(data.summary, aes(x = cell_line_label, y = n, fill = cancer_type)) +
      geom_bar(stat = 'identity') +
      facet_grid(classification ~ .) +
      scale_fill_manual(values = pal) +
      scale_y_continuous(breaks = pretty_breaks(10)) +
      theme_pubr(base_size = 14) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
            panel.border = element_blank(),
            strip.background = element_blank()) +
      labs(fill = "Cancer type",  x = '', y = 'Number of gene pairs')
  
  return(p)
}




