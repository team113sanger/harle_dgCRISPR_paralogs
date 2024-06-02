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
#* --               plot_n_bassik_hits_per_cell_line()                     -- *#
#* --                                                                      -- *#
################################################################################

# Stacked plot of read pair mapping proportions
plot_n_bassik_hits_per_cell_line <- function(data = NULL, pal = NULL) {
  p <-
    ggplot(data, aes(x = cell_line_label, y = n, fill = cancer_type)) +
      geom_bar(stat = 'identity') +
      scale_y_continuous(breaks = pretty_breaks(10)) +
      scale_fill_manual(values = pal) +
      theme_pubr(base_size = 14) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10)) +
      labs(fill = "Cancer Type", x = '', y = 'Total number of Bassik hits')
  return(p)
}

################################################################################
#* --                                                                      -- *#
#* --               plot_stacked_bassik_hits_per_cell_line()               -- *#
#* --                                                                      -- *#
################################################################################

plot_stacked_bassik_hits_per_cell_line <- function(data = NULL, pal = NULL) {
  p <-
    ggplot(data, aes(y = reorder(sorted_gene_pair, n), x = n, fill = cancer_type)) +
      geom_bar(position="stack", stat="identity") +
      scale_fill_manual(values = pal) +
      scale_x_continuous(position = "top") +
      theme_pubr(base_size = 14) +
      theme(axis.text.y = element_text(size = 10)) +
      labs(fill = "Cancer Type", y = '', x = 'Number of Bassik hits')
  return(p)
}

################################################################################
#* --                                                                      -- *#
#* --                 plot_bassik_hits_mean_gi_per_cell_line()             -- *#
#* --                                                                      -- *#
################################################################################

plot_bassik_hits_mean_gi_per_cell_line <- function(data = NULL, pal = NULL) {
  p <-
    ggplot(bassik_mean_gi, aes(x = mean_norm_gi, y = reorder(sorted_gene_pair, -mean_norm_gi), color = cancer_type)) +
      geom_point(size = 2, alpha = 0.6) +
      geom_vline(xintercept = -0.5, linetype = 'dotted') +
      scale_color_manual(values = cancer_type_pal) +
      scale_x_continuous(breaks = pretty_breaks(10), limits = c(-4, 0)) +
      theme_pubr(base_size = 14) +
      theme(axis.text.y = element_text(size = 10)) +
      labs(fill = "Cancer Type", y = '', x = 'Mean norm GI')
  return(p)
}