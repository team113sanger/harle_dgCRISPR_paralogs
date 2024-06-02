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
#* --                     collate_pyc_stats_json()                         -- *#
#* --                                                                      -- *#
################################################################################

# Collates JSON statistics files produced by pyCROQUET

collate_pyc_stats_json <- function(path = NULL, sample_metadata = NULL) {
  # Get list of pyCROQUET statistics JSON files
  pyc_stats_files <- list.files(path = path, pattern = "*.stats.json", full.names = T)
  # Create empty dataframe
  pycroquet_lane_stats <- data.frame()

  # Loop over the statistics files
  for (fn in pyc_stats_files) {
    # Read in JSON
    tmp_json <- jsonlite::fromJSON(fn, flatten = T)
    # Extract pair classifications from nested JSON
    tmp_classes <- unlist(tmp_json$pair_classifications)
    # Get only required columns
    tmp_json <- unlist(tmp_json[c('sample_name', 'total_reads', 'total_pairs', 'mapped_to_guide_reads', 'unmapped_reads')])
    # Combine into a data frame
    tmp_data <- data.frame(t(c(tmp_json, tmp_classes)))
    # Modify data frame
    tmp_data <- tmp_data %>%
      mutate_at(vars(-sample_name), as.numeric)

    # Add to dataframe for all samples
    if ( 0 == nrow(pycroquet_lane_stats)) {
      pycroquet_lane_stats <- tmp_data
    } else {
      pycroquet_lane_stats <- rbind(pycroquet_lane_stats, tmp_data)
    }
  }
  # Add sample metadata
  pycroquet_lane_stats <- pycroquet_lane_stats %>%
    left_join(sample_metadata, by = c('sample_name' = 'sanger_sample_name')) %>%
    rename(sanger_sample_name = sample_name)
  return(pycroquet_lane_stats)
}

################################################################################
#* --                                                                      -- *#
#* --                     get_pyc_sample_stats()                           -- *#
#* --                                                                      -- *#
################################################################################

# Summarise pyCROQUET statistics by sample

get_pyc_sample_stats <- function(data) {
  sample_data <- pycroquet_statistics %>%
    group_by(sanger_sample_name, sample_label, cell_line_label, cancer_type, qc_group, qc_group_with_controls, replicate) %>%
    summarise('total_reads' = sum(total_reads),
              'total_pairs' = sum(total_pairs),
              'mapped_to_guide_reads' = sum(mapped_to_guide_reads),
              'mapped' = sum(match),
              'unmapped' = sum(f_multi_3p, f_multi_5p, f_open_3p, f_open_5p, r_multi_3p, r_multi_5p, r_open_3p, r_open_5p, aberrant_match, ambiguous, no_match, swap),
              'swap' = sum(swap), .groups = 'keep') %>%
    ungroup()
  return(sample_data)
}

################################################################################
#* --                                                                      -- *#
#* --                     get_pyc_sample_stats()                           -- *#
#* --                                                                      -- *#
################################################################################

# Bar plot of total read pairs per sample
plot_total_read_pairs <- function(data = NULL) {
  p <-
    ggplot(data, aes(x = sample_label, y = total_pairs, fill = cancer_type, color = cancer_type)) +
      geom_bar(stat = 'identity') +
      scale_y_continuous(labels = label_number(suffix = " M", scale = 1e-6), breaks = pretty_breaks(10)) +
      scale_fill_npg(alpha = 0.3) +
      scale_color_npg(alpha = 0.6) +
      facet_grid(. ~ qc_group_with_controls, scales = 'free_x', space = 'free') +
      theme_pubr(base_size = 14) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
            panel.border = element_blank(),
            strip.background = element_blank(),
            strip.text = element_blank()) +
      labs(color = "Cancer Type", fill = "Cancer Type", x = '', y = 'Number of read pairs (millions)')
  return(p)
}

################################################################################
#* --                                                                      -- *#
#* --                 prepare_mapping_statistics()                         -- *#
#* --                                                                      -- *#
################################################################################

# Prepare data frame of mapping statistics from pyCROQUET results
prepare_mapping_statistics <- function(data = NULL, sample_metadata = NULL) {
  mapping_statistic <- data %>%
    select(sample_label, cancer_type, qc_group_with_controls, mapped, unmapped) %>%
    gather(category, n_pairs, -sample_label, -cancer_type, -qc_group_with_controls)
  mapping_statistic$sample_label <- factor(mapping_statistic$sample_label, levels = sample_metadata$sample_label)
  mapping_statistic$category <- factor(mapping_statistic$category, levels = c('unmapped', 'mapped'))
  return(mapping_statistic)
}

################################################################################
#* --                                                                      -- *#
#* --                     plot_mapped_read_pair()                          -- *#
#* --                                                                      -- *#
################################################################################

# Plot number of mapped and unmapped reads
plot_mapped_read_pair <- function(data = NULL) {
  p <-
    ggplot(data, aes(x = sample_label, y = n_pairs, fill = category)) +
    geom_col(color = 'gray5', linewidth = 0.1, alpha = 0.8) +
    scale_y_continuous(labels = label_number(suffix = " M", scale = 1e-6), breaks = pretty_breaks(10)) +
    scale_fill_brewer(palette = "Blues") +
    facet_grid(. ~ qc_group_with_controls, scales = 'free_x', space = 'free') +
    theme_pubr(base_size = 14) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
          panel.border = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank()) +
    labs(color = "Category",  x = '', y = 'Number of read pairs (millions)')
  return(p)
}

################################################################################
#* --                                                                      -- *#
#* --                     get_library_statistics()                         -- *#
#* --                                                                      -- *#
################################################################################

# Calculate library statistics
get_library_statistics <- function(data = NULL) {
 stats <- data %>%
  group_by(sample_label) %>%
  summarise(.groups = 'keep',
            'n' = n(),
            'min' = min(norm_count),
            'max' = max(norm_count),
            'median' = median(norm_count),
            'mean' = mean(norm_count),
            'low_counts' = sum(norm_count < opt$low),
            'gini_index' = calculate_gini_index(norm_count))
  return(stats)
}

################################################################################
#* --                                                                      -- *#
#* --                     calculate_gini_index()                           -- *#
#* --                                                                      -- *#
################################################################################

# Calculate Gini index
# Adapted from https://github.com/cran/ineq/blob/master/R/ineq.R
calculate_gini_index <- function(x = NULL) {
  x <- as.numeric(na.omit(x))
  n <- length(x)
  x <- sort(x)
  G <- sum(x * 1L:n)
  G <- 2 * G/sum(x) - (n + 1L)
  gini <- G/n
  return(gini)
}

################################################################################
#* --                                                                      -- *#
#* --               plot_mapped_read_pair_proportion()                     -- *#
#* --                                                                      -- *#
################################################################################

# Stacked plot of read pair mapping proportions
plot_mapped_read_pair_proportion <- function(data = NULL) {
  p <-
    ggplot(data, aes(x = sample_label, y = n_pairs, fill = category)) +
      geom_col(position = 'fill', color = 'gray5', linewidth = 0.1, alpha = 0.8) +
      scale_y_continuous(labels = scales::percent, breaks = pretty_breaks(10)) +
      scale_fill_brewer(palette = "Blues") +
      facet_grid(. ~ qc_group_with_controls, scales = 'free_x', space = 'free') +
      theme_pubr(base_size = 14) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
            panel.border = element_blank(),
            strip.background = element_blank(),
            strip.text = element_blank()) +
      labs(color = "Category",  x = '', y = 'Proportion of read pairs (%)')
  return(p)
}

################################################################################
#* --                                                                      -- *#
#* --                  plot_library_statistic()                            -- *#
#* --                                                                      -- *#
################################################################################

# Generic plotting function for library statistics
plot_library_statistic <- function(data = NULL, variable = NULL, ylab = NULL) {
  p <-
    ggplot(data, aes(x = sample_label, y = .data[[variable]], fill = cancer_type, color = cancer_type)) +
    geom_bar(stat = 'identity') +
    scale_y_continuous(breaks = pretty_breaks(10)) +
    scale_fill_npg(alpha = 0.3) +
    scale_color_npg(alpha = 0.6) +
    facet_grid(. ~ qc_group_with_controls, scales = 'free_x', space = 'free') +
    theme_pubr(base_size = 14) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
          panel.border = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank()) +
    labs(color = "Cancer Type", fill = "Cancer Type", x = '', y = ylab)
  return(p)
}

################################################################################
#* --                                                                      -- *#
#* --                  plot_control_correlation()                          -- *#
#* --                                                                      -- *#
################################################################################

plot_control_correlation <- function(data = NULL, path = NULL) {
  # Function for generating lower plot
  ggpairs_lower_scatter_basic <- function(data, mapping, ...){
    p <- ggplot(data = data, mapping = mapping) +
      geom_point(color = 'gray10', alpha = 0.3, size = 0.4) +
      geom_abline(color = 'firebrick', linewidth = 0.5)
    p
  }
  # Build and save ggpairs plot
  png(path, height = 2400, width = 2400, units = "px", res = 200)
  control_corplot <-
    ggpairs(data = data, columns = 1:n_controls,
            lower = list(continuous = ggpairs_lower_scatter_basic)) +
    theme_pubr(base_size = 10) +
    theme(strip.text.x = element_text(size = 8),
          strip.text.y = element_text(size = 8),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1))

  print(control_corplot)
  suppressMessages(dev.off())
  return(control_corplot)
}

################################################################################
#* --                                                                      -- *#
#* --                       plot_control_violin()                          -- *#
#* --                                                                      -- *#
################################################################################

plot_control_violin <- function(data = NULL, pal = NULL) {
  p <-
    ggplot(data, aes(x = sample_label, y = counts)) +
      geom_violin(aes(fill = cancer_type), trim = FALSE) +
      scale_y_continuous(breaks = pretty_breaks(12)) +
      scale_fill_manual(values = pal) +
      geom_boxplot(width = 0.1) +
      labs(fill = "Cancer Type", x = '', y = 'Normalised counts') +
      theme_pubr(base_size = 14) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10))
  return(p)
}

################################################################################
#* --                                                                      -- *#
#* --                     plot_control_essentials()                        -- *#
#* --                                                                      -- *#
################################################################################

plot_control_essentials <- function(data = NULL, palette = NULL) {
  p <-
    ggplot(data %>% filter(counts < 2000 & sample_label == 'control_mean'), aes(x = counts)) +
      geom_density(aes(fill = sgrna_group), alpha = 0.3) +
      facet_grid(sgrna_group ~ ., scales = 'free_x') +
      scale_y_continuous(breaks = pretty_breaks(4)) +
      scale_x_continuous(breaks = pretty_breaks(20)) +
      scale_fill_manual(values = palette) +
      labs(fill = "Group", x = '', y = 'Normalised counts < 2000 (control_mean)') +
      theme_pubr(base_size = 14) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10))
  return(p)
}

################################################################################
#* --                                                                      -- *#
#* --             prepare_correlation_data_for_boxplot()                   -- *#
#* --                                                                      -- *#
################################################################################

prepare_correlation_data_for_boxplot <- function(data = NULL) {
  # Filter to remove the control correlation between cell lines and control_mean
  data <- data %>%
    filter(grepl('control', sample_label.x) & grepl('control', sample_label.y) |
             !grepl('control', sample_label.x) & !grepl('control', sample_label.y)) %>%
    filter(sample_label.x != 'control_mean' & sample_label.y != 'control_mean') %>%
    filter(!(grepl('control', sample_label.x) & cell_line_label.x != cell_line_label.y))

  # Add a new category which can be used to determine if samples being compared were within the same cell line
  data <- data %>%
    mutate('correlation_category' = case_when(grepl('control', sample_label.x) & cell_line_label.x == cell_line_label.y  ~ paste(cell_line_label.x, '(control)'),
                                              cell_line_label.x != cell_line_label.y ~ 'Between cell lines',
                                              cell_line_label.x == cell_line_label.y ~ 'Within cell lines'))

  # Add to the bottom the same data frame with the correlation category as the cell line
  # Allows us to see both the per-cell line correlations and the broad summaries within and between cell lines (across all cell lines)
  data_for_cor_boxplot <-
    rbind(data %>% select(sample_label.x, sample_label.y, 'key' = correlation_category, r),
          data %>%
            filter(correlation_category == 'Within cell lines') %>%
            select(sample_label.x, sample_label.y, 'key' = cell_line_label.x, r))

  # Set the cancer_type for summary data and for controls to split these out in plot(s)
  data_for_cor_boxplot <- data_for_cor_boxplot %>%
    left_join(sample_mapping %>% select(cell_line_label, cancer_type), by = c('key' = 'cell_line_label'), relationship = "many-to-many")
  data_for_cor_boxplot$cancer_type <- factor(data_for_cor_boxplot$cancer_type, levels = c(levels(sample_mapping$cancer_type), 'Control', 'Summary'))
  for (i in 1:nrow(data_for_cor_boxplot)) {
    if (is.na(data_for_cor_boxplot$cancer_type[i])) {
      data_for_cor_boxplot$cancer_type[i] <- 'Summary'
    }
    if (grepl('control', data_for_cor_boxplot$sample_label.x[i]) == 1) {
      data_for_cor_boxplot$cancer_type[i] <- 'Control'
    }
  }
  data_for_cor_boxplot$cancer_type <- factor(data_for_cor_boxplot$cancer_type, levels = c(levels(sample_mapping$cancer_type), 'Control', 'Summary'))
  data_for_cor_boxplot$key <- factor(data_for_cor_boxplot$key, levels = unique(data_for_cor_boxplot[order(data_for_cor_boxplot$cancer_type), 'key']))
  data_for_cor_boxplot <- unique(data_for_cor_boxplot)

  # Add in label column
  data_for_cor_boxplot <- data_for_cor_boxplot %>%
    rowwise() %>%
    mutate(label = paste(sep = " : ", sample_label.x, sample_label.y))

  return(data_for_cor_boxplot)
}

################################################################################
#* --                                                                      -- *#
#* --                    plot_sample_correlation()                         -- *#
#* --                                                                      -- *#
################################################################################

plot_sample_correlation <- function(data = NULL) {
  p <-
    ggplot(data, aes(x = key, y = r, fill = cancer_type)) +
      geom_point(color = 'black', alpha = 0.7) +
      geom_boxplot(alpha = 0.7) +
      geom_text_repel(aes(label = label), data = boxplot_cor_data %>% filter(!key %in% c('Within cell lines', 'Between cell lines')) %>% filter(r < 0.55) %>% filter(key == 'WM3702'), nudge_x = 7, size = 3) +
      geom_text_repel(aes(label = label), data = boxplot_cor_data %>% filter(!key %in% c('Within cell lines', 'Between cell lines')) %>% filter(r < 0.55) %>% filter(key == 'COLO 792'), nudge_x = -7, size = 3) +
      geom_text_repel(aes(label = label), data = boxplot_cor_data %>% filter(!key %in% c('Within cell lines', 'Between cell lines')) %>% filter(r < 0.55) %>% filter(key == 'SK-MEL-5'), nudge_x = 7, size = 3) +
      geom_text_repel(aes(label = label), data = boxplot_cor_data %>% filter(!key %in% c('Within cell lines', 'Between cell lines')) %>% filter(r < 0.55) %>% filter(key == 'SK-MEL-28'), nudge_x = -7, size = 3) +
      geom_hline(yintercept = 0.55, linewidth = 0.5, color = 'firebrick', linetype = 'dashed') +
      facet_grid(. ~ cancer_type, scales = 'free_x', space = 'free') +
      scale_fill_manual(values = cancer_type_pal) +
      scale_y_continuous(breaks = pretty_breaks(16), limits = c(0.2, 1)) +
      labs(fill = "Group", x = '', y = "Spearman's correlation coefficient (r)") +
      theme_pubr(base_size = 14) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),
            panel.border = element_blank(),
            strip.background = element_blank(),
            strip.text = element_blank())
  return(p)
}

################################################################################
#* --                                                                      -- *#
#* --                            plot_pca()                                -- *#
#* --                                                                      -- *#
################################################################################

plot_pca <- function(data = NULL, outliers = TRUE) {

  p <-
    ggplot(data, aes(x = PC1, y = PC2, color = cancer_type, label = sample_label)) +
      geom_point(size = 3, alpha = 0.8) +
      geom_text_repel(data = subset(scores, grepl('control', sample_label)), nudge_x = 0.3, size = 5, force = 12, show.legend = F, alpha = 0.8)
  if (outliers == TRUE) {
    p <- p +
      geom_text_repel(data = subset(scores, PC1 > 1 & PC2 < 0), nudge_x = 0.3, size = 5, force = 12, show.legend = F, alpha = 0.8)
  }
  p <- p +
      scale_color_manual(values = cancer_type_pal, name = 'Cancer type') +
      scale_x_continuous(breaks = pretty_breaks(18), limits = c(-1.2, 3.4)) +
      scale_y_continuous(breaks = pretty_breaks(10), limits = c(-1.2, 1)) +
      labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
           y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
      theme_pubr(base_size = 14) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
            axis.text.y = element_text(size = 10))
  return(p)
}

################################################################################
#* --                                                                      -- *#
#* --                  plot_ess_noness_lfc_density()                       -- *#
#* --                                                                      -- *#
################################################################################

plot_ess_noness_lfc_density <- function(data = NULL) {
  p <-
    ggplot(data, aes(x = lfc, fill = sgrna_group)) +
    geom_density(alpha = 0.5) +
    facet_grid(cell_line_label ~ replicate, ) +
    scale_fill_manual(values = sgrna_group_pal) +
    scale_y_continuous(breaks = pretty_breaks(6), limits = c(0, 1.4)) +
    geom_vline(xintercept = 0, linewidth = 0.5, linetype = 'dashed', color = 'gray40', alpha = 0.7) +
    theme(strip.text.y.right = element_text(angle = 0)) +
    labs(fill = "Guide group", x = 'LFC', y = "Density") +
    theme_pubr(base_size = 14) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),
          axis.text.y = element_text(size = 10),
          panel.border = element_blank(),
          strip.background = element_blank())
  return(p)
}

################################################################################
#* --                                                                      -- *#
#* --                           calculate_nnmd()                           -- *#
#* --                                                                      -- *#
################################################################################

calculate_nnmd <- function(data = NULL) {
  results <- data %>%
    group_by(sample_label, sgrna_group) %>%
    summarise(.groups = 'keep',
              'mean' = mean(lfc),
              'median' = median(lfc),
              'sd' = sd(lfc)) %>%
    gather(key, value, -sample_label, -sgrna_group) %>%
    pivot_wider(names_from = c(key,sgrna_group), values_from = c(value)) %>%
    mutate('NNMD' = (`mean_Essential` - `mean_Non-essential`) / `sd_Non-essential`) %>%
    arrange(-NNMD)
  return(results)
}

################################################################################
#* --                                                                      -- *#
#* --                           plot_nnmd()                                -- *#
#* --                                                                      -- *#
################################################################################

plot_nnmd <- function(data = NULL) {
  p <-
    ggplot(data, aes(x = sample_label, y = NNMD, fill = cancer_type)) +
    geom_col(color = 'gray5', linewidth = 0.1, alpha = 0.8) +
    scale_y_reverse(breaks = pretty_breaks(12)) +
    scale_fill_manual(values = cancer_type_pal) +
    facet_grid(. ~ cell_line_label, scales = 'free_x', space = 'free') +
    geom_hline(yintercept = -2, color = 'gray30', linetype = 'dashed', linewidth = 0.5) +
    theme_pubr(base_size = 14) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
          panel.border = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank()) +
    labs(color = "Category",  x = '', y = 'NNMD')
  return(p)
}

################################################################################
#* --                                                                      -- *#
#* --                      plot_scaled_lfc_boxplot()                       -- *#
#* --                                                                      -- *#
################################################################################

plot_scaled_lfc_boxplot <- function(data = NULL, pal = NULL) {
  p <-
    ggplot(data %>% filter(sgrna_group %in% c('Essential', 'Non-essential')), aes(x = cell_line_label, y = LFC, fill = sgrna_group)) +
      geom_boxplot(outlier.size = 0.1) +
      facet_grid(type ~ cancer_type, scales = 'free_x') +
      scale_y_continuous(breaks = pretty_breaks(12)) +
      scale_fill_manual(values = pal) +
      labs(fill = "Cancer Type", x = '', y = 'LFC') +
      theme_pubr(base_size = 14) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
            panel.border = element_blank(),
            strip.background = element_blank())
    return(p)
}

################################################################################
#* --                                                                      -- *#
#* --                plot_scaled_lfc_violin_by_guide_source()              -- *#
#* --                                                                      -- *#
################################################################################

plot_scaled_lfc_violin_by_guide_source <- function(data = NULL) {
  p <- 
    ggplot(data, aes(x = sgrna_group, y = scaled_LFC)) +
      geom_hline(yintercept = 0, alpha = 0.5) +
      geom_violin(aes(fill = sgrna_group)) +
      scale_y_continuous(breaks = pretty_breaks(10)) +
      scale_fill_nejm(alpha = 0.3, name = 'Guide source') +
      geom_boxplot(width = 0.1) +
      labs(y = 'Scaled LFC', x = '') +
      theme_pubr(base_size = 14) +
      theme(axis.text.x = element_text(size = 12))
  return(p)
}

################################################################################
#* --                                                                      -- *#
#* --               plot_scaled_lfc_with_guide_type_distribution()         -- *#
#* --                                                                      -- *#
################################################################################

plot_scaled_lfc_with_guide_type_distribution <- function(data = NULL, annotation_colnames = NULL) {
  data.mean <- data %>%
    mutate('mean_lfc' = rowMeans(.[,(length(annotation_colnames)+1):ncol(data)])) %>%
    select(all_of(annotation_colnames), 'mean_lfc') %>%
    arrange(mean_lfc) %>%
    mutate('known_ess' = ifelse(sgrna_group == 'Essential', 1, 0)) %>%
    mutate('known_noness' = ifelse(sgrna_group == 'Non-essential', 1, 0)) %>%
    mutate('unknown_single' = ifelse(known_ess == 0 & known_noness == 0 & !is.na(singles_target_gene), 1, 0)) %>%
    mutate('unknown_double' = ifelse(known_ess == 0 & known_noness == 0 & unknown_single == 0, 1, 0))

  b <- ggplot(data.mean, aes(x = reorder(id, mean_lfc), y = mean_lfc)) +
        geom_bar(stat = 'identity', fill = 'gray50') +
        geom_hline(yintercept = 0) +
        scale_y_continuous(breaks = pretty_breaks(10), limits = c(-2, 1)) +
        labs(x = '', y = 'Mean scaled LFC') +
        theme_pubr() +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank())

  known_ess <- ggplot(data.mean, aes(x = reorder(id, mean_lfc), y = known_ess, fill = known_ess)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = 'gray40', high = 'firebrick') +
    scale_y_continuous(breaks = c(0,1)) +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA)) +
    ggtitle("Known essential")

  known_noness <- ggplot(data.mean, aes(x = reorder(id, mean_lfc), y = known_noness, fill = known_noness)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = 'gray40', high = 'firebrick') +
    scale_y_continuous(breaks = c(0,1)) +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA)) +
    ggtitle("Known non-essential")

  unknown_single <- ggplot(data.mean, aes(x = reorder(id, mean_lfc), y = unknown_single, fill = unknown_single)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = 'gray40', high = 'navyblue') +
    scale_y_continuous(breaks = c(0,1)) +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA)) +
    ggtitle("Unknown single")

  unknown_dual <- ggplot(data.mean, aes(x = reorder(id, mean_lfc), y = unknown_double, fill = unknown_double)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = 'gray40', high = 'navyblue') +
    scale_y_continuous(breaks = c(0,1)) +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA)) +
    ggtitle("Unknown dual")

  p <- grid.arrange(b, known_ess, known_noness, unknown_single, unknown_dual, nrow = 5, heights = c(7, 0.75, 0.75, 0.75, 0.75))
  return(p)
}