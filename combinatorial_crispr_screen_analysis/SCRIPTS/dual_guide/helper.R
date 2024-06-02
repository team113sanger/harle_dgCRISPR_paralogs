################################################################################
#* --                                                                      -- *#
#* --                     get_pred_vs_obs_y12()                            -- *#
#* --                                                                      -- *#
################################################################################

get_pred_vs_obs_y12 <- function(sn = NULL, gm = NULL, fcs = NULL) {
  # Create empty data frame for pred vs obs values
  pred_vs_obs_y12 <- data.frame(id = character(), sorted_gene_pair = character(), sample = character(),
                                obs_y12 = numeric(), pred_y12 = numeric(),
                                y1 = numeric(), y2 = numeric())

  # Create empty data frame for missing guides per sample
  missing_data <- data.frame(id = character(), sorted_gene_pair = character(), sample = character())

  for (i in 1:nrow(gm)) {
    # Get dual and respective single guide ids
    g12 <- as.vector(gm$id[i])
    g1  <- as.vector(gm$g1[i])
    g2  <- as.vector(gm$g2[i])
    gp  <- as.vector(gm$sorted_gene_pair[i])

    # Get observed fold changes for left guide (single)
    y1 <- fcs %>%
      filter(id == g1 & sorted_gene_pair == gp & sample == sn) %>%
      select(fc) %>%
      unlist() %>% as.vector()

    # Get observed fold changes for right guide (single)
    y2 <- fcs %>%
      filter(id == g2 & sorted_gene_pair == gp & sample == sn) %>%
      select(fc) %>%
      unlist() %>% as.vector()

    # Get observed fold changes for dual guide (dual)
    obs_y12 <- fcs %>%
      filter(id == g12 & sorted_gene_pair == gp & sample == sn) %>%
      select(fc) %>%
      unlist() %>% as.vector()

    # Check there are fold changes for the dual guide and both related single guides
    if ( length( y1 ) == 1 && length( y2 ) == 1 && length( obs_y12 ) == 1 ) {
      pred_y12 <- y1 + y2
      pred_vs_obs_y12.tmp <- data.frame(  'id' = g12, 'sorted_gene_pair' = gp, 'sample' = sn,
                                          'obs_y12' = obs_y12, 'pred_y12' = pred_y12,
                                          'y1' = y1, 'y2' = y2 )
      pred_vs_obs_y12 <- rbind(pred_vs_obs_y12, pred_vs_obs_y12.tmp)
    } else {
      # If one of the guides isn't present (filtered), then add to missing data
      missing_data.tmp <- data.frame('id' = g12, 'sorted_gene_pair' = gp, 'sample' = sn)
      missing_data <- rbind(missing_data, missing_data.tmp)
    }
    if(i %% 200 == 0) { print(i) }
  }
  return(list("pred_vs_obs_y12" = pred_vs_obs_y12, "missing_data" = missing_data))
}

################################################################################
#* --                                                                      -- *#
#* --               get_pred_vs_obs_y12_for_all_samples()                  -- *#
#* --                                                                      -- *#
################################################################################

get_pred_vs_obs_y12_for_all_samples <- function(samples, gm, fcs) {
  pred_vs_obs_y12 <- list()
  missing_data    <- list()
  sample.results  <- list()

  for (sn in samples) {
    print(paste("Processing:", sn))
    sample.results <- get_pred_vs_obs_y12(sn, gm, fcs)
    pred_vs_obs_y12[[sn]] <- sample.results[['pred_vs_obs_y12']]
    missing_data[[sn]] <- sample.results[['missing_data']]
  }
  return(list('pred_vs_obs_y12' = pred_vs_obs_y12, 'missing_data' = missing_data))
}
