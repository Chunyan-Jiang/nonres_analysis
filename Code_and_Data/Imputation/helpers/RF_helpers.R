# -------------------------------------------------------------------
# Helper functions for simulation-based imputation evaluation.
#
# Purpose:
#   This script defines helper functions for BS and BE Random Forest (RF)
#   simulation-based imputation.
#
# Contents:
#   1) simulate_mask_row():
#        - Simulates artificial missingness by masking observed values,
#          allowing both single-point and consecutive-block missing patterns.
#        - Adds {var}_masked, is_masked, and mask_rel_pos indicators to the data.
#
#   2) impute_masked_rf_rolling():
#        - Implements rolling Random Forest imputation for masked blocks.
#        - Trains a classification RF on non-masked observations and imputes
#          masked values sequentially by relative position within each block.
#        - Stores class probabilities for calibration.
#
#   3) simulate_imputation_evaluation():
#        - Wraps masking and RF imputation into a repeated simulation loop.
#        - Computes evaluation metrics (kappa, accuracy, MAE, rank correlations)
#          exclusively on artificially masked observations.
#        - Stores full imputed datasets with probability outputs for calibration.
#
# Usage:
#   These helper functions are sourced by RF simulation
# AI Usage: Chatgpt is used to debug and add comments
# -------------------------------------------------------------------




################################################################################
#                         SIMULATE MISSINGNESS FUNCTION                        #
################################################################################
simulate_mask_row <- function(df,
                              var = "vg_statebus",
                              id = "idnum",
                              time = "shift_date",
                              block_length = 2,
                              sampling_rate = 0.1) {
  # Simulate missingness by masking observed values.
  # Returns:
  #   - {var}_masked: copy of var with simulated missing values (NA)
  #   - is_masked: TRUE for values masked by this simulation (evaluation subset)
  #   - mask_rel_pos: relative position within a masked block (1..block_length)
  #   - original dataframe
  # Notes:
  #   - block_length controls consecutive missingness length
  
  # Step 0: rank
  df <- df %>%
    arrange(.data[[id]], .data[[time]]) %>%
    group_by(.data[[id]]) %>%
    mutate(row_in_group = row_number()) %>%
    ungroup()
  
  n_total <- sum(!is.na(df[[var]]))
  masked_var <- paste0(var, "_masked")
  
  message("Target masking: ~", round(sampling_rate * 100), "% of data.")
  
  # =====================================================================
  # Case 1: block_length == 1
  # =====================================================================
  if (block_length == 1) {
    candidates <- which(!is.na(df[[var]]))
    n_mask <- round(length(candidates) * sampling_rate)
    mask_idx <- sample(candidates, n_mask)
    
    df <- df %>%
      mutate(
        "{masked_var}" := .data[[var]],
        is_masked      = row_number() %in% mask_idx,
        mask_rel_pos   = ifelse(is_masked, 1L, NA_integer_),
        "{masked_var}" := ifelse(is_masked, NA, .data[[var]])
      )
    
    message("Masked ", n_mask, " records (",
            round(100 * n_mask / n_total, 2), "% of total).")
    return(df)
  }
  
  # =====================================================================
  # Case 2: block_length >= 2
  # =====================================================================
  
  # 1. Identify valid block start positions per ID
  candidates <- df %>%
    group_by(.data[[id]]) %>%
    mutate(
      # can_start = TRUE if the next block_length rows are all non-missing
      can_start = zoo::rollapplyr(
        !is.na(.data[[var]]),
        width   = block_length,
        FUN     = all,
        align   = "left",
        partial = FALSE,
        fill    = FALSE
      )
    ) %>%
    ungroup() %>%
    filter(can_start) %>%
    select(all_of(id), row_in_group)
  
  message("Found ", nrow(candidates), " valid block start positions.")
  
  # Number of blocks to sample
  n_blocks <- round(n_total * sampling_rate / block_length)
  n_sample <- min(nrow(candidates), n_blocks)
  
  if (n_sample == 0) {
    warning("No valid block start positions found.")
    df <- df %>%
      mutate("{masked_var}" := .data[[var]], is_masked = FALSE)
    return(df)
  }
  
  # 2. Randomly sample block start points (no overlap check)
  # Important: keep sample order (affects overwrite priority)
  sampled <- candidates[sample(nrow(candidates), n_sample), ] %>%
    mutate(
      start         = row_in_group,
      end           = row_in_group + block_length - 1,
      mask_block_id = row_number()         # later blocks overwrite earlier ones
    )
  
  message("Using ", n_sample, " blocks for masking.")
  
  # 3. Expand blocks into row ranges and assign relative position
  block_rows <- sampled %>%
    rowwise() %>%
    mutate(row_seq = list(start:end)) %>%
    unnest(row_seq) %>%
    mutate(
      row_in_group = row_seq,
      mask_rel_pos = row_seq - start + 1
    ) %>%
    select(all_of(id), row_in_group, mask_block_id, mask_rel_pos)
  
  # 4. If a row is covered by multiple blocks, keep the one with largest mask_block_id
  block_rows <- block_rows %>%
    group_by(!!sym(id), row_in_group) %>%
    slice_max(mask_block_id, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(-mask_block_id)
  
  # 5. Merge into df and generate masked output
  df <- df %>%
    left_join(block_rows, by = c(id, "row_in_group")) %>%
    mutate(
      is_masked     = !is.na(mask_rel_pos),
      "{masked_var}" := ifelse(is_masked, NA, .data[[var]])
    )
  
  n_masked <- sum(df$is_masked, na.rm = TRUE)
  message("Masked ", n_masked, " records (",
          round(100 * n_masked / n_total, 2), "% of total).")
  
  return(df)
}

################################################################################
#                         IMPUTATION FUNCTION                                  #
################################################################################
impute_masked_rf_rolling <- function(df,
                                     var = "vg_statebus",
                                     masked_var = "vg_statebus_masked",
                                     mask_flag = "is_masked",
                                     relpos_var = "mask_rel_pos",
                                     block_length = 2,
                                     id_var = "idnum",
                                     lag_list = c('vg_statebus_lag1',
                                                  'vg_statebus_lag2',
                                                  'vg_statebus_lag3',
                                                  'vg_statebus_lag4'),
                                     covariates = c('idnum','vg_westeast','sector_total','online',
                                                    'vg_statebus_lag1','vg_statebus_lag2','vg_statebus_lag3','vg_statebus_lag4',
                                                    'vg_comexp_lag1','vg_comexp_lag2','vg_comexp_lag3','vg_comexp_lag4',
                                                    'calendar_time','is_dec','is_aug'),
                                     ...) {
  
  # Rolling RF imputation for masked blocks:
  # - Train a classification RF model on non-masked rows (observed targets)
  # - Impute masked rows sequentially by relative position within the block (mask_rel_pos = 1..L)
  # - After imputing position p, update lagged predictors within each id to support imputing p+1
  #
  # Output:
  # - {var}_imputed: imputed target values (character/factor levels from RF classes)
  # - prob: list-column storing per-class predicted probabilities for imputed rows
  
  
  message("Start Rolling RF Imputation for: ", var)
  
  imputed_var <- paste0(var, "_imputed")
  
  # initialize imputed and prob list-column
  df[[imputed_var]] <- df[[masked_var]]
  df$prob <- vector("list", nrow(df))
  
  df_roll <- df
  df_roll[[imputed_var]] <- df_roll[[masked_var]]
  df_roll$prob <- vector("list", nrow(df_roll))
  
  #--------------------------------------------------
  # 1. Train RF (now explicitly probability=TRUE)
  #--------------------------------------------------
  train_df <- df %>%
    dplyr::filter(!.data[[mask_flag]]) %>%
    dplyr::select(dplyr::all_of(c(var, covariates))) %>%
    tidyr::drop_na()
  
  formula_rf <- stats::as.formula(paste(var, "~", paste(covariates, collapse = "+")))
  
  rf_model <- ranger::ranger(
    formula = formula_rf,
    data    = train_df,
    # Train RF with probability=TRUE to store class probabilities (used for calibration).
    probability = TRUE,
    ...
  )
  
  # impute masked blocks from position 1 to block_length, updating lags after each step.
  #--------------------------------------------------
  # 2. Rolling imputation loop:   
  #--------------------------------------------------
  for (p in 1:block_length) {
    
    message("  Imputing for relative position: ", p)
    
    rows_p <- which(df[[relpos_var]] == p)
    if (length(rows_p) == 0) next
    
    test_df <- df_roll[rows_p, covariates, drop = FALSE]
    
    # get probabilities instead of only predictions
    pred_obj  <- predict(rf_model, data = test_df, type = "response")
    pred_probs <- pred_obj$predictions            # matrix (n x K)
    
    # mode class remains your imputed value
    pred_class <- colnames(pred_probs)[max.col(pred_probs)]
    
    # Write prediction class back (unchanged original logic)
    df_roll[[imputed_var]][rows_p] <- pred_class
    
    # write probabilities to list-column
    df_roll$prob[rows_p] <- split(as.data.frame(pred_probs), seq_len(nrow(pred_probs)))
    
    # Update lag variables after each imputation step to maintain temporal consistency within each id.
    df_roll <- df_roll %>%
      dplyr::group_by(.data[[id_var]])
    
    for (k in seq_along(lag_list)) {
      lag_name <- lag_list[k]
      df_roll <- df_roll %>%
        dplyr::mutate(
          !!rlang::sym(lag_name) := dplyr::lag(.data[[imputed_var]], k)
        )
    }
    
    df_roll <- df_roll %>% dplyr::ungroup()
    
    for (lag_name in lag_list) {
      na_idx <- is.na(df_roll[[lag_name]])
      if (any(na_idx)) {
        df_roll[[lag_name]][na_idx] <- df[[lag_name]][na_idx]
      }
    }
  }
  
  #--------------------------------------------------
  # 3. Copy results back
  #--------------------------------------------------
  df[[imputed_var]] <- df_roll[[imputed_var]]
  df$prob <- df_roll$prob    # added
  
  message("Rolling RF imputation completed. Output variable: ", imputed_var)
  
  return(df)
}

################################################################################
#                         SIMULATION FUNCTION                                  #
################################################################################
simulate_imputation_evaluation <- function(
    df,
    var = "vg_statebus",
    id = "idnum",
    time = "shift_date",
    mask_fun = simulate_mask_row,
    mask_args = list(block_length = 2, sampling_rate = 0.1),
    impute_fun = impute_masked_rf_rolling,
    impute_args = list(masked_var = "vg_statebus_masked",
                       mask_flag = "is_masked"),
    n_sim = 30,
    seed = 123
) {
  
  results_all <- list()
  stored_dfs <- list()
  
  for (s in seq_len(n_sim)) {
    message("\n===== Simulation ", s, " / ", n_sim, " =====")
    # Reproducible randomness:
    # Use a fixed base seed and vary by simulation index (seed + s) so each run is distinct but reproducible.
    set.seed(seed + s)
    
    # Step 1: Masking
    masked_data <- do.call(
      mask_fun,
      c(list(df = df, var = var, id = id, time = time), mask_args)
    )
    
    # Step 2: Imputation (with prob)
    imputed_data <- do.call(
      impute_fun,
      c(list(df = masked_data, var = var), impute_args)
    )
    
    # store full dataframe with probabilities
    stored_dfs[[s]] <- imputed_data %>% dplyr::mutate(sim = s)
    
    # Step 3: evaluation
    # Evaluation is computed ONLY on rows that were artificially masked (is_masked == TRUE).
    imputed_col <- paste0(var, "_imputed")
    mask_flag <- impute_args$mask_flag %||% "is_masked"
    
    compare_df <- imputed_data %>%
      dplyr::filter(.data[[mask_flag]]) %>%
      dplyr::select(all_of(var), all_of(imputed_col)) %>%
      dplyr::filter(complete.cases(across(everything())))
    
    t <- compare_df[[var]]
    p <- compare_df[[imputed_col]]
    
    metrics_summary <- tibble::tibble(
      sim = s,
      # Variable names prefixed with 'mean_' represent pooled metrics for the current simulation iteration. 
      mean_kappa    = as.numeric(psych::cohen.kappa(table(t, p))$kappa),
      mean_accuracy = mean(t == p),
      mean_mae      = mean(abs(as.numeric(t) - as.numeric(p))),
      mean_rho      = suppressWarnings(cor(as.numeric(t), as.numeric(p), method = "spearman")),
      mean_tau      = suppressWarnings(cor(as.numeric(t), as.numeric(p), method = "kendall"))
    )
    
    results_all[[s]] <- metrics_summary
  }
  
  # Step 4: summary
  results_df <- dplyr::bind_rows(results_all)
  
  results_summary <- results_df %>%
    dplyr::summarise(
      avg_kappa    = mean(mean_kappa, na.rm = TRUE),
      avg_accuracy = mean(mean_accuracy, na.rm = TRUE),
      avg_mae      = mean(mean_mae, na.rm = TRUE),
      avg_rho      = mean(mean_rho, na.rm = TRUE),
      avg_tau      = mean(mean_tau, na.rm = TRUE)
    )
  
  list(
    all_runs = results_df,          # metrics of each simulation
    summary = results_summary,      # averages
    imputed_data_list = stored_dfs  # full df of each simulation (with prob)
  )
}