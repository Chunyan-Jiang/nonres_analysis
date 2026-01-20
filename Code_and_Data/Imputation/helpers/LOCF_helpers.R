# -------------------------------------------------------------------
# Helper functions for simulation-based imputation evaluation.
#
# Purpose:
#   This script defines helper functions for BS and BE LOCF simulation
#
# Contents:
#   1) simulate_mask_row():
#        - Simulates artificial missingness by masking observed values,
#          allowing both single-point and consecutive-block missing patterns.
#        - Adds {var}_masked and is_masked indicators to the data.
#
#   2) impute_masked_basic():
#        - Implements simple baseline imputation methods (LOCF, NOCB, Mode)
#          applied within each id in time order.
#
#   3) simulate_imputation_evaluation():
#        - Wraps the masking and imputation steps into a repeated simulation loop.
#        - Computes evaluation metrics (kappa, accuracy, MAE, rank correlations)
#          exclusively on artificially masked observations.
#
# Usage: These helper functions are sourced by LOCF Simulation
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
  message("Target masking: ~", round(sampling_rate * 100), "% of data.")
  
  # =========================================================
  # Case 1: block_length == 1 
  # =========================================================
  if (block_length == 1) {
    candidates <- which(!is.na(df[[var]]))
    n_mask <- round(length(candidates) * sampling_rate)
    mask_idx <- sample(candidates, n_mask)
    
    df <- df %>%
      mutate(
        "{var}_masked" := .data[[var]],
        is_masked = row_number() %in% mask_idx,
        "{var}_masked" := ifelse(is_masked, NA, .data[[var]])
      )
    
    message("Masked ", n_mask, " records (",
            round(100 * n_mask / n_total, 2), "% of total).")
    return(df)
  }
  
  # =========================================================
  # Case 2: block_length >= 2
  # =========================================================
  # Step 1: find all possible start point
  candidates <- df %>%
    group_by(.data[[id]]) %>%
    mutate(
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
  
  n_blocks <- round(n_total * sampling_rate / block_length)
  message("Found ", nrow(candidates), " valid block start positions.")
  
  if (nrow(candidates) == 0) {
    warning("No valid block start positions found.")
    df <- df %>%
      mutate("{var}_masked" := .data[[var]], is_masked = FALSE)
    return(df)
  }
  
  # Step 2: random sampling
  n_sample <- min(nrow(candidates), n_blocks)
  candidates_sampled <- candidates[sample(nrow(candidates), n_sample), ]
  
  kept <- candidates_sampled %>%
    mutate(start=row_in_group,
           end=row_in_group+block_length-1)
  
  # Step 3: generate block mask
  sampled <- kept %>%
    rename("{id}":=!!sym(id))
  
  message("Using ", nrow(sampled), " blocks for masking.")
  
  df <- df %>%
    mutate("{var}_masked" := .data[[var]], is_masked = FALSE)
  
  mask_df <- sampled %>%
    group_by(.data[[id]]) %>%
    summarise(rows_to_mask = list(unlist(mapply(`:`, start, end))), .groups = "drop")
  
  df <- df %>%
    left_join(mask_df, by = id) %>%
    mutate(
      rows_to_mask = lapply(rows_to_mask, function(x)
        if (is.null(x) || all(is.na(x))) integer(0) else x),
      is_masked = Map(function(i, rs) i %in% rs, row_in_group, rows_to_mask) |> unlist(),
      "{var}_masked" := ifelse(is_masked, NA, .data[[var]])
    ) %>%
    select(-rows_to_mask)
  
  # Step 4: summary
  n_masked <- sum(df$is_masked, na.rm = TRUE)
  message("Masked ", n_masked, " records (",
          round(100 * n_masked / n_total, 2), "% of total).")
  
  return(df)
}



################################################################################
#                         IMPUTATION FUNCTION                                  #
################################################################################

impute_masked_basic <- function(df,
                                var = "vg_statebus",
                                masked_var = "vg_statebus_masked",
                                mask_flag = "is_masked",
                                id = "idnum",
                                time = "shift_date",
                                method = c("locf", "nocb", "mode")) {
  # Basic imputation methods applied within each id (time-ordered):
  # - locf: carry last observed value forward
  # - nocb: carry next observed value backward
  # - mode: fill missing with the global mode within the id (based on available values)
  # Output column:
  # - {var}_imputed
  method <- match.arg(method)
  
  var_imputed <- paste0(var, "_imputed")
  
  df <- df %>%
    arrange(.data[[id]], .data[[time]]) %>%
    group_by(.data[[id]]) %>%
    mutate(
      "{var_imputed}" := case_when(
        # LOCF
        method == "locf" ~ na.locf(.data[[masked_var]], na.rm = FALSE),
        
        # NOCB
        method == "nocb" ~ na.locf(.data[[masked_var]], fromLast = TRUE, na.rm = FALSE),
        
        # MODE
        method == "mode" ~ {
          v <- .data[[masked_var]]
          mode_val <- as.numeric(names(which.max(table(v[!is.na(v)]))))
          ifelse(is.na(v), mode_val, v)
        }
      )
    ) %>%
    ungroup()
  
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
    # Mask param ↓
    mask_fun = simulate_mask_row,
    mask_args = list(block_length = 2, sampling_rate = 0.1),
    # Imputation param ↓
    impute_fun = impute_masked_basic, 
    impute_args = list(masked_var = "vg_statebus_masked",
                       mask_flag = "is_masked",
                       method='locf'),
    # Simulation ↓
    n_sim = 30,
    seed = 123
) {
  
  results_all <- list()
  
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
    
    # Step 2: Imputation
    imputed_data <- do.call(
      impute_fun,
      c(list(df = masked_data, var = var, id = id, time = time), impute_args)
    )
    
    # Step 3: calculate metrics
    # Evaluation is computed ONLY on rows that were artificially masked (is_masked == TRUE).
    imputed_col <- paste0(var, "_imputed")
    mask_flag <- impute_args$mask_flag %||% "is_masked"
    
    compare_df <- imputed_data %>%
      filter(.data[[mask_flag]]) %>%
      select(all_of(var), all_of(imputed_col)) %>%
      filter(complete.cases(across(everything())))
    
    t <- compare_df[[var]]
    p <- compare_df[[imputed_col]]
    
    metrics_summary <- tibble(
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
  results_df <- bind_rows(results_all)
  results_summary <- results_df %>%
    summarise(
      avg_kappa    = mean(mean_kappa, na.rm = TRUE),
      avg_accuracy = mean(mean_accuracy, na.rm = TRUE),
      avg_mae      = mean(mean_mae, na.rm = TRUE),
      avg_rho      = mean(mean_rho, na.rm = TRUE),
      avg_tau      = mean(mean_tau, na.rm = TRUE)
    )
  
  list(
    all_runs = results_df,
    summary = results_summary
  )
}