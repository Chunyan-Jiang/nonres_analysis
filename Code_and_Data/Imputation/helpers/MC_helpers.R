# -------------------------------------------------------------------
# Helper functions for simulation-based imputation evaluation.
#
# Purpose:
#   This script defines helper functions for BS and BE Markov Chain (MC)
#   simulation-based imputation.
#
# Contents:
#   1) simulate_mask_row():
#        - Simulates artificial missingness by masking observed values,
#          allowing both single-point and consecutive-block missing patterns.
#        - Adds {var}_masked and is_masked indicators to the data.
#
#   2) impute_masked_mc():
#        - Implements homogeneous Markov Chain imputation applied within each id
#          in time order.
#        - Builds conditional distributions based on lag history of length h
#          and performs rolling imputation for masked observations.
#
#   3) simulate_imputation_evaluation():
#        - Wraps masking and MC imputation into a repeated simulation loop.
#        - Computes evaluation metrics (kappa, accuracy, MAE, rank correlations)
#          exclusively on artificially masked observations.
#
# Usage:
#   These helper functions are sourced by MC simulation
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
impute_masked_mc <- function(df,
                             var = 'vg_statebus',
                             id = 'idnum',
                             time = 'shift_date',
                             mask_flag = 'is_masked',
                             h = 1,
                             lag_list = c('vg_statebus_lag1',
                                          'vg_statebus_lag2',
                                          'vg_statebus_lag3',
                                          'vg_statebus_lag4'),
                             strategy = c("random", "mode")) {
  # Markov Chain (MC) imputation applied within each id (time-ordered):
  # - Builds a conditional probability table based on history length 'h' (lags).
  # - Dynamically updates future lag columns during iteration (rolling imputation).
  # Strategies:
  # - random: sample value based on conditional probabilities (stochastic).
  # - mode: choose the most probable value given the history (deterministic).
  # Output column:
  # - {var}_imputed
  
  # 1. Setup and validation -------------------------------------------------
  library(data.table)
  strategy <- match.arg(strategy)
  if (!strategy %in% c("random", "mode")) {
    stop("strategy must be either 'random' or 'mode'")
  }
  
  DT <- as.data.table(df)
  
  # Ensure chronological order by id and time
  setorderv(DT, c(id, time))
  
  # 2. Build conditional probability table ----------------------------------
  lag_use <- lag_list[1:h]
  cond_cols <- c(lag_use, var)
  
  # Frequency of each lag pattern
  cond_probs <- DT[!is.na(get(var)), .N, by = cond_cols]
  cond_probs[, p := N / sum(N), by = lag_use]
  
  # Create a pattern key for lag combinations
  cond_probs[, pattern_key := do.call(paste, c(.SD, sep = "_")), .SDcols = lag_use]
  cond_probs_map <- cond_probs[, .(pattern_key, value = get(var), p)]
  setkey(cond_probs_map, pattern_key)
  
  # 3. Global fallback distribution -----------------------------------------
  global_probs <- DT[!is.na(get(var)), .N, by = var]
  global_probs[, p := N / sum(N)]
  
  # 4. Select only IDs with masked rows -------------------------------------
  masked_ids <- unique(DT[get(mask_flag) == 1, get(id)])
  DT_subset <- DT[get(id) %in% masked_ids]
  
  # 5. Rolling imputation (by id) -------------------------------------------
  imputed_result <- DT_subset[, {
    x <- get(var)
    is_mask <- get(mask_flag)
    
    # Iterate within this id (time-ordered)
    for (i in seq_len(.N)) {
      if (is_mask[i] == 1) {  # only impute masked rows
        # Build the lag pattern for the current position
        lag_vals <- if (i > h) sapply(1:h, function(k) x[i - k]) else rep(NA, h)
        pattern_key <- if (all(!is.na(lag_vals))) paste(lag_vals, collapse = "_") else NA_character_
        prob <- if (!is.na(pattern_key)) cond_probs_map[pattern_key] else data.table()
        
        # Impute value using conditional or fallback distribution
        if (i <= h || any(is.na(lag_vals)) || nrow(prob) == 0 || length(prob$value)!=length(prob$p) || anyNA(prob$p)) {
          # Fallback to global distribution
          if (strategy == "random") {
            x[i] <- sample(global_probs[[var]], 1, prob = global_probs$p)
          } else {
            x[i] <- global_probs[[var]][which.max(global_probs$p)]
          }
        } else {
          # Use conditional distribution
          if (strategy == "random") {
            x[i] <- sample(prob$value, 1, prob = prob$p)
          } else {
            x[i] <- prob$value[which.max(prob$p)]
          }
        }
        
        # Update future lag columns within the same id
        for (lag_idx in seq_len(h)) {
          future_row <- i + lag_idx
          if (future_row <= .N) {  # do NOT go beyond current id
            lag_col <- lag_use[lag_idx]
            DT_subset[[lag_col]][future_row] <- x[i]
          } else {
            break
          }
        }
      }
    }
    
    # Return results for this id
    out <- data.table(get(id), get(time), x)
    setnames(out, c(id, time, "var_imputed"))
    out
  }, by = id]
  
  # 6. Merge imputed values back to the original data -----------------------
  setkeyv(DT, c(id, time))
  setkeyv(imputed_result, c(id, time))
  DT[imputed_result, var_imputed := i.var_imputed]
  
  # 7. Convert back to data.frame and rename output column ------------------
  df_out <- as.data.frame(DT)
  new_colname <- paste0(var, "_imputed")
  names(df_out)[names(df_out) == "var_imputed"] <- new_colname
  
  # 8. Return final result --------------------------------------------------
  return(df_out)
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
      mean_kappa    = as.numeric(psych::cohen.kappa(cbind(as.character(t), as.character(p)))$kappa),
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