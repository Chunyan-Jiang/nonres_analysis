
# ------------------------------------------------------------
# Save GAM diagnostics and cluster-robust inference results
#
# This function performs a comprehensive diagnostic workflow
# for a fitted GAM:
#   1. Saves the GAM summary.
#   2. Saves the output of gam.check().
#   3. Saves concurvity diagnostics.
#   4. Computes a cluster-robust covariance matrix using vcovCL.
#   5. Saves cluster-robust coefficient tests.
#
#
# Args:
#   model       : A fitted mgcv::gam object.
#   data        : Data frame used for model fitting.
#   cluster_var : Name of the clustering variable (character).
#   outdir      : Output directory for diagnostic files.
#                 Defaults to the current working directory.
#
# Returns:
#   Writes the following outputs to disk:
#     - Standard GAM summary.
#     - gam.check() diagnostics.
#     - Concurvity diagnostics.
#     - Cluster-robust coefficient test table.
#     - Cluster-robust covariance matrix (RDS file).
#
#   Invisibly returns a list containing:
#     - vcovCL    : Cluster-robust covariance matrix.
#     - coeftest : Cluster-robust coefficient test table.

# ------------------------------------------------------------
save_gam_full_diagnostics <- function(
    model,
    data,
    cluster_var,
    outdir = "."
) {
  
  # Create output directory if it does not exist
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  # Use the model object name for consistent file naming
  model_name <- deparse(substitute(model))
  
  # Helper function to safely capture console output
  write_block <- function(expr, filename) {
    txt <- capture.output(expr)
    writeLines(txt, file.path(outdir, filename))
  }
  
  # Save GAM summary
  write_block(
    summary(model),
    paste0(model_name, "_original_res.txt")
  )
  
  # Save gam.check diagnostics
  write_block(
    gam.check(model),
    paste0(model_name, "_gamcheck.txt")
  )
  
  # Save concurvity diagnostics
  write_block(
    concurvity(model, full = TRUE),
    paste0(model_name, "_concurvity.txt")
  )
  
  # Compute cluster-robust covariance matrix at the cluster level
  V_star <- vcovCL(
    model,
    cluster = data[[cluster_var]]
  )
  
  # Save covariance matrix
  saveRDS(
    V_star,
    file = file.path(outdir, paste0(model_name, "_vcovCL.rds"))
  )
  
  # Compute cluster-robust coefficient tests
  ct_robust <- coeftest(model, vcov. = V_star)
  
  # Save cluster-robust coefficient table
  write_block(
    ct_robust,
    paste0(model_name, "_cr.txt")
  )
  
  # Return objects invisibly for optional downstream use
  invisible(list(
    vcovCL    = V_star,
    coeftest = ct_robust
  ))
}


# ------------------------------------------------------------
# Compute cluster-robust confidence intervals for a GAM smooth
#
# This function extracts a single smooth term from a fitted
# GAM and computes pointwise confidence intervals using a
# cluster-robust covariance matrix.
#
# The implementation is based on the linear predictor matrix
# ("lpmatrix") representation of GAMs.
#
# Args:
#   model        : A fitted mgcv::gam object.
#   newdata      : Data frame at which the smooth is evaluated.
#   smooth_label : Label of the smooth term (e.g. "s(x)").
#   V            : Cluster-robust covariance matrix.
#   level        : Confidence level (default: 0.95).
#
# Returns:
#   A data frame combining newdata with:
#     - fit : Estimated smooth contribution.
#     - se  : Cluster-robust standard error.
#     - lwr : Lower confidence bound.
#     - upr : Upper confidence bound.
# ------------------------------------------------------------
smooth_ci <- function(
    model,
    newdata,
    smooth_label,
    V,
    level = 0.95
) {
  
  # Linear predictor matrix for the new data
  Xp <- predict(model, newdata = newdata, type = "lpmatrix")
  
  # Model coefficient vector
  beta <- coef(model)
  
  # Identify the requested smooth term
  sm_idx <- which(
    vapply(model$smooth, function(s) s$label, "") == smooth_label
  )
  
  if (length(sm_idx) != 1)
    stop("smooth_label not found or not unique.")
  
  sm   <- model$smooth[[sm_idx]]
  cols <- sm$first.para : sm$last.para
  
  # Extract smooth-specific design matrix and parameters
  Xs <- Xp[, cols, drop = FALSE]
  b  <- beta[cols]
  Vs <- V[cols, cols, drop = FALSE]
  
  # Smooth estimate
  fit <- as.numeric(Xs %*% b)
  
  # Cluster-robust standard error
  se <- sqrt(rowSums((Xs %*% Vs) * Xs))
  
  # Normal-approximation critical value
  crit <- qnorm(1 - (1 - level) / 2)
  
  # Return smooth estimate with confidence intervals
  cbind(
    newdata,
    fit = fit,
    se  = se,
    lwr = fit - crit * se,
    upr = fit + crit * se
  )
}
