#!/usr/bin/env Rscript
# ==============================================================================
# 02_kaplan_meier.R - Kaplan-Meier survival analysis
# ==============================================================================
# Purpose:
#   - Load fitted models from 01_fit_models.R
#   - Perform log-rank tests for PRS quantiles
#   - Calculate survival percentiles (median, 25th, 75th, custom)
#   - Calculate RMST (Restricted Mean Survival Time)
#   - Log-rank trend test across ordered quantiles
#   - Pairwise log-rank comparisons with multiple testing correction
#   - Incidence rates (events per 1000 person-years)
#   - Create risk tables
#   - Stratified analyses (by sex, age groups, etc.)
#
# Inputs:
#   - fitted_models.rds (from 01_fit_models.R)
#   - data_processed.rds (from 01_fit_models.R)
#
# Outputs:
#   - prefix_logrank_tests.csv
#   - prefix_percentile_survival.csv (customizable percentiles)
#   - prefix_rmst.csv
#   - prefix_logrank_trend.csv
#   - prefix_pairwise_comparisons.csv
#   - prefix_incidence_rates.csv
#   - prefix_risk_table.csv
#   - prefix_stratified_results.csv (if stratification used)
#
# Author: Ying Wang
# ==============================================================================

# ==============================================================================
# PACKAGE MANAGEMENT
# ==============================================================================

# Dependency policy: CHECK, do not install. Silent runtime installation fails in
# biobanks without internet access and lets sites drift to different package
# versions, which is exactly what a multi-site consortium must avoid. Set
# PRS_ALLOW_INSTALL=1 to opt into installation for interactive/dev use.
.PRS_MIN_VERSIONS <- c(optparse="1.7.0", data.table="1.14.0", survival="3.5",
                       riskRegression="2023.03.22", ggplot2="3.4.0", survminer="0.4.9",
                       rms="6.0", pec="2022.05.04", boot="1.3", scales="1.2.0",
                       RColorBrewer="1.1", rlang="1.0.0", gridExtra="2.3", ggpubr="0.6.0",
                       timeROC="0.4", survAUC="1.1", prodlim="2019.11.13", splines="0")
load_package <- function(pkg) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  if (!ok) {
    if (identical(Sys.getenv("PRS_ALLOW_INSTALL"), "1")) {
      cat("Installing missing package (PRS_ALLOW_INSTALL=1):", pkg, "\n")
      install.packages(pkg, repos = "https://cloud.r-project.org/", quiet = TRUE)
    } else {
      stop("Required package not installed: ", pkg,
           "\n  Install it, then re-run:\n    install.packages(\"", pkg, "\")",
           "\n  (or set PRS_ALLOW_INSTALL=1 to allow automatic installation)")
    }
  }
  minv <- .PRS_MIN_VERSIONS[[pkg]]
  if (!is.null(minv) && nzchar(minv) && minv != "0") {
    have <- tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) NA_character_)
    if (!is.na(have) && utils::compareVersion(have, minv) < 0)
      stop("Package '", pkg, "' is version ", have, " but >= ", minv,
           " is required for reproducible consortium results.",
           "\n  Update it:  install.packages(\"", pkg, "\")")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

required_packages <- c("optparse", "data.table", "survival", "splines", "ggplot2", "survminer",
                       "scales", "RColorBrewer")  # Cairo omitted (optional; grDevices fallback in save_km)
invisible(lapply(required_packages, load_package))

# Disclosure helper (identical to prs_risk_utils::suppress_small_cells; inlined
# because this script does not source that file). Suppresses a row when ANY count
# column is a small non-zero integer (1..threshold-1), NA-ing counts + estimates.
suppress_small_cells <- function(x, count_cols, threshold = 10L, estimate_cols = character()) {
  count_cols    <- intersect(count_cols, names(x))
  estimate_cols <- intersect(estimate_cols, names(x))
  if (!length(count_cols)) { x[, minimum_cell_count_pass := TRUE]; return(x[]) }
  small <- Reduce(`|`, lapply(count_cols, function(cc) {
    v <- suppressWarnings(as.numeric(x[[cc]])); is.finite(v) & v > 0 & v < threshold
  }))
  small[is.na(small)] <- FALSE
  x[, minimum_cell_count_pass := !small]
  if (any(small)) for (cc in unique(c(count_cols, estimate_cols))) x[small, (cc) := NA]
  x[]
}

# ==============================================================================
# COMMAND-LINE ARGUMENTS
# ==============================================================================

option_list <- list(
  # Input files
  make_option("--models_file", type="character", default=NULL,
              help="Path to fitted_models.rds from 01_fit_models.R [REQUIRED]"),
  make_option("--data_file", type="character", default=NULL,
              help="Path to data_processed.rds from 01_fit_models.R [REQUIRED]"),
  
  # Analysis options
  make_option("--quantiles", type="character", default="all",
              help="Which quantiles to analyze: 'all', 'Q5', 'Q10', 'extreme_1', etc. Comma-separated [default: all]"),
  make_option("--risk_table_times", type="character", default=NULL,
              help="Times for risk table (comma-separated, e.g., '5,10,15') [default: auto]"),
  
  # NEW: Percentile survival options
  make_option("--survival_percentiles", type="character", default="25,50,75",
              help="Percentiles for survival time estimation (comma-separated, e.g., '10,25,50,75,90') [default: 25,50,75]"),
  
  # NEW: RMST options
  make_option("--calculate_rmst", action="store_true", default=FALSE,
              help="Calculate Restricted Mean Survival Time [default: FALSE]"),
  make_option("--rmst_tau", type="character", default=NULL,
              help="Comma-separated time horizon(s) for RMST (e.g., '10,15,20'). If NULL, uses max follow-up [default: NULL]"),
  make_option("--rmst_bootstrap", action="store_true", default=FALSE,
              help="Use bootstrap to estimate RMST standard errors (slower) [default: FALSE]"),
  make_option("--rmst_n_boot", type="integer", default=200,
              help="Number of bootstrap samples for RMST SE [default: 200]"),
  
  # NEW: Trend test
  make_option("--trend_test", action="store_true", default=FALSE,
              help="Perform log-rank trend test for ordered quantiles [default: FALSE]"),
  
  # NEW: Pairwise comparisons
  make_option("--pairwise_comparisons", action="store_true", default=FALSE,
              help="Perform pairwise log-rank tests between all quantile groups [default: FALSE]"),
  make_option("--pairwise_method", type="character", default="holm",
              help="Multiple testing correction method: 'bonferroni', 'holm', 'BH', 'BY', 'none' [default: holm]"),
  
  # NEW: Incidence rates
  make_option("--incidence_rates", action="store_true", default=FALSE,
              help="Calculate incidence rates (events per 1000 person-years) [default: FALSE]"),
  
  # Stratification
  make_option("--stratify_by", type="character", default=NULL,
              help="Variable for stratified analyses (e.g., 'sex', 'age_strata') [default: NULL]"),
  
  # Output
  make_option("--score_role", type="character", default="all",
              help=paste("Which PRS drives the displays: 'all' (default), or a comma-separated",
                         "subset of onset,outcome,progression. The legacy behaviour was the",
                         "primary score only, which is the ONSET PRS.")),
  make_option("--min_cell_count", type="integer", default=10,
              help="Disclosure: blank risk-table counts for groups below this at-risk [default: 10]"),
  make_option("--role_quantiles", type="character", default="Q5,Q10",
              help=paste("Quantile schemes emitted per role [default: Q5,Q10]. 'all' emits every",
                         "scheme for every role. Extreme-tail schemes are always kept on the",
                         "primary score only, to bound output volume.")),
  make_option("--outdir", type="character", default="output",
              help="Output directory [default: output]"),
  make_option("--prefix", type="character", default="",
              help="Prefix for output files [default: auto]"),
  make_option("--verbose", action="store_true", default=TRUE,
              help="Print detailed output [default: TRUE]"),
  
  # KM plot
  make_option("--generate_plots", action="store_true", default=FALSE,
              help="Generate KM plots [default: FALSE]"),
  make_option("--plot_formats", type="character", default="png",
              help="Comma-separated formats to save (pdf,png) [default: png]"),
  make_option("--plot_width", type="double", default=7,
              help="Plot width in inches [default: 7]"),
  make_option("--plot_height", type="double", default=5,
              help="Plot height in inches [default: 5]"),
  make_option("--left_trunc_col", type="character", default=NULL,
              help="Optional column for delayed entry (left truncation). If provided, uses Surv(entry, time, status)"),
  make_option("--break_time_by", type="double", default=NA,
              help="Time axis major break spacing (same units as TIME) [default: auto]"),
  make_option("--min_at_risk_prop", type="double", default=0.05,
              help="Stop plotting once any stratum has < this proportion at risk [default: 0.05]"),
  make_option("--min_at_risk_n", type="integer", default=10,
              help="Stop plotting once any stratum has fewer than this many at risk [default: 10]"),
  make_option("--palette", type="character", default="Dark2",
              help="ggplot2 palette name or comma-separated hex colors [default: Dark2]"),
  make_option("--title_suffix", type="character", default="",
              help="Optional suffix appended to plot titles [default: empty]"),
  make_option("--show_conf_int", action="store_true", default=TRUE,
              help="Show 95% CI bands [default: TRUE]"),
  make_option("--show_pval", action="store_true", default=TRUE,
              help="Annotate log-rank p-value [default: TRUE]")
)

opt <- parse_args(OptionParser(
  option_list=option_list,
  description="\nPerform Kaplan-Meier survival analysis from fitted models.",
  epilogue="
Examples:
  # Basic analysis
  Rscript 02_kaplan_meier.R --models_file output/fitted_models.rds \\
    --data_file output/data_processed.rds \\
    --survival_percentiles '10,25,50,75,90' \\
    --calculate_rmst --rmst_tau '10,15' \\
    --trend_test --pairwise_comparisons \\
    --incidence_rates

  # With plots and all metrics
  Rscript 02_kaplan_meier.R --models_file output/fitted_models.rds \\
    --data_file output/data_processed.rds \\
    --generate_plots --survival_percentiles '25,50,75' \\
    --calculate_rmst --trend_test --pairwise_comparisons --incidence_rates

  # With stratification
  Rscript 02_kaplan_meier.R --models_file output/fitted_models.rds \\
    --data_file output/data_processed.rds --stratify_by sex \\
    --survival_percentiles '25,50,75' --calculate_rmst
"
))

# Validate
if (is.null(opt$models_file)) stop("--models_file is required")
if (is.null(opt$data_file)) stop("--data_file is required")
if (!file.exists(opt$models_file)) stop("Models file not found: ", opt$models_file)
if (!file.exists(opt$data_file)) stop("Data file not found: ", opt$data_file)

dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

cat("==============================================================================\n")
cat("ENHANCED KAPLAN-MEIER SURVIVAL ANALYSIS\n")
cat("==============================================================================\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("==============================================================================\n\n")

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("LOADING DATA\n")
cat("----------------------------------------------------------------------\n")

results <- readRDS(opt$models_file)
DT <- readRDS(opt$data_file)

cat("✓ Loaded fitted models\n")
cat("✓ Loaded processed data (N=", nrow(DT), ")\n\n", sep="")

metadata <- results$metadata

# ------------------------------------------------------------------------------
# Which PRS drives the displays?
# ------------------------------------------------------------------------------
# metadata$quantiles is built from the PRIMARY score, which is the ONSET PRS as
# the driver supplies it — so by default every KM curve here described the onset
# score, not the outcome or progression score the manuscript is about.
# quantiles_by_role (written by 01) carries the same schemes for all three roles;
# its keys are prefixed with the role so the existing "PRS_<key>" column lookup
# resolves to the PRS_<role>_<scheme> columns without further changes.
select_role_quantiles <- function(metadata, score_role, role_quantiles) {
  qb <- metadata$quantiles_by_role
  if (is.null(qb) || !length(qb)) {
    cat("  Note: no per-role quantiles in metadata (older 01 output); using the primary score.\n")
    return(metadata$quantiles)
  }
  roles <- if (identical(score_role, "all")) names(qb)
           else intersect(trimws(strsplit(score_role, ",")[[1]]), names(qb))
  if (!length(roles))
    stop("--score_role '", score_role, "' matched no role. Available: ",
         paste(names(qb), collapse = ", "))
  keep_q <- if (identical(role_quantiles, "all")) NULL
            else trimws(strsplit(role_quantiles, ",")[[1]])
  out <- list()
  # Extreme-tail schemes stay on the primary score only: crossing them with three
  # roles multiplies the output set several-fold for little added information.
  legacy <- metadata$quantiles
  if (!is.null(legacy)) for (k in grep("^extreme_", names(legacy), value = TRUE)) out[[k]] <- legacy[[k]]
  for (r in roles) {
    rl  <- qb[[r]]
    sel <- if (is.null(keep_q)) names(rl) else intersect(names(rl), keep_q)
    for (k in sel) out[[paste0(r, "_", k)]] <- rl[[k]]
  }
  cat("  Score roles displayed:", paste(roles, collapse = ", "),
      "| stratifications:", length(out), "\n")
  out
}
quantile_list <- select_role_quantiles(metadata, opt$score_role, opt$role_quantiles)

timeCol <- "TIME"
statusCol <- "PHENO"

# Construct prefix
if (nchar(opt$prefix) > 0) {
  prefix <- paste0(opt$prefix, "_")
} else {
  base_prefix_parts <- c()
  if (!is.null(metadata$options$cohort) && metadata$options$cohort != "") {
    base_prefix_parts <- c(base_prefix_parts, metadata$options$cohort)
  }
  if (!is.null(metadata$options$ancestry) && metadata$options$ancestry != "") {
    base_prefix_parts <- c(base_prefix_parts, metadata$options$ancestry)
  }
  if (!is.null(metadata$options$pheno_name) && metadata$options$pheno_name != "") {
    base_prefix_parts <- c(base_prefix_parts, metadata$options$pheno_name)
  }
  if (length(base_prefix_parts) > 0) {
    base_prefix <- paste(base_prefix_parts, collapse="_")
    prefix <- paste0(base_prefix, "_02_kaplan_meier_")
  } else {
    prefix <- "02_kaplan_meier_"
  }
}
cat("Output prefix:", prefix, "\n\n")

log_file <- file.path(opt$outdir, paste0(prefix, "log.txt"))
log_con <- file(log_file, open="wt")
sink(log_con, type="output", split=opt$verbose)
sink(log_con, type="message")

# Parse options
survival_percentiles <- as.numeric(trimws(strsplit(opt$survival_percentiles, ",")[[1]]))
if (!is.null(opt$rmst_tau)) {
  rmst_tau <- as.numeric(trimws(strsplit(opt$rmst_tau, ",")[[1]]))
  # --rmst_tau is the CANDIDATE set. RMST is only interpretable at horizons with
  # adequate follow-up, so intersect with the stage-01 supported horizons
  # (eval_time_points, with the usual fallbacks) before integrating. The
  # per-group support rule inside calculate_rmst() is the second gate.
  ts <- metadata$time_support
  supported <- NULL
  if (!is.null(ts)) {
    supported <- ts$eval_time_points
    if (is.null(supported) || !length(supported)) supported <- ts$requested_valid_time_points
    if (is.null(supported) || !length(supported)) supported <- ts$t_pass_auto
  }
  if (!is.null(supported) && length(supported)) {
    keep <- rmst_tau[rmst_tau <= max(supported, na.rm=TRUE)]
    dropped <- setdiff(rmst_tau, keep)
    if (length(dropped))
      cat("  RMST: tau", paste(dropped, collapse=", "),
          "exceed the stage-01 supported horizons (max", round(max(supported,na.rm=TRUE),2), "); dropped.\n")
    rmst_tau <- keep
  }
  if (!length(rmst_tau)) { cat("  RMST: no requested tau within supported horizons; disabling RMST.\n")
                           rmst_tau <- NULL; opt$calculate_rmst <- FALSE }
} else {
  rmst_tau <- NULL
}

cat("Configuration:\n")
cat("  N observations:", format(nrow(DT), big.mark=","), "\n")
cat("  N events:", format(sum(DT[[statusCol]]), big.mark=","), 
    sprintf(" (%.1f%%)\n", 100*mean(DT[[statusCol]])))
cat("  Median follow-up:", sprintf("%.2f years\n", median(DT[[timeCol]])))
cat("  Available quantiles:", paste(names(quantile_list), collapse=", "), "\n")
cat("  Survival percentiles:", paste(survival_percentiles, collapse=", "), "\n")
cat("  Calculate RMST:", opt$calculate_rmst, "\n")
if (opt$calculate_rmst) {
  if (is.null(rmst_tau)) {
    cat("  RMST tau: max follow-up\n")
  } else {
    cat("  RMST tau:", paste(rmst_tau, collapse=", "), "\n")
  }
  cat("  RMST bootstrap:", opt$rmst_bootstrap, "\n")
  if (opt$rmst_bootstrap) {
    cat("  RMST n_boot:", opt$rmst_n_boot, "\n")
  }
}
cat("  Trend test:", opt$trend_test, "\n")
cat("  Pairwise comparisons:", opt$pairwise_comparisons, "\n")
if (opt$pairwise_comparisons) {
  cat("  P-value adjustment:", opt$pairwise_method, "\n")
}
cat("  Incidence rates:", opt$incidence_rates, "\n")
cat("  Generate plots:", opt$generate_plots, "\n")
if (opt$generate_plots) {
  cat("  Plot formats:", opt$plot_formats, "\n")
}
cat("\n")

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

perform_logrank_test <- function(data, time_col, status_col, group_col) {
  formula_str <- sprintf("Surv(%s, %s) ~ %s", time_col, status_col, group_col)
  survdiff_obj <- survdiff(as.formula(formula_str), data=data)
  
  chisq <- survdiff_obj$chisq
  df <- length(survdiff_obj$n) - 1
  p_value <- pchisq(chisq, df, lower.tail=FALSE)
  
  return(list(
    chisq = chisq,
    df = df,
    p_value = p_value,
    n_groups = length(survdiff_obj$n),
    n_per_group = survdiff_obj$n,
    events_per_group = survdiff_obj$obs
  ))
}

# NEW: Calculate percentile survival times
calculate_percentile_survival <- function(data, time_col, status_col, group_col, percentiles = c(25, 50, 75)) {
  formula_str <- sprintf("Surv(%s, %s) ~ %s", time_col, status_col, group_col)
  fit <- survfit(as.formula(formula_str), data=data)
  
  # quantile.survfit(probs = p) already returns the time at which the EVENT-time
  # distribution reaches p, i.e. where survival falls to 1 - p. Converting to
  # (100 - p)/100 first therefore asked for the COMPLEMENT: the requested 25th
  # percentile was returned as the 75th and vice versa (the median, p = 0.5, was
  # unaffected, which is why this went unnoticed). Verified against an
  # exponential with a closed-form answer: for rate 0.1 the true 25th percentile
  # is 2.877, which probs = 0.25 returns and probs = 0.75 reports as 13.825.
  survival_probs <- percentiles / 100
  
  results_list <- list()
  
  # Extract strata information
  if (is.null(fit$strata)) {
    # Single group
    group_names <- "Overall"
    n_groups <- 1
  } else {
    group_names <- names(fit$strata)
    n_groups <- length(fit$strata)
  }
  
  for (i in seq_along(percentiles)) {
    perc <- percentiles[i]
    surv_prob <- survival_probs[i]
    
    # Use quantile method from survfit
    quant_times <- quantile(fit, probs = surv_prob, conf.int = TRUE)
    
    if (is.null(fit$strata)) {
      # Single group
      result <- data.table(
        group = "Overall",
        percentile = perc,
        time = as.numeric(quant_times$quantile),
        lcl = as.numeric(quant_times$lower),
        ucl = as.numeric(quant_times$upper)
      )
    } else {
      # Multiple groups
      result <- data.table(
        group = gsub("^.*=", "", group_names),
        percentile = perc,
        time = as.numeric(quant_times$quantile),
        lcl = as.numeric(quant_times$lower),
        ucl = as.numeric(quant_times$upper)
      )
    }
    
    results_list[[i]] <- result
  }
  
  return(rbindlist(results_list))
}

# NEW: Calculate RMST manually from KM curve
calculate_rmst_from_curve <- function(times, surv, tau) {
  # Area under the survival curve up to tau. A Kaplan-Meier curve is a STEP
  # function: S is constant on [t_i, t_{i+1}) at the value it takes at t_i. The
  # correct area is therefore the rectangle sum sum(diff(t) * S_at_interval_start),
  # NOT the trapezoid sum(diff(t) * (s_i + s_{i+1})/2), which cuts the corner off
  # every step and biases RMST downward. The rectangle form reproduces
  # summary(survfit, rmean = tau)$table[["rmean"]] exactly.

  # Restrict to times <= tau; prepend time 0 with survival 1
  idx <- times <= tau
  t <- c(0, times[idx])
  s <- c(1, surv[idx])

  # Extend to tau, carrying the last survival probability forward
  if (tail(t, 1) < tau) {
    t <- c(t, tau)
    s <- c(s, tail(s, 1))
  }

  # Step (rectangle) integral: value carried from the START of each interval
  rmst <- sum(diff(t) * head(s, -1))

  return(rmst)
}

# NEW: Calculate RMST with bootstrap SE
calculate_rmst <- function(data, time_col, status_col, group_col, tau = NULL, n_boot = 200,
                           reference_label = NULL, min_at_risk = 10L) {
  # If tau is NULL, use maximum follow-up time
  if (is.null(tau)) {
    tau <- max(data[[time_col]], na.rm = TRUE)
  }

  # SUPPORT-FILTER tau. RMST integrates S(t) to tau, so tau must be reached with
  # adequate data in EVERY compared group — otherwise a group with 12y of
  # follow-up would carry its last estimate flat out to 35y and the restricted
  # means would not be comparable. Drop any tau exceeding the minimum
  # group-specific maximum follow-up, or where any group has < min_at_risk still
  # at risk at tau.
  gsplit <- split(data[[time_col]], data[[group_col]])
  gsplit <- gsplit[lengths(gsplit) > 0]
  min_group_maxfu <- if (length(gsplit)) min(vapply(gsplit, max, numeric(1), na.rm = TRUE)) else Inf
  at_risk_ok <- function(tt) all(vapply(gsplit, function(x) sum(x >= tt), integer(1)) >= min_at_risk)
  tau_all <- tau
  tau <- Filter(function(tt) tt <= min_group_maxfu && at_risk_ok(tt), tau_all)
  dropped <- setdiff(tau_all, tau)
  if (length(dropped))
    cat(sprintf("  RMST: dropped tau %s (beyond min group follow-up %.2f or < %d at risk in a group).\n",
                paste(dropped, collapse=", "), min_group_maxfu, min_at_risk))
  if (!length(tau)) {
    cat("  RMST: no tau has adequate follow-up in all groups; skipping.\n")
    return(data.table())
  }

  results_list <- list()

  for (t in tau) {
    boot_rmst_keep <- NULL   # reset per tau: never reuse a previous tau's matrix
    # Fit survfit
    formula_str <- sprintf("Surv(%s, %s) ~ %s", time_col, status_col, group_col)
    fit <- survfit(as.formula(formula_str), data=data)
    # Fit survfit
    formula_str <- sprintf("Surv(%s, %s) ~ %s", time_col, status_col, group_col)
    fit <- survfit(as.formula(formula_str), data=data)
    
    # Calculate RMST for each group
    if (is.null(fit$strata)) {
      # Single group
      rmst_val <- calculate_rmst_from_curve(fit$time, fit$surv, t)
      result <- data.table(
        group = "Overall",
        tau = t,
        rmst = rmst_val,
        rmst_se = NA_real_  # SE calculation requires bootstrap
      )
    } else {
      # Multiple groups
      strata_names <- names(fit$strata)
      group_results <- list()
      
      # Calculate cumulative indices for each stratum
      cumsum_strata <- c(0, cumsum(fit$strata))
      
      for (i in seq_along(strata_names)) {
        idx_start <- cumsum_strata[i] + 1
        idx_end <- cumsum_strata[i + 1]
        
        times_i <- fit$time[idx_start:idx_end]
        surv_i <- fit$surv[idx_start:idx_end]
        
        rmst_val <- calculate_rmst_from_curve(times_i, surv_i, t)
        
        group_results[[i]] <- data.table(
          group = gsub("^.*=", "", strata_names[i]),
          tau = t,
          rmst = rmst_val,
          rmst_se = NA_real_
        )
      }
      result <- rbindlist(group_results)
    }
    
    # Bootstrap for SE (simplified - just use percentile method for CI)
    # For computational efficiency, use analytical approximation
    # SE ≈ sqrt(Var(RMST)) can be estimated from survival curve variance
    # For now, use Greenwood's formula approximation
    
    # Simple bootstrap approach for SE estimation
    if (nrow(result) > 1 && n_boot > 0) {
      tryCatch({
        # Columns are keyed by GROUP LABEL, not position. Previously the loop
        # walked seq_along(strata_names) against the replicate's own strata
        # order, so a replicate missing a group (small strata, extreme tails)
        # shifted every later column into the wrong group's SE.
        boot_rmst <- matrix(NA_real_, nrow = n_boot, ncol = nrow(result),
                            dimnames = list(NULL, as.character(result$group)))

        for (b in 1:n_boot) {
          # Bootstrap sample
          boot_idx <- sample(1:nrow(data), replace = TRUE)
          boot_data <- data[boot_idx, ]

          # Fit and calculate RMST
          boot_fit <- survfit(as.formula(formula_str), data=boot_data)

          if (!is.null(boot_fit$strata)) {
            b_names  <- gsub("^.*=", "", names(boot_fit$strata))
            b_cumsum <- c(0, cumsum(boot_fit$strata))
            for (i in seq_along(b_names)) {
              col <- match(b_names[i], colnames(boot_rmst))
              if (is.na(col)) next                      # group not in the original fit
              idx_start <- b_cumsum[i] + 1
              idx_end   <- b_cumsum[i + 1]
              if (idx_end >= idx_start) {
                times_i <- boot_fit$time[idx_start:idx_end]
                surv_i  <- boot_fit$surv[idx_start:idx_end]
                boot_rmst[b, col] <- calculate_rmst_from_curve(times_i, surv_i, t)
              }
            }
          }
        }

        # Calculate SE from bootstrap
        result[, rmst_se := apply(boot_rmst, 2, sd, na.rm = TRUE)]
        # Keep the replicate-level matrix so the DIFFERENCE can be formed within
        # each replicate below (paired), rather than combining two marginal SEs.
        boot_rmst_keep <- boot_rmst
      }, error = function(e) {
        # If bootstrap fails, leave SE as NA
        result[, rmst_se := NA_real_]
      })
    }
    
    # Add confidence intervals
    result[, rmst_lcl := rmst - 1.96 * rmst_se]
    result[, rmst_ucl := rmst + 1.96 * rmst_se]
    
    # Calculate RMST differences if reference group exists
    if (nrow(result) > 1) {
      # Resolve the reference BY LABEL. It was row 1, which silently assumes the
      # reference group sorts first; with per-role quantiles and extreme-tail
      # schemes ("Middle" as reference) that assumption does not hold.
      ref_row <- if (!is.null(reference_label) &&
                     !is.na(match(as.character(reference_label), as.character(result$group))))
                   match(as.character(reference_label), as.character(result$group)) else 1L
      ref_rmst <- result[ref_row, rmst]
      ref_se   <- result[ref_row, rmst_se]

      result[, rmst_diff := rmst - ref_rmst]

      # PAIRED bootstrap difference: form group - reference INSIDE each replicate,
      # then take the SE/CI/p from that distribution. sqrt(se^2 + ref_se^2) assumed
      # the two estimates were independent, but they come from the same resampled
      # data and are correlated.
      if (exists("boot_rmst_keep") && !is.null(boot_rmst_keep) &&
          ncol(boot_rmst_keep) == nrow(result)) {
        boot_diff <- sweep(boot_rmst_keep, 1, boot_rmst_keep[, ref_row], "-")
        d_se <- apply(boot_diff, 2, sd, na.rm = TRUE)
        d_lo <- apply(boot_diff, 2, function(z) as.numeric(stats::quantile(z, 0.025, na.rm = TRUE)))
        d_hi <- apply(boot_diff, 2, function(z) as.numeric(stats::quantile(z, 0.975, na.rm = TRUE)))
        result[, rmst_diff_se  := d_se]
        result[, rmst_diff_lcl := d_lo]
        result[, rmst_diff_ucl := d_hi]
        result[, rmst_diff_pval := 2 * pnorm(-abs(rmst_diff / rmst_diff_se))]
        result[ref_row, `:=`(rmst_diff_se = NA_real_, rmst_diff_lcl = NA_real_,
                             rmst_diff_ucl = NA_real_, rmst_diff_pval = NA_real_)]
      } else if (!is.na(ref_se)) {
        result[, rmst_diff_se := sqrt(rmst_se^2 + ref_se^2)]
        result[, rmst_diff_lcl := rmst_diff - 1.96 * rmst_diff_se]
        result[, rmst_diff_ucl := rmst_diff + 1.96 * rmst_diff_se]
        result[, rmst_diff_pval := 2 * pnorm(-abs(rmst_diff / rmst_diff_se))]
      } else {
        result[, rmst_diff_se := NA_real_]
        result[, rmst_diff_lcl := NA_real_]
        result[, rmst_diff_ucl := NA_real_]
        result[, rmst_diff_pval := NA_real_]
      }
    }
    
    results_list[[length(results_list) + 1]] <- result
  }
  
  return(rbindlist(results_list))
}

# NEW: Log-rank trend test
perform_trend_test <- function(data, time_col, status_col, group_col) {
  # Assumes group_col is ordered (e.g., Q1, Q2, Q3, Q4, Q5)
  # Assign scores to groups
  groups <- sort(unique(data[[group_col]]))
  scores <- seq_along(groups)
  
  # Create score variable
  data_temp <- copy(data)
  data_temp[, score := match(get(group_col), groups)]
  
  # Perform trend test using survdiff with rho=1 (or rho=0 for log-rank)
  # Actually, for trend test, we use the score as a continuous variable
  formula_str <- sprintf("Surv(%s, %s) ~ score", time_col, status_col)
  
  # Use Cox model to get trend test
  cox_fit <- coxph(as.formula(formula_str), data=data_temp)
  
  # Extract results
  coef <- summary(cox_fit)$coefficients["score", "coef"]
  se <- summary(cox_fit)$coefficients["score", "se(coef)"]
  z <- summary(cox_fit)$coefficients["score", "z"]
  p_value <- summary(cox_fit)$coefficients["score", "Pr(>|z|)"]
  
  # Also perform survdiff for comparison
  formula_str2 <- sprintf("Surv(%s, %s) ~ %s", time_col, status_col, group_col)
  survdiff_obj <- survdiff(as.formula(formula_str2), data=data, rho=0)
  
  # Score test for trend
  obs <- survdiff_obj$obs
  exp <- survdiff_obj$exp
  var_matrix <- survdiff_obj$var
  
  # Linear trend statistic
  trend_stat <- sum(scores * (obs - exp))^2 / sum(scores^2 * diag(var_matrix))
  trend_p <- pchisq(trend_stat, df=1, lower.tail=FALSE)
  
  return(list(
    n_groups = length(groups),
    groups = paste(groups, collapse=", "),
    trend_chisq = trend_stat,
    trend_p = trend_p,
    cox_coef = coef,
    cox_se = se,
    cox_z = z,
    cox_p = p_value
  ))
}

# NEW: Pairwise log-rank comparisons
perform_pairwise_comparisons <- function(data, time_col, status_col, group_col, method = "holm") {
  groups <- sort(unique(data[[group_col]]))
  n_groups <- length(groups)
  
  if (n_groups < 2) {
    return(NULL)
  }
  
  # Generate all pairwise combinations
  pairs <- combn(groups, 2, simplify = FALSE)
  
  results_list <- list()
  
  for (pair in pairs) {
    group1 <- pair[1]
    group2 <- pair[2]
    
    # Subset data
    subset_data <- data[get(group_col) %in% c(group1, group2)]
    
    # Perform log-rank test
    formula_str <- sprintf("Surv(%s, %s) ~ %s", time_col, status_col, group_col)
    survdiff_obj <- survdiff(as.formula(formula_str), data=subset_data)
    
    chisq <- survdiff_obj$chisq
    p_value <- pchisq(chisq, df=1, lower.tail=FALSE)
    
    result <- data.table(
      group1 = group1,
      group2 = group2,
      n1 = survdiff_obj$n[1],
      n2 = survdiff_obj$n[2],
      events1 = survdiff_obj$obs[1],
      events2 = survdiff_obj$obs[2],
      chisq = chisq,
      p_value = p_value
    )
    
    results_list[[length(results_list) + 1]] <- result
  }
  
  results <- rbindlist(results_list)
  
  # Apply multiple testing correction
  if (method != "none") {
    results[, p_adjusted := p.adjust(p_value, method = method)]
  } else {
    results[, p_adjusted := p_value]
  }
  
  return(results)
}

# NEW: Calculate incidence rates
calculate_incidence_rates <- function(data, time_col, status_col, group_col) {
  # Calculate person-years and events by group
  data_temp <- copy(data)
  
  results <- data_temp[, .(
    n = .N,
    events = sum(get(status_col)),
    person_years = sum(get(time_col)),
    mean_follow_up = mean(get(time_col))
  ), by = group_col]
  
  setnames(results, group_col, "group")
  
  # Calculate incidence rate per 1000 person-years
  results[, incidence_rate := (events / person_years) * 1000]
  
  # Calculate 95% CI using Poisson approximation
  results[, ir_lcl := (qchisq(0.025, 2*events) / (2*person_years)) * 1000]
  results[, ir_ucl := (qchisq(0.975, 2*(events+1)) / (2*person_years)) * 1000]
  
  # Calculate rate ratios (using first group as reference)
  if (nrow(results) > 1) {
    ref_rate <- results[1, incidence_rate]
    results[, rate_ratio := incidence_rate / ref_rate]
    
    # Calculate RR confidence intervals using log transformation
    results[, rr_lcl := exp(log(rate_ratio) - 1.96 * sqrt(1/events + 1/results[1, events]))]
    results[, rr_ucl := exp(log(rate_ratio) + 1.96 * sqrt(1/events + 1/results[1, events]))]
  }
  
  return(results)
}

create_risk_table <- function(data, time_col, status_col, group_col, times, min_cell_count = 10L) {
  formula_str <- sprintf("Surv(%s, %s) ~ %s", time_col, status_col, group_col)
  fit <- survfit(as.formula(formula_str), data=data)

  risk_summary <- summary(fit, times=times)

  result <- data.table(
    group = rep(gsub("^.*=", "", names(fit$strata)), times=length(times)),
    time = rep(times, each=length(fit$strata)),
    n_risk = risk_summary$n.risk,
    n_event = risk_summary$n.event,
    n_censor = risk_summary$n.censor,
    surv = risk_summary$surv,
    surv_lower = risk_summary$lower,
    surv_upper = risk_summary$upper
  )
  # Disclosure: suppress when ANY reported count (at-risk, events, censored) is a
  # small non-zero integer, not only n_risk — n_risk=100 with n_event=2 still
  # discloses the 2. Blanks the counts and the survival estimates.
  suppress_small_cells(result,
                       count_cols = c("n_risk","n_event","n_censor"),
                       threshold = min_cell_count,
                       estimate_cols = c("surv","surv_lower","surv_upper"))
  return(result)
}

# Build a Surv object, supporting optional left truncation
make_surv <- function(data, time_col, status_col, entry_col=NULL) {
  if (!is.null(entry_col) && entry_col %in% names(data)) {
    Surv(time=data[[entry_col]], time2=data[[time_col]], event=data[[status_col]])
  } else {
    Surv(time=data[[time_col]], event=data[[status_col]])
  }
}

# Determine a sensible truncation time where groups are still well represented
compute_trunc_time <- function(fit, min_prop=0.05, min_n=10) {
  # fit is a survfit object (by groups)
  ss <- summary(fit)
  # n.risk comes concatenated across strata; we compute per-time min proportion & count
  # Retrieve group sizes
  strata_sizes <- as.numeric(sapply(strsplit(names(fit$strata), "="), function(x) tail(x, 1)))
  # Above fails; safer: get sizes directly from summary table if available
  if (!is.null(fit$n)) strata_sizes <- fit$n
  # Build a data.frame with time, strata index, n.risk
  # 'ss$strata' repeats names; map them to an index
  strata_names <- names(fit$strata)
  if (is.null(strata_names)) {
    # Single stratum
    prop_ok <- (ss$n.risk / max(ss$n.risk)) >= min_prop & ss$n.risk >= min_n
    return(max(ss$time[prop_ok], na.rm=TRUE))
  }
  strata_index <- match(ss$strata, strata_names)
  max_n_per_stratum <- tapply(ss$n.risk, strata_index, max, na.rm=TRUE)
  # For each row, compute proportion within its stratum
  prop <- mapply(function(nr, si) nr / max_n_per_stratum[si], ss$n.risk, strata_index)
  # Find the last time where **all** strata meet thresholds
  ok <- prop >= min_prop & ss$n.risk >= min_n
  # times are per row; need the largest time such that all strata have ok at or before time
  # Conservative: take min across strata of their last ok time, then overall max of those mins
  last_ok_by_stratum <- tapply(ss$time[ok], strata_index[ok], function(v) if (length(v)) max(v) else 0)
  trunc_time <- min(unlist(last_ok_by_stratum), na.rm=TRUE)
  if (!is.finite(trunc_time)) trunc_time <- max(ss$time, na.rm=TRUE)
  trunc_time
}

# Parse palette option
get_palette <- function(pal_str, k) {
  if (grepl("#", pal_str)) {
    cols <- trimws(strsplit(pal_str, ",")[[1]])
    if (length(cols) < k) cols <- rep(cols, length.out=k)
    return(cols[1:k])
  } else {
    return(scales::hue_pal()(k)) # fallback
  }
}

plot_km_with_risktable <- function(data, time_col, status_col, group_col,
                                   entry_col = NULL, title = "", breaks_by = NULL,
                                   conf_int = TRUE, pval_text = NULL, palette = "Dark2",
                                   trunc_prop = 0.05, trunc_n = 10) {
  
  stopifnot(time_col %in% names(data), status_col %in% names(data), group_col %in% names(data))
  if (!is.null(entry_col) && !entry_col %in% names(data)) entry_col <- NULL
  
  # Build a minimal data.frame with fixed names used in the literal formula
  df <- data.frame(
    TIME   = data[[time_col]],
    STATUS = as.integer(data[[status_col]]),   # force 0/1
    GROUP  = factor(data[[group_col]])
  )
  if (!is.null(entry_col)) df$ENTRY <- data[[entry_col]]
  
  # ---- FIT KM with a *literal* formula (no variables, no as.formula) ----
  if (!is.null(entry_col)) {
    fit <- survival::survfit(Surv(ENTRY, TIME, STATUS) ~ GROUP, data = df)
  } else {
    fit <- survival::survfit(Surv(TIME, STATUS) ~ GROUP, data = df)
  }
  
  # Palette
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) stop("Please install RColorBrewer")
  K <- if (is.null(fit$strata)) 1L else length(fit$strata)
  pal_cols <- if (grepl("#", palette)) {
    cols <- trimws(strsplit(palette, ",")[[1]]); rep(cols, length.out = K)
  } else if (palette %in% rownames(RColorBrewer::brewer.pal.info)) {
    # Get max colors available for this palette
    max_colors <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
    
    if (K <= max_colors) {
      # If we need fewer colors than available, request them directly
      RColorBrewer::brewer.pal(max(3, K), palette)
    } else {
      # If we need more colors than available, use colorRampPalette to interpolate
      base_colors <- RColorBrewer::brewer.pal(max_colors, palette)
      colorRampPalette(base_colors)(K)
    }
  } else {
    scales::hue_pal()(K)
  }
  
  # Robust truncation time
  compute_trunc_time <- function(fit, min_prop = 0.05, min_n = 10) {
    ss <- summary(fit)
    if (is.null(fit$strata)) {
      n0 <- max(ss$n.risk, na.rm = TRUE)
      ok <- (ss$n.risk >= min_n) & (ss$n.risk / n0 >= min_prop)
      return(ifelse(any(ok), max(ss$time[ok], na.rm = TRUE), max(ss$time, na.rm = TRUE)))
    }
    strata_names <- names(fit$strata)
    idx <- as.integer(factor(ss$strata, levels = strata_names))
    n0_per <- as.numeric(fit$n)
    prop <- ss$n.risk / n0_per[idx]
    ok <- (ss$n.risk >= min_n) & (prop >= min_prop)
    if (!any(ok)) return(max(ss$time, na.rm = TRUE))
    last_ok_by <- tapply(ss$time[ok], idx[ok], function(x) if (length(x)) max(x) else 0)
    if (any(unlist(last_ok_by) > 0)) min(unlist(last_ok_by)) else max(ss$time, na.rm = TRUE)
  }
  t_cut <- compute_trunc_time(fit, min_prop = trunc_prop, min_n = trunc_n)
  
  # ---- PLOT (pass the same df so the risk table code has data) ----
  gp <- survminer::ggsurvplot(
    fit,
    data = df,
    conf.int = conf_int,
    risk.table = TRUE,
    risk.table.height = 0.22,
    break.time.by = breaks_by,
    xlim = c(0, t_cut),
    palette = pal_cols,
    censor = TRUE,
    legend.title = group_col,          # show original name
    legend.labs = levels(df$GROUP),
    ggtheme = ggplot2::theme_minimal(base_size = 12)
  )
  
  gp$plot <- gp$plot +
    ggplot2::labs(title = title, x = time_col, y = "Survival probability") +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   legend.position = "right")
  
  if (!is.null(pval_text)) {
    gp$plot <- gp$plot +
      ggplot2::annotate("text", x = t_cut * 0.6, y = 0.05, label = pval_text, hjust = 0)
  }
  
  gp
}

save_km <- function(gp, outbase, width=7, height=8.5, formats=c("png")) {
  # One figure (KM curves + risk table). Cairo preferred, grDevices fallback so
  # headless machines still produce figures; never aborts.
  have_cairo <- requireNamespace("Cairo", quietly = TRUE)
  for (fmt in formats) {
    if (!fmt %in% c("png","pdf")) next
    outfile <- paste0(outbase, ".", fmt)
    ok <- FALSE
    if (have_cairo) ok <- tryCatch({
      if (fmt == "pdf") Cairo::CairoPDF(outfile, width=width, height=height)
      else              Cairo::CairoPNG(outfile, width=width, height=height, units="in", dpi=300)
      print(gp); grDevices::dev.off(); TRUE
    }, error = function(e) { try(grDevices::dev.off(), silent=TRUE); FALSE })
    if (!ok) tryCatch({
      if (fmt == "pdf") grDevices::pdf(outfile, width=width, height=height)
      else              grDevices::png(outfile, width=width, height=height, units="in", res=300)
      print(gp); grDevices::dev.off()
      if (!have_cairo) message("Cairo unavailable; wrote ", basename(outfile), " via grDevices fallback")
    }, error = function(e) warning("Could not save KM plot ", basename(outfile), ": ", conditionMessage(e)))
  }
}


# ==============================================================================
# QUANTILE SELECTION
# ==============================================================================

cat("QUANTILE SELECTION\n")
cat("----------------------------------------------------------------------\n")

if (opt$quantiles == "all") {
  quantiles_to_analyze <- names(quantile_list)
} else {
  quantiles_to_analyze <- trimws(strsplit(opt$quantiles, ",")[[1]])
  
  invalid <- setdiff(quantiles_to_analyze, names(quantile_list))
  if (length(invalid) > 0) {
    stop("Invalid quantiles: ", paste(invalid, collapse=", "), "\n",
         "Available: ", paste(names(quantile_list), collapse=", "))
  }
}

cat("Analyzing:\n")
for (q in quantiles_to_analyze) {
  q_obj <- quantile_list[[q]]
  if (grepl("^extreme_", q)) {
    cat("  -", q, "(", q_obj$n_groups, "groups)\n")
  } else {
    cat("  -", q, "(", q_obj$n_quantiles, "quantiles)\n")
  }
}
cat("\n")

# ==============================================================================
# TIME PARAMETERS
# ==============================================================================

max_time <- max(DT[[timeCol]])
if (is.na(opt$break_time_by)) {
  # Choose a sensible default based on the data
  time_breaks <- pretty(c(0, max_time), n = 5)
  opt$break_time_by <- time_breaks[2]
} else {
  time_breaks <- pretty(c(0, max_time), n = 5)
}

if (!is.null(opt$risk_table_times)) {
  risk_times <- as.numeric(trimws(strsplit(opt$risk_table_times, ",")[[1]]))
} else {
  # Robust horizon selection. Raw max(time) is fragile: a few outlier long censoring
  # times push the horizons past the informative follow-up, so the KM risk table can
  # come back mostly/entirely NA (e.g. a short prospective run capped at 20/40/60y by
  # pretty()). Cap at the 90th-percentile follow-up so outliers cannot extend the grid,
  # and PREFER the requested evaluation grid (the driver passes --time-points via
  # --rmst_tau) so the risk table uses the SAME horizons as discrimination/calibration.
  robust_max <- suppressWarnings(as.numeric(stats::quantile(DT[[timeCol]], 0.90, na.rm = TRUE)))
  if (!is.finite(robust_max) || robust_max <= 0) robust_max <- max_time
  grid <- NULL
  if (!is.null(opt$rmst_tau) && nzchar(opt$rmst_tau))
    grid <- suppressWarnings(as.numeric(trimws(strsplit(opt$rmst_tau, ",")[[1]])))
  if (is.null(grid) || !any(is.finite(grid))) grid <- time_breaks   # fall back to pretty() breaks
  risk_times <- sort(unique(grid[is.finite(grid) & grid > 0 & grid <= robust_max]))
  if (!length(risk_times)) {                       # data shorter than the whole grid
    pb <- pretty(c(0, robust_max), n = 5); risk_times <- pb[pb > 0 & pb <= robust_max]
  }
  if (!length(risk_times)) risk_times <- robust_max
  cat(sprintf("  (risk-table horizons: robust 90th-pct follow-up cap = %.1f; using %s)\n",
              robust_max, paste(risk_times, collapse=", ")))
}

cat("TIME PARAMETERS\n")
cat("----------------------------------------------------------------------\n")
cat("  Max time:", max_time, "\n")
cat("  Risk table times:", paste(risk_times, collapse=", "), "\n\n")

# ==============================================================================
# MAIN ANALYSIS
# ==============================================================================

cat("==============================================================================\n")
cat("KAPLAN-MEIER ANALYSIS\n")
cat("==============================================================================\n\n")

all_logrank <- list()
all_percentile_surv <- list()
all_rmst <- list()
all_trend_tests <- list()
all_pairwise <- list()
all_incidence_rates <- list()
all_risk_tables <- list()

for (q_name in quantiles_to_analyze) {
  cat("Analyzing:", q_name, "\n")
  cat("----------------------------------------------------------------------\n")
  
  if (grepl("^extreme_", q_name)) {
    quant_col <- paste0("PRS_", q_name)
  } else {
    quant_col <- paste0("PRS_", q_name)
  }
  
  if (!quant_col %in% names(DT)) {
    cat("  ⚠ Column not found. Skipping.\n\n")
    next
  }
  
  # Log-rank test
  cat("  Performing log-rank test...\n")
  logrank <- perform_logrank_test(DT, timeCol, statusCol, quant_col)
  
  logrank_dt <- data.table(
    quantile = q_name,
    chisq = logrank$chisq,
    df = logrank$df,
    p_value = logrank$p_value,
    n_groups = logrank$n_groups
  )
  
  all_logrank[[q_name]] <- logrank_dt
  
  cat("    χ²=", sprintf("%.3f", logrank$chisq), 
      ", df=", logrank$df,
      ", p=", format.pval(logrank$p_value, digits=3), "\n", sep="")
  
  # Percentile survival
  cat("  Calculating percentile survival times...\n")
  percentile_surv <- calculate_percentile_survival(DT, timeCol, statusCol, quant_col, 
                                                   percentiles = survival_percentiles)
  percentile_surv[, quantile := q_name]
  all_percentile_surv[[q_name]] <- percentile_surv
  
  print(percentile_surv[, .(group, percentile, time, lcl, ucl)])
  
  # RMST
  if (opt$calculate_rmst) {
    cat("  Calculating RMST...\n")
    n_boot <- if (opt$rmst_bootstrap) opt$rmst_n_boot else 0
    # Reference resolved by LABEL: quantile objects carry $reference ("Q1" for
    # standard schemes, "Middle" for extreme-tail schemes), which is not always
    # the first row of the RMST table.
    qo <- quantile_list[[q_name]]
    ref_lab <- if (!is.null(qo) && !is.null(qo$reference)) as.character(qo$reference) else NULL
    rmst_result <- calculate_rmst(DT, timeCol, statusCol, quant_col, tau = rmst_tau, n_boot = n_boot,
                                  reference_label = ref_lab)
    rmst_result[, quantile := q_name]
    all_rmst[[q_name]] <- rmst_result
    
    print(rmst_result[, .(group, tau, rmst, rmst_lcl, rmst_ucl)])
    if ("rmst_diff" %in% names(rmst_result)) {
      cat("  RMST differences vs reference:\n")
      print(rmst_result[, .(group, tau, rmst_diff, rmst_diff_lcl, rmst_diff_ucl, rmst_diff_pval)])
    }
  }
  
  # Trend test
  if (opt$trend_test) {
    cat("  Performing trend test...\n")
    trend_result <- perform_trend_test(DT, timeCol, statusCol, quant_col)
    trend_dt <- data.table(
      quantile = q_name,
      n_groups = trend_result$n_groups,
      groups = trend_result$groups,
      trend_chisq = trend_result$trend_chisq,
      trend_p = trend_result$trend_p,
      cox_coef = trend_result$cox_coef,
      cox_se = trend_result$cox_se,
      cox_z = trend_result$cox_z,
      cox_p = trend_result$cox_p
    )
    all_trend_tests[[q_name]] <- trend_dt
    
    cat("    Trend χ²=", sprintf("%.3f", trend_result$trend_chisq),
        ", p=", format.pval(trend_result$trend_p, digits=3), "\n", sep="")
  }
  
  # Pairwise comparisons
  if (opt$pairwise_comparisons) {
    cat("  Performing pairwise comparisons...\n")
    pairwise_result <- perform_pairwise_comparisons(DT, timeCol, statusCol, quant_col, 
                                                    method = opt$pairwise_method)
    if (!is.null(pairwise_result)) {
      pairwise_result[, quantile := q_name]
      all_pairwise[[q_name]] <- pairwise_result
      
      # Print significant comparisons
      sig_pairs <- pairwise_result[p_adjusted < 0.05]
      if (nrow(sig_pairs) > 0) {
        cat("    Significant pairs (adjusted p < 0.05):\n")
        for (i in 1:nrow(sig_pairs)) {
          cat("      ", sig_pairs[i, group1], " vs ", sig_pairs[i, group2], 
              ": p=", format.pval(sig_pairs[i, p_adjusted], digits=3), "\n", sep="")
        }
      } else {
        cat("    No significant pairwise differences after correction\n")
      }
    }
  }
  
  # Incidence rates
  if (opt$incidence_rates) {
    cat("  Calculating incidence rates...\n")
    incidence_result <- calculate_incidence_rates(DT, timeCol, statusCol, quant_col)
    incidence_result[, quantile := q_name]
    all_incidence_rates[[q_name]] <- incidence_result
    
    print(incidence_result[, .(group, n, events, person_years, incidence_rate, ir_lcl, ir_ucl)])
    if ("rate_ratio" %in% names(incidence_result)) {
      cat("  Rate ratios vs reference:\n")
      print(incidence_result[, .(group, rate_ratio, rr_lcl, rr_ucl)])
    }
  }
  
  # Risk table
  cat("  Creating risk table...\n")
  risk_table <- create_risk_table(DT, timeCol, statusCol, quant_col, risk_times, min_cell_count = opt$min_cell_count)
  risk_table[, quantile := q_name]
  all_risk_tables[[q_name]] <- risk_table
  
  cat("\n")
}

# Combine and save
# fill=TRUE throughout: with per-role stratifications a sparse role x quantile
# combination can yield a narrower table, and an un-filled rbindlist aborts the
# entire stage rather than reporting that one stratum.
logrank_combined <- rbindlist(all_logrank, fill=TRUE)
percentile_surv_combined <- rbindlist(all_percentile_surv, fill=TRUE)
risk_table_combined <- rbindlist(all_risk_tables, fill=TRUE)

fwrite(logrank_combined, file.path(opt$outdir, paste0(prefix, "logrank_tests.csv")))
fwrite(percentile_surv_combined, file.path(opt$outdir, paste0(prefix, "percentile_survival.csv")))
fwrite(risk_table_combined, file.path(opt$outdir, paste0(prefix, "risk_table.csv")))

cat("✓ Saved main results\n")

# Save optional results
if (length(all_rmst) > 0) {
  rmst_combined <- rbindlist(all_rmst, fill=TRUE)
  fwrite(rmst_combined, file.path(opt$outdir, paste0(prefix, "rmst.csv")))
  cat("✓ Saved RMST results\n")
}

if (length(all_trend_tests) > 0) {
  trend_combined <- rbindlist(all_trend_tests, fill=TRUE)
  fwrite(trend_combined, file.path(opt$outdir, paste0(prefix, "logrank_trend.csv")))
  cat("✓ Saved trend test results\n")
}

if (length(all_pairwise) > 0) {
  pairwise_combined <- rbindlist(all_pairwise, fill=TRUE)
  fwrite(pairwise_combined, file.path(opt$outdir, paste0(prefix, "pairwise_comparisons.csv")))
  cat("✓ Saved pairwise comparison results\n")
}

if (length(all_incidence_rates) > 0) {
  incidence_combined <- rbindlist(all_incidence_rates, fill=TRUE)
  # Disclosure: suppress a group with a small event count (and its derived rates,
  # which would otherwise back-reveal the count from person-years).
  suppress_small_cells(incidence_combined, count_cols = "events",
                       threshold = opt$min_cell_count,
                       estimate_cols = intersect(c("person_years","mean_follow_up","incidence_rate",
                                                   "ir_lcl","ir_ucl","rate_ratio","rr_lcl","rr_ucl"),
                                                 names(incidence_combined)))
  fwrite(incidence_combined, file.path(opt$outdir, paste0(prefix, "incidence_rates.csv")))
  cat("✓ Saved incidence rate results\n")
}

cat("\n")

# ==============================================================================
# STRATIFIED ANALYSES
# ==============================================================================

if (!is.null(opt$stratify_by)) {
  cat("==============================================================================\n")
  cat("STRATIFIED ANALYSES\n")
  cat("==============================================================================\n\n")
  
  strata_vars <- trimws(strsplit(opt$stratify_by, ",")[[1]])
  strata_vars <- intersect(strata_vars, names(DT))
  
  if (length(strata_vars) == 0) {
    cat("No valid stratification variables.\n\n")
  } else {
    cat("Stratifying by:", paste(strata_vars, collapse=", "), "\n\n")
    
    all_stratified <- list()
    
    for (strata_var in strata_vars) {
      cat("Variable:", strata_var, "\n")
      
      strata_levels <- unique(DT[[strata_var]])
      strata_levels <- strata_levels[!is.na(strata_levels)]
      
      if (length(strata_levels) < 2 || length(strata_levels) > 5) {
        cat("  ⚠ Skipping (", length(strata_levels), " levels)\n\n", sep="")
        next
      }
      
      for (q_name in quantiles_to_analyze) {
        quant_col <- ifelse(grepl("^extreme_", q_name), 
                            paste0("PRS_", q_name), 
                            paste0("PRS_", q_name))
        
        if (!quant_col %in% names(DT)) next
        
        for (level in strata_levels) {
          subset_data <- DT[get(strata_var) == level]
          
          if (nrow(subset_data) < 50 || sum(subset_data[[statusCol]]) < 10) {
            cat("  ", level, " ", q_name, ": insufficient data\n", sep="")
            next
          }
          
          logrank_strata <- perform_logrank_test(subset_data, timeCol, statusCol, quant_col)
          percentile_strata <- calculate_percentile_survival(subset_data, timeCol, statusCol, quant_col,
                                                             percentiles = survival_percentiles)
          
          cat("  ", level, " ", q_name, ": χ²=", sprintf("%.2f", logrank_strata$chisq),
              ", p=", format.pval(logrank_strata$p_value, digits=3), "\n", sep="")
          
          # Save stratified results
          strat_result <- data.table(
            strata_var = strata_var,
            strata_level = as.character(level),
            quantile = q_name,
            group = percentile_strata$group,
            percentile = percentile_strata$percentile,
            survival_time = percentile_strata$time,
            survival_lcl = percentile_strata$lcl,
            survival_ucl = percentile_strata$ucl,
            logrank_chisq = logrank_strata$chisq,
            logrank_df = logrank_strata$df,
            logrank_p = logrank_strata$p_value
          )
          
          all_stratified[[paste(strata_var, level, q_name, sep="_")]] <- strat_result
        }
      }
      cat("\n")
    }
    
    if (length(all_stratified) > 0) {
      stratified_combined <- rbindlist(all_stratified, fill=TRUE)
      # Disclosure: suppress small per-stratum count cells.
      suppress_small_cells(stratified_combined,
                           count_cols = intersect(c("n_risk","n_event","n_censor","N","Events","events"),
                                                  names(stratified_combined)),
                           threshold = opt$min_cell_count)
      fwrite(stratified_combined, file.path(opt$outdir, paste0(prefix, "stratified_results.csv")))
      cat("✓ Saved stratified results\n\n")
    }
  }
}

cat("==============================================================================\n")
cat("EXTRACTING KM CURVE DATA\n")
cat("==============================================================================\n\n")

# Extract and save KM curve data for all quantiles
all_km_data <- list()

for (q_name in quantiles_to_analyze) {
  quant_col <- paste0("PRS_", q_name)
  if (!quant_col %in% names(DT)) next
  
  cat("Extracting KM data for", q_name, "...\n")
  
  # Fit KM model
  formula_str <- sprintf("Surv(%s, %s) ~ %s", timeCol, statusCol, quant_col)
  km_fit <- survfit(as.formula(formula_str), data=DT)
  
  # Extract survival curve data
  km_summary <- summary(km_fit)
  
  # Create data table with all curve information
  n_strata <- length(km_fit$strata)
  
  if (n_strata > 0) {
    # Multiple strata
    strata_names <- names(km_fit$strata)
    strata_sizes <- km_fit$strata
    
    km_data <- data.table(
      quantile = q_name,
      stratum = rep(strata_names, strata_sizes),
      time = km_summary$time,
      n_risk = km_summary$n.risk,
      n_event = km_summary$n.event,
      n_censor = km_summary$n.censor,
      surv = km_summary$surv,
      surv_lower = km_summary$lower,
      surv_upper = km_summary$upper,
      cumulative_events = cumsum(km_summary$n.event)
    )
  } else {
    # Single stratum
    km_data <- data.table(
      quantile = q_name,
      stratum = "Overall",
      time = km_summary$time,
      n_risk = km_summary$n.risk,
      n_event = km_summary$n.event,
      n_censor = km_summary$n.censor,
      surv = km_summary$surv,
      surv_lower = km_summary$lower,
      surv_upper = km_summary$upper,
      cumulative_events = cumsum(km_summary$n.event)
    )
  }
  
  all_km_data[[q_name]] <- km_data
  
  # Save individual quantile KM data
  fwrite(km_data, file.path(opt$outdir, paste0(prefix, "km_curve_data_", q_name, ".csv")))
  cat("  ✓ Saved KM curve data for", q_name, "\n")
}

# Combine all KM data
if (length(all_km_data) > 0) {
  km_combined <- rbindlist(all_km_data, fill=TRUE)
  fwrite(km_combined, file.path(opt$outdir, paste0(prefix, "km_curve_data_all.csv")))
  cat("✓ Saved combined KM curve data\n\n")
}

cat("==============================================================================\n")
cat("PLOTTING KM CURVES\n")
cat("==============================================================================\n\n")

if (!opt$generate_plots) {
  cat("Plotting disabled (use --generate_plots to enable)\n")
  cat("KM curve data saved to CSV files\n\n")
} else {
  cat("Generating KM plots...\n\n")
  
  plot_formats <- trimws(strsplit(opt$plot_formats, ",")[[1]])
  if (length(plot_formats) == 0) plot_formats <- c("png")
  
  for (q_name in quantiles_to_analyze) {
    quant_col <- paste0("PRS_", q_name)
    if (!quant_col %in% names(DT)) next
    
    cat("Plotting", q_name, "...\n")
    
    # Title and (optional) p-value annotation
    this_p <- tryCatch({
      logrank_combined[quantile == q_name]$p_value[1]
    }, error=function(e) NA_real_)
    pval_text <- if (opt$show_pval && is.finite(this_p)) {
      paste0("Log-rank p = ", format.pval(this_p, digits=3))
    } else NULL
    
    main_title <- paste0("KM by ", quant_col,
                         ifelse(nchar(opt$title_suffix)>0, paste0(" — ", opt$title_suffix), ""))
    
    gp <- plot_km_with_risktable(
      data = DT,
      time_col = timeCol,
      status_col = statusCol,
      group_col = quant_col,
      entry_col = opt$left_trunc_col,
      title = main_title,
      breaks_by = opt$break_time_by,
      conf_int = opt$show_conf_int,
      pval_text = pval_text,
      palette = opt$palette,
      trunc_prop = opt$min_at_risk_prop,
      trunc_n = opt$min_at_risk_n
    )
    
    outbase <- file.path(opt$outdir, paste0(prefix, "km_", q_name))
    save_km(gp, outbase, width=opt$plot_width, height=opt$plot_height, formats=plot_formats)
    
    cat("✓ Saved KM plots for", q_name, "→", outbase, ".{", paste(plot_formats, collapse=","), "}\n")
  }
  
  cat("\n")
  
  # Stratified plots
  if (!is.null(opt$stratify_by) && exists("strata_vars") && length(strata_vars) > 0) {
    for (strata_var in strata_vars) {
      for (q_name in quantiles_to_analyze) {
        quant_col <- paste0("PRS_", q_name)
        if (!quant_col %in% names(DT)) next
        for (level in unique(DT[[strata_var]])) {
          subset_data <- DT[get(strata_var) == level]
          if (nrow(subset_data) < 50 || sum(subset_data[[statusCol]]) < 5) next
          
          this_p <- tryCatch({
            logrank_combined[quantile==q_name]$p_value[1]
          }, error=function(e) NA_real_)
          pval_text <- if (opt$show_pval && is.finite(this_p)) {
            paste0("Log-rank p = ", format.pval(this_p, digits=3))
          } else NULL
          
          title2 <- paste0("KM by ", quant_col, " — ", strata_var, "=", level,
                           ifelse(nchar(opt$title_suffix)>0, paste0(" — ", opt$title_suffix), ""))
          
          gp <- plot_km_with_risktable(
            data = subset_data,
            time_col = timeCol,
            status_col = statusCol,
            group_col = quant_col,
            entry_col = opt$left_trunc_col,
            title = title2,
            breaks_by = opt$break_time_by,
            conf_int = opt$show_conf_int,
            pval_text = pval_text,
            palette = opt$palette,
            trunc_prop = opt$min_at_risk_prop,
            trunc_n = opt$min_at_risk_n
          )
          
          outbase <- file.path(opt$outdir, paste0(prefix, "km_", q_name, "_", strata_var, "_", level))
          save_km(gp, outbase, width=opt$plot_width, height=opt$plot_height, formats=plot_formats)
          cat("  ✓ Saved stratified KM:", strata_var, "=", level, "—", q_name, "\n")
        }
      }
    }
  }
}


# ==============================================================================
# SUMMARY
# ==============================================================================

cat("==============================================================================\n")
cat("SUMMARY\n")
cat("==============================================================================\n\n")

cat("✓ Analysis completed!\n\n")

cat("Results for", length(quantiles_to_analyze), "quantile(s):\n")
for (q in quantiles_to_analyze) {
  logrank_res <- logrank_combined[quantile == q]
  cat("  -", q, ": p=", format.pval(logrank_res$p_value, digits=3), "\n", sep="")
}

cat("\nOutput files:\n")
cat("  - Log-rank tests\n")
cat("  - Percentile survival times (", paste(survival_percentiles, collapse=", "), "th percentiles)\n")
if (opt$calculate_rmst) {
  cat("  - RMST (Restricted Mean Survival Time)\n")
}
if (opt$trend_test) {
  cat("  - Log-rank trend test\n")
}
if (opt$pairwise_comparisons) {
  cat("  - Pairwise comparisons (", opt$pairwise_method, " adjustment)\n", sep="")
}
if (opt$incidence_rates) {
  cat("  - Incidence rates (per 1000 person-years)\n")
}
cat("  - Risk tables\n")
cat("  - KM curve data (CSV)\n")
if (opt$generate_plots) {
  cat("  - KM plots (", paste(plot_formats, collapse=", "), ")\n", sep="")
} else {
  cat("  - KM plots: SKIPPED (use --generate_plots to enable)\n")
}
if (!is.null(opt$stratify_by)) cat("  - Stratified results\n")

cat("\n")
cat("==============================================================================\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("==============================================================================\n")

sink(type="message")
sink(type="output")
close(log_con)

cat("\nLog:", log_file, "\n")