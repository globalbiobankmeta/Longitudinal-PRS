#!/usr/bin/env Rscript
# ==============================================================================
# 04_reclassification.R - Reclassification and Improvement Metrics
# ==============================================================================

# ------------------------------- PACKAGES -------------------------------------
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
required_packages <- c("optparse", "data.table", "survival", "ggplot2", "gridExtra", "ggpubr")
invisible(lapply(required_packages, load_package))  # Cairo intentionally omitted (optional; see prs_risk_utils.R)

# Shared risk-prediction helper, located next to this script. Prefer PRS_SCRIPT_DIR
# (exported by the driver) so we never parse --file, whose spaces some macOS/OneDrive
# setups encode as ~+~ (which breaks normalizePath()). Fall back to --file (with ~+~
# decoded), then to getwd().
resolve_script_dir <- function() {
  env_dir <- Sys.getenv("PRS_SCRIPT_DIR", unset = "")
  if (nzchar(env_dir)) return(normalizePath(env_dir, winslash = "/", mustWork = TRUE))
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (!length(file_arg)) return(normalizePath(getwd(), winslash = "/", mustWork = TRUE))
  script_path <- sub("^--file=", "", file_arg[[1]])
  script_path <- gsub("~+~", " ", script_path, fixed = TRUE)   # decode OneDrive-style space encoding
  dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))
}
.this_dir <- resolve_script_dir()
.utils <- file.path(.this_dir, "prs_risk_utils.R")
if (!file.exists(.utils)) stop("Required helper not found: ", .utils,
                               "\n  Resolved script directory: ", .this_dir,
                               "\n  prs_risk_utils.R must sit beside this script.")
source(.utils)

# save_plot_cairo() is provided by prs_risk_utils.R (sourced above) with a
# headless grDevices fallback.

# --------------------------- COMMAND-LINE ARGS --------------------------------
option_list <- list(
  # Inputs
  make_option("--models_file", type="character", default=NULL,
              help="Path to fitted_models.rds from 01_fit_models.R [REQUIRED]"),
  make_option("--data_file", type="character", default=NULL,
              help="Path to data_processed.rds from 01_fit_models.R [REQUIRED]"),
  
  # NEW: use supported times from 01*
  make_option("--use_supported_times", action="store_true", default=FALSE,
              help="Use metadata$time_support$eval_time_points from 01*. Falls back to requested_valid_time_points, then t_pass_auto, then --time_points/median."),
  
  # (Fallback) explicit times
  make_option("--time_points", type="character", default=NULL,
              help="Comma-separated time points in years (e.g., '5,10'). If NULL and not using supported times, uses median follow-up."),
  
  # NEW: guard thresholds (applied to chosen times before analysis)
  make_option("--min_events", type="integer", default=10,
              help="Minimum cumulative events by t to keep t [default: 10]"),
  make_option("--min_at_risk", type="integer", default=50,
              help="Minimum at-risk at t to keep t [default: 50]"),
  make_option("--ghat_threshold", type="numeric", default=0.02,
              help="Minimum G-hat (censoring survival) at t to keep t [default: 0.02]"),
  
  # Risk thresholds (categorical NRI)
  make_option("--risk_thresholds", type="character", default="0.1,0.2",
              help="Comma-separated thresholds for categorical NRI (e.g., '0.1,0.2') [default: 0.1,0.2]"),
  
  # Model choices
  make_option("--reference_model", type="character", default="base",
              help="Reference model: 'base' or 'clinical' [default: base]"),
  make_option("--comparison_model", type="character", default="minimal",
              help="Comparison model: 'minimal' or 'full' [default: minimal]"),
  
  # Analysis options
  make_option("--bootstrap", action="store_true", default=FALSE,
              help="Bootstrap CIs for IDI/NRI [default: FALSE]"),
  make_option("--n_bootstrap", type="integer", default=500,
              help="Bootstrap replicates [default: 500]"),
  make_option("--conf_level", type="numeric", default=0.95,
              help="Confidence level [default: 0.95]"),
  
  # Output / plotting
  make_option("--outdir", type="character", default="output",
              help="Output directory [default: output]"),
  make_option("--prefix", type="character", default="",
              help="Prefix for output files [default: auto]"),
  make_option("--plot_width", type="numeric", default=10,
              help="Plot width [inches] [default: 10]"),
  make_option("--plot_height", type="numeric", default=8,
              help="Plot height [inches] [default: 8]"),
  make_option("--verbose", action="store_true", default=TRUE,
              help="Print detailed output [default: TRUE]"),
  make_option("--make_plots", action="store_true", default=FALSE,
              help="Generate reclassification plots [default: FALSE]"),
  make_option("--plot_formats", type="character", default="pdf",
              help="Comma-separated formats to save (pdf,png) [default: pdf]")
)

opt <- parse_args(OptionParser(
  option_list=option_list,
  description="\nCalculate reclassification metrics (NRI, IDI) from fitted models.",
  epilogue="
Examples:
  Rscript 04_reclassification.R --models_file output/fitted_models.rds --data_file output/data_processed.rds
  Rscript 04_reclassification.R --models_file output/fitted_models.rds --data_file output/data_processed.rds --use_supported_times
  Rscript 04_reclassification.R --models_file output/fitted_models.rds --data_file output/data_processed.rds --time_points '5,10,15'
"
))

# ------------------------------ VALIDATE I/O ----------------------------------
if (is.null(opt$models_file)) stop("--models_file is required")
if (is.null(opt$data_file))   stop("--data_file is required")
if (!file.exists(opt$models_file)) stop("Models file not found: ", opt$models_file)
if (!file.exists(opt$data_file))   stop("Data file not found: ", opt$data_file)
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

cat("==============================================================================\n")
cat("RECLASSIFICATION ANALYSIS\n")
cat("==============================================================================\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("==============================================================================\n\n")

# ---------------------------- LOAD DATA/MODELS --------------------------------
cat("LOADING DATA\n")
cat("----------------------------------------------------------------------\n")
results <- readRDS(opt$models_file)
DT      <- readRDS(opt$data_file)
cat("✓ Loaded fitted models\n")
cat("✓ Loaded processed data (N=", nrow(DT), ")\n\n", sep="")

metadata <- results$metadata

# Prefix
if (nchar(opt$prefix) > 0) {
  prefix <- paste0(opt$prefix, "_04_reclassification_")
} else {
  base_prefix_parts <- c()
  if (!is.null(metadata$options$cohort)     && metadata$options$cohort     != "") base_prefix_parts <- c(base_prefix_parts, metadata$options$cohort)
  if (!is.null(metadata$options$ancestry)   && metadata$options$ancestry   != "") base_prefix_parts <- c(base_prefix_parts, metadata$options$ancestry)
  if (!is.null(metadata$options$pheno_name) && metadata$options$pheno_name != "") base_prefix_parts <- c(base_prefix_parts, metadata$options$pheno_name)
  prefix <- if (length(base_prefix_parts) > 0) paste0(paste(base_prefix_parts, collapse="_"), "_04_reclassification_") else "04_reclassification_"
}
cat("Output prefix:", prefix, "\n\n")

# Logging
log_file <- file.path(opt$outdir, paste0(prefix, "log.txt"))
log_con  <- file(log_file, open="wt")
sink(log_con, type="output", split=opt$verbose)
sink(log_con, type="message")

# ----------------------------- SELECT MODELS ----------------------------------
# Model selection accepts ANY model name present in the fitted-models object, so
# manuscript comparisons such as M5 vs M2 (--reference_model=Outcome
# --comparison_model=Outcome_Progression) and M6 vs M7 are expressible. The
# previous switch() hard-coded base/clinical/minimal/full and could not reach
# the multi-PRS models at all.
available_models <- names(results$models)
model_map <- results$metadata$model_map   # may be NULL for single-PRS runs

resolve_model <- function(name, what) {
  if (is.null(name) || !nzchar(name)) stop("--", what, " is required")
  if (name %in% available_models) return(results$models[[name]])
  # tolerate manuscript ids (M0-M7) and case-insensitive names
  if (!is.null(model_map) && name %in% model_map$model_id) {
    nm <- model_map[model_id == name, model][1]
    if (!is.na(nm) && nm %in% available_models) return(results$models[[nm]])
    if (identical(nm, "Base") && "base" %in% available_models) return(results$models$base)
  }
  hit <- available_models[tolower(available_models) == tolower(name)]
  if (length(hit) == 1L) return(results$models[[hit]])
  stop("Model '", name, "' not found for --", what, ".\n  Available models: ",
       paste(available_models, collapse=", "),
       if (!is.null(model_map)) paste0("\n  Available ids: ",
         paste(stats::na.omit(model_map$model_id), collapse=", ")) else "")
}

ref_model  <- resolve_model(opt$reference_model,  "reference_model")
comp_model <- resolve_model(opt$comparison_model, "comparison_model")

# Manuscript ids for the two chosen models, for labelling the outputs.
id_for <- function(name) {
  if (is.null(model_map)) return(NA_character_)
  if (name %in% model_map$model_id) return(name)
  hit <- model_map[tolower(model) == tolower(name), model_id]
  if (length(hit)) hit[1] else NA_character_
}
ref_model_id  <- id_for(opt$reference_model)
comp_model_id <- id_for(opt$comparison_model)

cat("Configuration:\n")
cat("  Reference model:", opt$reference_model,
    if (!is.na(ref_model_id))  paste0(" (", ref_model_id, ")")  else "", "\n")
cat("  Comparison model:", opt$comparison_model,
    if (!is.na(comp_model_id)) paste0(" (", comp_model_id, ")") else "", "\n")
cat("  N observations:", format(nrow(DT), big.mark=","), "\n")
cat("  N events:", format(sum(DT$PHENO), big.mark=","), sprintf(" (%.1f%%)\n", 100*mean(DT$PHENO)))

# ------------------------------- HELPERS --------------------------------------
# Predicted absolute risk now comes from the shared helper (prs_risk_utils.R),
# which indexes the baseline hazard on the model's own time scale and returns
# risk conditional on entry. The previous local implementation indexed a
# start-stop model's baseline (age scale) with a follow-up-scale horizon, and
# used a nearest-neighbour time lookup that could select a time AFTER the
# requested horizon.
calculate_predicted_risk <- function(fit, newdata, time_point) {
  predict_risk_at_horizon(fit, newdata, time_point, start_col = "START")
}

# NOTE ON CENSORING (applies to all three estimators below).
# These are binary-outcome statistics. Applied to censored survival data with the
# EVENTUAL status they answer the wrong question: someone whose T2 occurs after
# the evaluated horizon is counted as an event at that horizon, and someone
# censored before the horizon is counted as a non-event. Both are wrong.
# `event` must therefore be the status AT THE HORIZON (D) and `w` the IPCW weight
# from ipcw_weights(); individuals censored before the horizon carry weight 0 and
# are excluded. With no censoring before the horizon all weights are equal and
# these reduce exactly to the classical uncensored formulas.
.wmean <- function(x, w) { ok <- is.finite(x) & is.finite(w); if (!any(ok) || sum(w[ok]) == 0) NA_real_ else sum(x[ok]*w[ok])/sum(w[ok]) }

# boot_time / boot_status / boot_horizon: when supplied, each bootstrap replicate
# RECOMPUTES the censoring KM and the IPCW weights on its own resampled rows, so
# the interval reflects uncertainty in the censoring distribution too. Resampling
# fixed weights (the previous behaviour) treats G(t) as if it were known exactly
# and understates the variance.
calculate_IDI <- function(pred_base, pred_new, event, bootstrap_ci=FALSE, n_boot=500, conf_level=0.95,
                          w=NULL, boot_time=NULL, boot_status=NULL, boot_horizon=NULL,
                          boot_pred_base=NULL, boot_pred_new=NULL) {
  if (is.null(w)) w <- rep(1, length(event))
  # Bootstrap over the FULL sample when the full vectors are supplied: resample all
  # rows, recompute the censoring KM + IPCW on the resample, then drop unknown-
  # status/non-finite. Resampling the already-keep-filtered subset (which had
  # dropped everyone censored before the horizon) cannot reproduce the censoring
  # process. Returns the kept rows' base/new predictions, status and weights.
  .boot_draw <- function() {
    if (is.null(boot_time) || is.null(boot_status) || is.null(boot_horizon) ||
        is.null(boot_pred_base) || is.null(boot_pred_new)) {
      idx <- sample(length(event), replace=TRUE)
      return(list(pb=pred_base[idx], pn=pred_new[idx], D=event[idx], w=w[idx]))
    }
    idx <- sample(length(boot_time), replace=TRUE)
    ip <- tryCatch(ipcw_weights(boot_time[idx], boot_status[idx], boot_horizon),
                   error = function(e) NULL)
    pb <- boot_pred_base[idx]; pn <- boot_pred_new[idx]
    # IPCW failed on this resample: DROP the replicate (return NULL) rather than
    # substitute all-non-events, which would bias IDI/NRI toward zero.
    if (is.null(ip)) return(NULL)
    k <- ip$keep & is.finite(pb) & is.finite(pn)
    list(pb=pb[k], pn=pn[k], D=ip$D[k], w=ip$w[k])
  }
  we <- w * (event == 1); wn <- w * (event == 0)
  IS_base_events <- .wmean(pred_base, we); IS_new_events <- .wmean(pred_new, we)
  IS_base_nonev  <- .wmean(pred_base, wn); IS_new_nonev  <- .wmean(pred_new, wn)
  IDI <- (IS_new_events - IS_base_events) - (IS_new_nonev - IS_base_nonev)
  denom <- (IS_base_events - IS_base_nonev)
  rel_IDI <- IDI / denom
  result <- data.table(IDI=IDI, rel_IDI=rel_IDI,
                       IS_base=IS_base_events, IS_new=IS_new_events,
                       IP_base=IS_base_nonev,  IP_new=IS_new_nonev)
  if (bootstrap_ci) {
    cat("  Calculating bootstrap confidence intervals (", n_boot, " samples)...\n", sep="")
    boot_IDI <- boot_rel <- rep(NA_real_, n_boot)
    for (i in seq_len(n_boot)) {
      dr <- .boot_draw()
      if (is.null(dr)) next            # IPCW failed on this resample; drop it
      pb <- dr$pb; pn <- dr$pn; e <- dr$D; ww <- dr$w
      wwe <- ww * (e == 1); wwn <- ww * (e == 0)
      ISb_e <- .wmean(pb, wwe); ISn_e <- .wmean(pn, wwe)
      ISb_n <- .wmean(pb, wwn); ISn_n <- .wmean(pn, wwn)
      boot_IDI[i] <- (ISn_e - ISb_e) - (ISn_n - ISb_n)
      dsb <- (ISb_e - ISb_n); boot_rel[i] <- if (is.finite(dsb) && dsb != 0) boot_IDI[i]/dsb else NA
    }
    result$IDI_n_boot_success <- sum(!is.na(boot_IDI))
    if (result$IDI_n_boot_success < n_boot)
      cat(sprintf("  Note: %d/%d IDI bootstrap replicates dropped (IPCW failure on resample).\n",
                  n_boot - result$IDI_n_boot_success, n_boot))
    alpha <- 1 - conf_level
    result$IDI_lower <- quantile(boot_IDI, alpha/2, na.rm=TRUE)
    result$IDI_upper <- quantile(boot_IDI, 1-alpha/2, na.rm=TRUE)
    result$IDI_se    <- sd(boot_IDI, na.rm=TRUE)
    result$IDI_z     <- result$IDI / result$IDI_se
    result$IDI_p     <- 2 * pnorm(-abs(result$IDI_z))
    result$rel_IDI_lower <- quantile(boot_rel, alpha/2, na.rm=TRUE)
    result$rel_IDI_upper <- quantile(boot_rel, 1-alpha/2, na.rm=TRUE)
    result$rel_IDI_se    <- sd(boot_rel, na.rm=TRUE)
    result$rel_IDI_z     <- result$rel_IDI / result$rel_IDI_se
    result$rel_IDI_p     <- 2 * pnorm(-abs(result$rel_IDI_z))
  }
  result
}

calculate_NRI_continuous <- function(pred_base, pred_new, event, bootstrap_ci=FALSE, n_boot=500, conf_level=0.95,
                                     w=NULL, boot_time=NULL, boot_status=NULL, boot_horizon=NULL,
                                     boot_pred_base=NULL, boot_pred_new=NULL) {
  if (is.null(w)) w <- rep(1, length(event))
  # See calculate_IDI(): bootstrap over the FULL sample, recomputing IPCW per
  # replicate, when the full vectors are supplied.
  .boot_draw <- function() {
    if (is.null(boot_time) || is.null(boot_status) || is.null(boot_horizon) ||
        is.null(boot_pred_base) || is.null(boot_pred_new)) {
      idx <- sample(length(event), replace=TRUE)
      return(list(pb=pred_base[idx], pn=pred_new[idx], D=event[idx], w=w[idx]))
    }
    idx <- sample(length(boot_time), replace=TRUE)
    ip <- tryCatch(ipcw_weights(boot_time[idx], boot_status[idx], boot_horizon),
                   error = function(e) NULL)
    pb <- boot_pred_base[idx]; pn <- boot_pred_new[idx]
    # IPCW failed on this resample: DROP the replicate (return NULL) rather than
    # substitute all-non-events, which would bias IDI/NRI toward zero.
    if (is.null(ip)) return(NULL)
    k <- ip$keep & is.finite(pb) & is.finite(pn)
    list(pb=pb[k], pn=pn[k], D=ip$D[k], w=ip$w[k])
  }
  # Counts become IPCW-weighted sums; denominators the weighted group totals.
  .cnt <- function(sel, grp, ww) sum(ww[sel & grp], na.rm=TRUE)
  ev <- event == 1; nev <- event == 0
  up <- pred_new > pred_base; dn <- pred_new < pred_base
  up_e <- .cnt(up, ev, w); down_e <- .cnt(dn, ev, w); n_e <- sum(w[ev], na.rm=TRUE)
  up_n <- .cnt(up, nev, w); down_n <- .cnt(dn, nev, w); n_n <- sum(w[nev], na.rm=TRUE)
  NRI_e <- (up_e - down_e)/n_e; NRI_n <- (down_n - up_n)/n_n; NRI <- NRI_e + NRI_n
  res <- data.table(NRI=NRI, NRI_events=NRI_e, NRI_nonevents=NRI_n,
                    up_events=up_e, down_events=down_e, up_nonevents=up_n, down_nonevents=down_n)
  if (bootstrap_ci) {
    cat("  Calculating bootstrap confidence intervals (", n_boot, " samples)...\n", sep="")
    boot <- rep(NA_real_, n_boot)
    for (i in seq_len(n_boot)) {
      dr <- .boot_draw()
      if (is.null(dr)) next            # IPCW failed on this resample; drop it
      pb <- dr$pb; pn <- dr$pn; e <- dr$D; ww <- dr$w
      eev <- e == 1; enev <- e == 0; u <- pn > pb; d <- pn < pb
      ue <- .cnt(u, eev, ww); de <- .cnt(d, eev, ww); ne <- sum(ww[eev], na.rm=TRUE)
      un <- .cnt(u, enev, ww); dnn <- .cnt(d, enev, ww); nn <- sum(ww[enev], na.rm=TRUE)
      boot[i] <- (ue - de)/ne + (dnn - un)/nn
    }
    res$NRI_n_boot_success <- sum(!is.na(boot))
    if (res$NRI_n_boot_success < n_boot)
      cat(sprintf("  Note: %d/%d NRI bootstrap replicates dropped (IPCW failure on resample).\n",
                  n_boot - res$NRI_n_boot_success, n_boot))
    alpha <- 1 - conf_level
    res$NRI_lower <- quantile(boot, alpha/2, na.rm=TRUE)
    res$NRI_upper <- quantile(boot, 1-alpha/2, na.rm=TRUE)
    res$NRI_se    <- sd(boot, na.rm=TRUE)
    res$NRI_z     <- res$NRI / res$NRI_se
    res$NRI_p     <- 2 * pnorm(-abs(res$NRI_z))
  }
  res
}

calculate_NRI_categorical <- function(pred_base, pred_new, event, thresholds, w=NULL) {
  if (is.null(w)) w <- rep(1, length(event))
  thr <- c(-Inf, thresholds, Inf)
  # Ensure all categories are present as levels, even if empty
  all_labels <- 1:(length(thr) - 1)

  cb  <- factor(cut(pred_base, thr, labels=FALSE), levels=all_labels)
  cn  <- factor(cut(pred_new, thr, labels=FALSE), levels=all_labels)

  # IPCW-weighted reclassification tables: cell = sum of weights, not raw count.
  # xtabs keeps the full level grid so upper/lower.tri stay aligned.
  .wtab <- function(sel) {
    d <- data.frame(cb = cb[sel], cn = cn[sel], w = w[sel])
    as.matrix(stats::xtabs(w ~ cb + cn, data = d, drop.unused.levels = FALSE))
  }
  T_e <- .wtab(event == 1)
  T_n <- .wtab(event == 0)
  # T[i, j] = count moving from reference category i to new category j.
  # UPWARD reclassification is j > i, i.e. the UPPER triangle.
  # These two were previously swapped, which inverted the sign of the
  # categorical NRI: a model that improved reclassification was reported as
  # harming it. Verified against a table with known movements.
  up_e   <- sum(T_e[upper.tri(T_e)]);   down_e <- sum(T_e[lower.tri(T_e)])
  up_n   <- sum(T_n[upper.tri(T_n)]);   down_n <- sum(T_n[lower.tri(T_n)])
  n_e <- sum(w[event==1], na.rm=TRUE); n_n <- sum(w[event==0], na.rm=TRUE)
  data <- list(
    NRI       = (up_e - down_e)/n_e + (down_n - up_n)/n_n,
    NRI_events= (up_e - down_e)/n_e,
    NRI_nonevents=(down_n - up_n)/n_n,
    reclass_table_events=T_e, reclass_table_nonevents=T_n, thresholds=thresholds
  )
  data
}

# --------------------------- TIME POINT SELECTION ------------------------------
cat("\nTIME POINT SELECTION\n")
cat("----------------------------------------------------------------------\n")

parse_times_cli <- function(x) {
  if (is.null(x) || nchar(x)==0) return(numeric(0))
  t <- suppressWarnings(as.numeric(trimws(strsplit(x, ",")[[1]])))
  sort(unique(t[is.finite(t) & t > 0]))
}

candidate_times <- numeric(0)
if (isTRUE(opt$use_supported_times)) {
  ts <- metadata$time_support
  if (!is.null(ts$eval_time_points) && length(ts$eval_time_points) > 0) {
    candidate_times <- as.numeric(ts$eval_time_points)
    cat("Using eval_time_points from 01*: ", paste(candidate_times, collapse=", "), "\n", sep="")
  } else if (!is.null(ts$requested_valid_time_points) && length(ts$requested_valid_time_points) > 0) {
    candidate_times <- as.numeric(ts$requested_valid_time_points)
    cat("Using requested_valid_time_points from 01*: ", paste(candidate_times, collapse=", "), "\n", sep="")
  } else if (!is.null(ts$t_pass_auto) && length(ts$t_pass_auto) > 0) {
    candidate_times <- as.numeric(ts$t_pass_auto)
    cat("Using t_pass_auto from 01*: ", paste(candidate_times, collapse=", "), "\n", sep="")
  } else {
    cat("No supported times in metadata; falling back.\n")
  }
}

if (length(candidate_times) == 0) {
  if (!is.null(opt$time_points)) {
    candidate_times <- parse_times_cli(opt$time_points)
    cat("Using --time_points: ", paste(candidate_times, collapse=", "), "\n", sep="")
  } else {
    candidate_times <- median(DT$TIME)
    cat("Using median follow-up time: ", sprintf("%.2f", candidate_times), "\n", sep="")
  }
}

max_time <- max(DT$TIME, na.rm=TRUE)
candidate_times <- sort(unique(as.numeric(candidate_times)))
candidate_times <- candidate_times[candidate_times > 0 & candidate_times <= max_time]
stopifnot(length(candidate_times) > 0)

# --------- Apply guards (G-hat, events, at-risk) & save audit CSV -------------
cens_surv <- survfit(Surv(DT$TIME, 1 - DT$PHENO) ~ 1)
audit <- data.table(time_years = candidate_times)
audit[, `:=`(
  ghat          = vapply(time_years, function(t) { out <- tryCatch(summary(cens_surv, times=t)$surv, error=function(e) NA_real_); as.numeric(out) }, numeric(1)),
  n_events_by_t = vapply(time_years, function(t) sum(DT$PHENO==1 & DT$TIME <= t), integer(1)),
  n_at_risk     = vapply(time_years, function(t) sum(DT$TIME >= t), integer(1))
)]
audit[, status := fifelse(!is.finite(ghat), "No estimate",
                          fifelse(ghat < opt$ghat_threshold, sprintf("FAIL (<%g)", opt$ghat_threshold),
                                  fifelse(n_at_risk < opt$min_at_risk, sprintf("FAIL (at_risk<%d)", opt$min_at_risk),
                                          fifelse(n_events_by_t < opt$min_events, sprintf("FAIL (events<%d)", opt$min_events), "PASS"))))]

time_points <- audit[status=="PASS", time_years]
cat("Times after guards: ", if (length(time_points)) paste(time_points, collapse=", ") else "NONE", "\n\n", sep="")
if (length(time_points)==0) stop("No time points passed the guards.")

audit_file <- file.path(opt$outdir, paste0(prefix, "reclassification_times_ghat_audit.csv"))
fwrite(audit, audit_file)
cat("✓ Saved time audit: ", basename(audit_file), "\n\n", sep="")

# ---------------------------- RISK THRESHOLDS ---------------------------------
risk_thresholds <- sort(as.numeric(trimws(strsplit(opt$risk_thresholds, ",")[[1]])))

# --------------------------- RESULT COLLECTORS --------------------------------
all_IDI_results <- list()
all_NRI_cont_results <- list()
all_NRI_cat_results <- list()
all_reclass_tables <- list()
all_risk_movements <- list()
all_plots <- list()

# ------------------------------ MAIN LOOP -------------------------------------
for (tp_idx in seq_along(time_points)) {
  time_point <- time_points[tp_idx]
  
  cat("==============================================================================\n")
  cat("ANALYSIS FOR TIME POINT: ", sprintf("%.2f years", time_point), "\n")
  cat("==============================================================================\n\n")
  
  cat("CALCULATING PREDICTED RISKS\n")
  cat("----------------------------------------------------------------------\n")
  cat("Reference model (", opt$reference_model, ")...\n", sep="")
  pred_base <- calculate_predicted_risk(ref_model, DT, time_point)
  cat("Comparison model (", opt$comparison_model, ")...\n", sep="")
  pred_new  <- calculate_predicted_risk(comp_model, DT, time_point)
  
  cat("  Mean risk (reference):", sprintf("%.4f\n", mean(pred_base)))
  cat("  Mean risk (comparison):", sprintf("%.4f\n", mean(pred_new)))
  cat("  Risk range (reference):", sprintf("%.4f - %.4f\n", min(pred_base), max(pred_base)))
  cat("  Risk range (comparison):", sprintf("%.4f - %.4f\n", min(pred_new), max(pred_new)))
  cat("\n")

  # -------------------- CENSORING-AWARE OUTCOME AT THIS HORIZON ---------------
  # Status at `time_point`, with IPCW weights. Using DT$PHENO (the EVENTUAL
  # status) here would count a T2 occurring after the horizon as an event at the
  # horizon, and treat someone censored before the horizon as a non-event.
  ipcw <- ipcw_weights(DT$TIME, DT$PHENO, time_point)
  D_h  <- ipcw$D
  w_h  <- ipcw$w
  # Individuals whose risk is NA (horizon outside baseline support) cannot be
  # classified. Without this they would stay in the weighted denominators while
  # contributing NA movements, biasing every rate downward.
  keep <- ipcw$keep & is.finite(pred_base) & is.finite(pred_new)
  n_nonfinite <- sum(ipcw$keep & !(is.finite(pred_base) & is.finite(pred_new)))
  if (n_nonfinite > 0)
    cat(sprintf("  %d individual(s) additionally excluded: predicted risk is NA at this horizon.\n",
                n_nonfinite))
  cat("CENSORING-AWARE OUTCOME AT HORIZON\n")
  cat("----------------------------------------------------------------------\n")
  cat("  Horizon:", sprintf("%.2f years", time_point), "\n")
  cat("  Censored before horizon (dropped, status unknowable):", ipcw$n_dropped, "\n")
  cat("  Censoring survival G(t):", sprintf("%.4f", ipcw$G_at_horizon), "\n")
  cat("  Contributing individuals:", sum(keep), "of", nrow(DT), "\n")
  cat("  Events by horizon (weighted):", sprintf("%.1f", sum(w_h[D_h == 1])),
      " | event-free (weighted):", sprintf("%.1f", sum(w_h[D_h == 0 & keep])), "\n")
  if (sum(D_h[keep] == 1) < 10)
    cat("  ⚠ Fewer than 10 events by this horizon; reclassification is unstable here.\n")
  cat("\n")
  # Restrict to individuals whose status at the horizon is known.
  pb_k <- pred_base[keep]; pn_k <- pred_new[keep]
  D_k  <- D_h[keep];       w_k  <- w_h[keep]
  ipcw_meta <- data.table(time_point = time_point,
                          n_contributing = sum(keep),
                          n_dropped_censored_before_horizon = ipcw$n_dropped,
                          G_at_horizon = ipcw$G_at_horizon,
                          weighted_events = sum(w_h[D_h == 1]),
                          weighted_nonevents = sum(w_h[D_h == 0 & keep]))

  # IDI
  cat("INTEGRATED DISCRIMINATION IMPROVEMENT (IDI, IPCW at horizon)\n")
  cat("----------------------------------------------------------------------\n")
  # Point estimate on the keep-subset; bootstrap resamples the FULL sample
  # (pred_base/pred_new over all rows, full TIME/PHENO) and recomputes IPCW.
  IDI_result <- calculate_IDI(pb_k, pn_k, D_k, bootstrap_ci=opt$bootstrap,
                              n_boot=opt$n_bootstrap, conf_level=opt$conf_level, w=w_k,
                              boot_time=DT$TIME, boot_status=DT$PHENO, boot_horizon=time_point,
                              boot_pred_base=pred_base, boot_pred_new=pred_new)
  cat("IDI:", sprintf("%.5f", IDI_result$IDI))
  if (opt$bootstrap) cat(sprintf(" (95%% CI: %.5f - %.5f, p=%s)", IDI_result$IDI_lower, IDI_result$IDI_upper, format.pval(IDI_result$IDI_p, digits=3)))
  cat("\n  Relative IDI:", sprintf("%.3f%%", 100*IDI_result$rel_IDI), "\n\n")
  IDI_result[, time_point := time_point]
  IDI_result <- cbind(IDI_result, ipcw_meta[, !"time_point"])
  all_IDI_results[[tp_idx]] <- IDI_result

  # Continuous NRI
  cat("NET RECLASSIFICATION IMPROVEMENT (CONTINUOUS, IPCW at horizon)\n")
  cat("----------------------------------------------------------------------\n")
  NRI_cont_result <- calculate_NRI_continuous(pb_k, pn_k, D_k,
                                              bootstrap_ci=opt$bootstrap,
                                              n_boot=opt$n_bootstrap,
                                              conf_level=opt$conf_level, w=w_k,
                                              boot_time=DT$TIME, boot_status=DT$PHENO,
                                              boot_horizon=time_point,
                                              boot_pred_base=pred_base, boot_pred_new=pred_new)
  cat("  Total NRI:", sprintf("%.4f", NRI_cont_result$NRI), "\n")
  cat("  NRI (events):", sprintf("%.4f", NRI_cont_result$NRI_events), "\n")
  cat("  NRI (non-events):", sprintf("%.4f", NRI_cont_result$NRI_nonevents), "\n\n")
  NRI_cont_result[, time_point := time_point]
  NRI_cont_result <- cbind(NRI_cont_result, ipcw_meta[, !"time_point"])
  all_NRI_cont_results[[tp_idx]] <- NRI_cont_result

  # Categorical NRI
  cat("NET RECLASSIFICATION IMPROVEMENT (CATEGORICAL, IPCW at horizon)\n")
  cat("----------------------------------------------------------------------\n")
  cat("Risk thresholds:", paste(sprintf("%.1f%%", 100*risk_thresholds), collapse=", "), "\n\n")
  NRI_cat_result <- calculate_NRI_categorical(pb_k, pn_k, D_k, thresholds=risk_thresholds, w=w_k)
  cat("  Total NRI:", sprintf("%.4f\n", NRI_cat_result$NRI))
  cat("  NRI (events):", sprintf("%.4f\n", NRI_cat_result$NRI_events))
  cat("  NRI (non-events):", sprintf("%.4f\n", NRI_cat_result$NRI_nonevents), "\n\n")
  
  reclass_events_dt <- as.data.table(NRI_cat_result$reclass_table_events); setnames(reclass_events_dt, c("Reference_Category","New_Category","Count"))
  reclass_events_dt[, `:=`(event_status="Event", time_point=time_point)]
  reclass_nonevents_dt <- as.data.table(NRI_cat_result$reclass_table_nonevents); setnames(reclass_nonevents_dt, c("Reference_Category","New_Category","Count"))
  reclass_nonevents_dt[, `:=`(event_status="No Event", time_point=time_point)]
  all_NRI_cat_results[[tp_idx]] <- data.table(time_point=time_point,
                                              thresholds=paste(risk_thresholds, collapse=","),
                                              NRI=NRI_cat_result$NRI,
                                              NRI_events=NRI_cat_result$NRI_events,
                                              NRI_nonevents=NRI_cat_result$NRI_nonevents)
  all_reclass_tables[[tp_idx]] <- rbind(reclass_events_dt, reclass_nonevents_dt)
  
  # Risk movement
  cat("CREATING RISK MOVEMENT DATA\n")
  cat("----------------------------------------------------------------------\n")
  # `event` is the EVENTUAL status (kept for traceability); `event_at_horizon` and
  # `ipcw_weight` are what the IDI/NRI above actually used. Rows with weight 0 were
  # censored before the horizon and contributed to neither.
  risk_movement <- data.table(
    ID = seq_len(nrow(DT)),
    risk_base = pred_base, risk_new = pred_new, event = DT$PHENO,
    event_at_horizon = D_h, ipcw_weight = w_h, contributed = keep,
    risk_change = pred_new - pred_base, time_point = time_point
  )
  risk_movement[, movement := ifelse(risk_change > 0, "Up", ifelse(risk_change < 0, "Down", "No change"))]
  risk_movement[, event_label := ifelse(event_at_horizon == 1, "Event", "No Event")]
  print(risk_movement[, .N, by=.(event_label, movement)]); cat("\n")
  all_risk_movements[[tp_idx]] <- risk_movement
  
  # Plots (opt-in)
  if (opt$make_plots) {
    p1 <- ggplot(risk_movement, aes(x=risk_base, y=risk_new, color=event_label)) +
      geom_abline(slope=1, intercept=0, linetype="dashed", color="gray50") +
      geom_point(alpha=0.5, size=1.5) +
      scale_color_manual(values=c("Event"="red", "No Event"="blue")) +
      labs(x=paste0("Predicted Risk (", opt$reference_model, ")"),
           y=paste0("Predicted Risk (", opt$comparison_model, ")"),
           color="", title=paste0("Risk Reclassification at ", sprintf("%.1f", time_point), " Years")) +
      theme_bw() + theme(legend.position="top")
    p2 <- ggplot(risk_movement, aes(x=risk_change, fill=event_label)) +
      geom_histogram(bins=50, alpha=0.7, position="identity") +
      geom_vline(xintercept=0, linetype="dashed", color="black") +
      scale_fill_manual(values=c("Event"="red", "No Event"="blue")) +
      labs(x="Change in Predicted Risk", y="Count", fill="", title="Distribution of Risk Changes") +
      theme_bw() + theme(legend.position="top")
    p3 <- ggplot(risk_movement, aes(x=event_label, y=risk_change, fill=event_label)) +
      geom_boxplot(alpha=0.7) +
      geom_hline(yintercept=0, linetype="dashed", color="black") +
      scale_fill_manual(values=c("Event"="red", "No Event"="blue")) +
      labs(x="", y="Change in Predicted Risk", title="Risk Change by Event Status") +
      theme_bw() + theme(legend.position="none")
    combined_plot <- ggpubr::ggarrange(p1, p2, p3, ncol=2, nrow=2, layout_matrix=rbind(c(1,1), c(2,3)))
    plot_file_tp <- file.path(opt$outdir, paste0(prefix, "reclassification_plot_", sprintf("%.1f", time_point), "yr.pdf"))
    plot_formats_parsed <- trimws(strsplit(opt$plot_formats, ",")[[1]])
    save_plot_cairo(combined_plot, plot_file_tp, width=opt$plot_width, height=opt$plot_height, formats=plot_formats_parsed)
    cat("✓ Saved:", basename(plot_file_tp), "\n\n")
    all_plots[[tp_idx]] <- combined_plot
  }
}

# ---------------------------- SAVE COMBINED OUTS ------------------------------
cat("==============================================================================\n")
cat("SAVING COMBINED RESULTS\n")
cat("==============================================================================\n\n")

# Annotate every table with both the internal model names and the manuscript ids
# (e.g. "M5 vs M2"), so results are directly citable.
annotate_comparison <- function(dt) {
  dt[, `:=`(
    comparison        = paste(opt$comparison_model, "vs", opt$reference_model),
    comparison_id     = if (!is.na(comp_model_id) && !is.na(ref_model_id))
                          paste(comp_model_id, "vs", ref_model_id) else NA_character_,
    model_comparison  = opt$comparison_model,
    model_reference   = opt$reference_model,
    model_id_comparison = comp_model_id,
    model_id_reference  = ref_model_id
  )]
  setcolorder(dt, c("comparison","comparison_id","time_point"))
  dt[]
}

IDI_combined <- annotate_comparison(rbindlist(all_IDI_results))
fwrite(IDI_combined, file.path(opt$outdir, paste0(prefix, "idi_results.csv")));           cat("✓ Saved IDI results\n")

NRI_cont_combined <- annotate_comparison(rbindlist(all_NRI_cont_results))
fwrite(NRI_cont_combined, file.path(opt$outdir, paste0(prefix, "nri_continuous.csv")));   cat("✓ Saved continuous NRI\n")

NRI_cat_combined <- annotate_comparison(rbindlist(all_NRI_cat_results))
fwrite(NRI_cat_combined, file.path(opt$outdir, paste0(prefix, "nri_categorical.csv")));   cat("✓ Saved categorical NRI\n")

reclass_combined <- rbindlist(all_reclass_tables)
fwrite(reclass_combined, file.path(opt$outdir, paste0(prefix, "reclassification_table.csv")));  cat("✓ Saved reclassification tables\n")

risk_movement_combined <- rbindlist(all_risk_movements)
fwrite(risk_movement_combined, file.path(opt$outdir, paste0(prefix, "risk_movement.csv")));     cat("✓ Saved risk movement data\n")

if (opt$make_plots && length(time_points) == 1 && length(all_plots)) {
  plot_file <- file.path(opt$outdir, paste0(prefix, "reclassification_plot.pdf"))
  plot_formats_parsed <- trimws(strsplit(opt$plot_formats, ",")[[1]])
  save_plot_cairo(all_plots[[1]], plot_file, width=opt$plot_width, height=opt$plot_height, formats=plot_formats_parsed)
  cat("✓ Saved:", basename(plot_file), "\n")
}

# -------------------------------- SUMMARY -------------------------------------
cat("\n==============================================================================\n")
cat("SUMMARY\n")
cat("==============================================================================\n\n")
cat("✓ Reclassification analysis completed!\n\n")
cat("Comparison:", opt$comparison_model, "vs", opt$reference_model, "\n")
cat("Time points analyzed:", paste(sprintf("%.2f", time_points), collapse=", "), "years\n\n")

cat("Key Results Summary:\n\n")
for (tp in time_points) {
  cat(sprintf("Time point: %.2f years\n", tp))
  cat("  IDI:", sprintf("%.5f", IDI_combined[time_point==tp]$IDI), "\n")
  cat("  Relative IDI:", sprintf("%.3f%%", 100*IDI_combined[time_point==tp]$rel_IDI), "\n")
  cat("  Continuous NRI:", sprintf("%.4f", NRI_cont_combined[time_point==tp]$NRI), "\n")
  cat("  Categorical NRI:", sprintf("%.4f", NRI_cat_combined[time_point==tp]$NRI), "\n\n")
}
cat("Outputs saved to:", opt$outdir, "\n")
cat("  - idi_results.csv, nri_continuous.csv, nri_categorical.csv\n")
cat("  - reclassification_table.csv, risk_movement.csv\n")
if (opt$make_plots) cat("  - reclassification_plot_*.pdf\n")
cat("\n==============================================================================\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("==============================================================================\n")

sink(type="message"); sink(type="output"); close(log_con)
cat("\nLog file:", log_file, "\n")
