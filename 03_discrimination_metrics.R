#!/usr/bin/env Rscript
# ==============================================================================
# 03_discrimination_metrics.R - Time-Dependent Discrimination & Integrated Metrics
# ==============================================================================
# Purpose:
#   - Load fitted models/data from 01_fit_models.R
#   - Compute time-dependent AUC(t) for each model
#   - Compute per-horizon iAUC(0,t) and IBS(0,t) via trapezoidal integration
#   - (When available) also save single-τ iAUC and IBS returned by Score()
#   - Generate Brier(t) prediction-error curves
#   - Optionally compute AUC(t) using timeROC (for comparison)
#   - Apply G-hat, events, and at-risk guards at each time point
#   - Save plots and a merged summary table
#
# NEW:
#   - --use_supported_times uses eval_time_points saved by 01* (metadata$time_support$eval_time_points)
#     falling back to requested_valid_time_points, then t_pass_auto, and re-applies guards here.
#
# Features:
#   - riskRegression::Score as the primary engine for AUC(t) & Brier(t)
#   - Optional timeROC engine (or both) for AUC(t) sensitivity checks
#   - Per-horizon integration: iAUC(0,t), IBS(0,t) at every valid time
#   - Backward-compatible: captures Score()-returned iAUC / IBS if present
#   - Handles Null model AUC(t) fallback (0.5) when not returned by engine
#   - Influence-function CIs by default; optional bootstrap
#
# Inputs:
#   - fitted_models.rds     (from 01_fit_models.R)
#   - data_processed.rds    (from 01_fit_models.R)
#
# Key options (CLI):
#   --use_supported_times   Use time points saved by 01* (recommended)
#   --time_points           Comma-separated years for evaluation (e.g., '1,2,5,10')
#   --min_events            Minimum cumulative events by t to keep t
#   --min_at_risk           Minimum at-risk at t to keep t
#   --ghat_threshold        Minimum censoring survival G-hat(t) to keep t
#   --auc_engine            'riskRegression' (default), 'timeROC', or 'both'
#   --compute_integrated    Also compute integrated metrics and curves (default TRUE)
#   --tau                   Max horizon for single-τ summaries (defaults to max valid t)
#   --bootstrap, --n_boot   Bootstrap CIs for integrated metrics (optional)
#   --skip_ci               Skip CIs for faster runs
#   --save_plots            Save AUC(t) and Brier(t) figures (default TRUE)
#   --outdir, --prefix      Output path/naming
#
# Outputs (prefix_*):
#   - ghat_censoring.csv                  # G-hat at requested/supported times + pass/fail flags
#   - rr_auc_by_time.csv                  # riskRegression AUC(t)
#   - rr_brier_curve.csv                  # Brier(t) prediction-error curves
#   - rr_iauc_by_time.csv                 # iAUC(0,t) at each valid time
#   - rr_ibs_by_time.csv                  # IBS(0,t) at each valid time
#   - rr_iauc.csv                         # Single-τ iAUC from Score() (if returned)
#   - rr_ibs.csv                          # Single-τ IBS  from Score() (if returned)
#   - discrimination_summary.csv          # Merge of AUC(t), Brier(t), iAUC(0,t), IBS(0,t), G-hat
#   - auc_time_curve.png                  # AUC(t) curves
#   - pe_curve.png                        # Prediction-error (Brier) curves
#   - discrimination_metrics.log          # Verbose run log
#   - tdauc_results.csv                   # (Optional) AUC(t) from timeROC if requested
#
# Author: Ying Wang (yiwang@broadinstitute.org)
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

# Cairo is NOT mandatory — it needs X11/XQuartz and is only used for figures.
# Keeping it out of required_packages lets the metric CSVs run headless.
required_packages <- c(
  "optparse", "data.table", "survival", "splines",
  "riskRegression", "ggplot2"
)
invisible(lapply(required_packages, load_package))

cat("Package versions:\n")
for (p in required_packages) cat(sprintf("  %-14s %s\n", paste0(p, ":"), as.character(utils::packageVersion(p))))

# Save a plot with Cairo if available (nicer, needs X11); otherwise fall back to
# base grDevices so a headless machine still produces figures. Never aborts.
save_plot_cairo <- function(plot_obj, filename, width, height, dpi=300, formats=c("png")) {
  have_cairo <- requireNamespace("Cairo", quietly = TRUE)
  for (fmt in formats) {
    if (!fmt %in% c("png","pdf")) { warning("Unsupported format: ", fmt); next }
    outfile <- sub("\\.(png|pdf)$", paste0(".", fmt), filename)
    ok <- FALSE
    if (have_cairo) ok <- tryCatch({
      if (fmt == "png") Cairo::CairoPNG(filename = outfile, width = width, height = height, units = "in", dpi = dpi)
      else              Cairo::CairoPDF(file = outfile, width = width, height = height)
      print(plot_obj); grDevices::dev.off(); TRUE
    }, error = function(e) { try(grDevices::dev.off(), silent = TRUE); FALSE })
    if (!ok) tryCatch({
      if (fmt == "png") grDevices::png(filename = outfile, width = width, height = height, units = "in", res = dpi)
      else              grDevices::pdf(file = outfile, width = width, height = height)
      print(plot_obj); grDevices::dev.off()
      if (!have_cairo) message("Cairo unavailable; wrote ", basename(outfile), " via grDevices fallback")
    }, error = function(e) warning("Could not save plot ", basename(outfile), ": ", conditionMessage(e)))
  }
}


# ==============================================================================
# COMMAND-LINE ARGUMENTS
# ==============================================================================
option_list <- list(
  make_option("--fitted_models", type="character",
              help="Path to fitted_models.rds from 01_fit_models.R [REQUIRED]"),
  make_option("--data_processed", type="character",
              help="Path to data_processed.rds from 01_fit_models.R [REQUIRED]"),
  make_option("--use_supported_times", action="store_true", default=FALSE,
              help="Use eval_time_points saved by 01* (metadata$time_support$eval_time_points). Overrides --time_points when present."),
  make_option("--time_points", type="character", default="1,2,5,10",
              help="Comma-separated time points for analysis (years) [default: 1,2,5,10]"),
  make_option("--min_events", type="integer", default=10,
              help="Minimum events required at time point [default: 10]"),
  make_option("--min_at_risk", type="integer", default=50,
              help="Minimum at-risk required at time point [default: 50]"),
  make_option("--ghat_threshold", type="numeric", default=0.02,
              help="Minimum G-hat (censoring survival) threshold [default: 0.02]"),
  make_option("--auc_engine", type="character", default="riskRegression",
              help="AUC calculation method: 'riskRegression' (default), 'timeROC', or 'both' [default: riskRegression]"),
  make_option("--score_all_models", action="store_true", default=FALSE,
              help=paste("Score every fitted model including the per-quantile helpers.",
                         "Default scores only the M0-M7 ladder, which keeps the pairwise",
                         "contrast tables interpretable. [default: FALSE]")),
  make_option("--contrast_pairs", type="character", default="M4:M1,M5:M2,M7:M6",
              help=paste("Comma-separated MODEL:REFERENCE pairs to flag as key paired contrasts",
                         "in rr_auc_contrasts.csv / rr_brier_contrasts.csv. Accepts manuscript",
                         "model ids (resolved via metadata$model_map) or raw model names.",
                         "[default: M4:M1,M5:M2,M7:M6 — the two co-primaries + confirmation]")),
  make_option("--skip_ci", action="store_true", default=FALSE,
              help="Skip confidence intervals (faster) [default: FALSE]"),
  make_option("--compute_integrated", action="store_true", default=TRUE,
              help="Compute iAUC/IBS & prediction-error curves [default: TRUE]"),
  make_option("--tau", type="numeric", default=NA,
              help="Integration horizon for iAUC/IBS (years). Default: max valid time."),
  make_option("--bootstrap", action="store_true", default=FALSE,
              help="Use bootstrap for CIs in integrated metrics (slower) [default: FALSE]"),
  make_option("--n_boot", type="integer", default=500,
              help="Bootstrap replicates if --bootstrap [default: 500]"),
  make_option("--save_plots", action="store_true", default=TRUE,
              help="Save prediction-error and AUC(t) curves [default: TRUE]"),
  make_option("--plot_formats", type="character", default="png",
              help="Comma-separated formats to save (pdf,png) [default: png]"),
  make_option("--plot_width", type="numeric", default=8,
              help="Plot width [inches] [default: 8]"),
  make_option("--plot_height", type="numeric", default=5,
              help="Plot height [inches] [default: 5]"),
  make_option("--outdir", type="character", default="output",
              help="Output directory [default: output]"),
  make_option("--prefix", type="character", default="",
              help="Prefix for output files [default: auto]"),
  make_option("--verbose", action="store_true", default=TRUE,
              help="Print detailed output [default: TRUE]")
)
opt <- parse_args(OptionParser(option_list=option_list))

stopifnot(!is.null(opt$fitted_models), !is.null(opt$data_processed))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

# Validate auc_engine option
auc_engine <- tolower(opt$auc_engine)
if (!auc_engine %in% c("riskregression", "timeroc", "both")) {
  stop("--auc_engine must be 'riskRegression', 'timeROC', or 'both'")
}

# Load timeROC only if needed
if (auc_engine %in% c("timeroc", "both")) {
  load_package("timeROC")
  cat(sprintf("  %-14s %s\n", "timeROC:", as.character(utils::packageVersion("timeROC"))))
}
cat("\n")

# ==============================================================================
# LOAD DATA
# ==============================================================================
cat("==============================================================================\n")
cat("TIME-DEPENDENT DISCRIMINATION + INTEGRATED METRICS\n")
cat("==============================================================================\n")
cat("AUC Engine:", opt$auc_engine, "\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

results <- readRDS(opt$fitted_models)
DT <- readRDS(opt$data_processed)

# Construct prefix (auto if not provided)
if (nchar(opt$prefix) > 0) {
  prefix <- paste0(opt$prefix, "_03_discrimination_metrics_")
} else {
  base_prefix_parts <- c()
  if (!is.null(results$metadata$options$cohort) && results$metadata$options$cohort != "") {
    base_prefix_parts <- c(base_prefix_parts, results$metadata$options$cohort)
  }
  if (!is.null(results$metadata$options$ancestry) && results$metadata$options$ancestry != "") {
    base_prefix_parts <- c(base_prefix_parts, results$metadata$options$ancestry)
  }
  if (!is.null(results$metadata$options$pheno_name) && results$metadata$options$pheno_name != "") {
    base_prefix_parts <- c(base_prefix_parts, results$metadata$options$pheno_name)
  }
  prefix <- if (length(base_prefix_parts) > 0) {
    paste0(paste(base_prefix_parts, collapse="_"), "_03_discrimination_metrics_")
  } else "03_discrimination_metrics_"
}
cat("Output prefix:", prefix, "\n")

# Log
log_file <- file.path(opt$outdir, paste0(prefix, "discrimination_metrics.log"))
log_con <- file(log_file, open="wt")
sink(log_con, type="output", split=opt$verbose)
sink(log_con, type="message")
cat("Log:", log_file, "\n\n")

models_list <- results$models
stopifnot(is.list(models_list), length(models_list) > 0)

cat("Models:", paste(names(models_list), collapse=", "), "\n")
cat("Sample size:", nrow(DT), "| Events:", sum(DT$PHENO, na.rm=TRUE), "\n\n")

# ==============================================================================
# HELPERS
# ==============================================================================
parse_times_cli <- function(x) {
  if (is.null(x) || nchar(x) == 0) return(numeric(0))
  tt <- suppressWarnings(as.numeric(trimws(strsplit(x, ",")[[1]])))
  tt[is.finite(tt)]
}

# ==============================================================================
# TIME SOURCE & G-HAT FILTERING
# ==============================================================================
# Choose requested candidate times
candidate_times <- NULL
if (isTRUE(opt$use_supported_times)) {
  cat("Using supported times saved by 01* (metadata$time_support)\n")
  ts <- results$metadata$time_support
  # Priority: eval_time_points -> requested_valid_time_points -> t_pass_auto
  if (!is.null(ts$eval_time_points) && length(ts$eval_time_points) > 0) {
    candidate_times <- as.numeric(ts$eval_time_points)
    cat("  Found eval_time_points in metadata:", paste(candidate_times, collapse=", "), "\n")
  } else if (!is.null(ts$requested_valid_time_points) && length(ts$requested_valid_time_points) > 0) {
    candidate_times <- as.numeric(ts$requested_valid_time_points)
    cat("  Using requested_valid_time_points from metadata:", paste(candidate_times, collapse=", "), "\n")
  } else if (!is.null(ts$t_pass_auto) && length(ts$t_pass_auto) > 0) {
    candidate_times <- as.numeric(ts$t_pass_auto)
    cat("  Using t_pass_auto from metadata:", paste(candidate_times, collapse=", "), "\n")
  } else {
    cat("  ⚠ No supported times found in metadata; falling back to --time_points\n")
    candidate_times <- parse_times_cli(opt$time_points)
  }
} else {
  candidate_times <- parse_times_cli(opt$time_points)
  cat("Using --time_points:", paste(candidate_times, collapse=", "), "\n")
}

# Basic sanity vs follow-up
max_time <- max(DT$TIME, na.rm = TRUE)
max_event_time <- max(DT$TIME[DT$PHENO == 1], na.rm = TRUE)

candidate_times <- sort(unique(candidate_times))
candidate_times <- candidate_times[candidate_times > 0]
beyond <- candidate_times[candidate_times > max_time]
if (length(beyond)) {
  cat("Removing time points beyond follow-up:", paste(beyond, collapse=", "), "\n")
  candidate_times <- setdiff(candidate_times, beyond)
}
stopifnot(length(candidate_times) > 0)

cat("Max follow-up:", round(max_time,2), "years | Max event time:", round(max_event_time,2), "years\n")
cat("Candidate time points:", paste(candidate_times, collapse=", "), "\n\n")

# Recompute G-hat and guards HERE (ensures current thresholds are respected)
cens_surv <- survfit(Surv(DT$TIME, 1 - DT$PHENO) ~ 1)
ghat_results <- data.table()
for (t in candidate_times) {
  ghat_t <- tryCatch({summary(cens_surv, times=t)$surv}, error=function(e) NA)
  n_at_risk <- sum(DT$TIME >= t)
  n_events_by_t <- sum(DT$PHENO == 1 & DT$TIME <= t)
  status <- if (is.na(ghat_t)) "No estimate"
  else if (ghat_t < opt$ghat_threshold) sprintf("FAIL (<%g)", opt$ghat_threshold)
  else if (n_at_risk < opt$min_at_risk) sprintf("FAIL (at_risk<%d)", opt$min_at_risk)
  else if (n_events_by_t < opt$min_events) sprintf("FAIL (events<%d)", opt$min_events)
  else "PASS"
  ghat_results <- rbind(ghat_results,
                        data.table(time_years = t,
                                   ghat = as.numeric(ghat_t),
                                   n_events_by_t = n_events_by_t,
                                   n_at_risk = n_at_risk,
                                   status = status))
}
fwrite(ghat_results, file.path(opt$outdir, paste0(prefix, "ghat_censoring.csv")))
valid_time_points <- ghat_results[status == "PASS", time_years]

if (length(valid_time_points) == 0) stop("No time points passed guards (G-hat, events, at-risk).")
cat("Valid times after guards:", paste(valid_time_points, collapse=", "), "\n\n")

# Integration horizon τ
tau <- if (is.na(opt$tau)) max(valid_time_points) else min(opt$tau, max(valid_time_points))
cat("Integration horizon (tau):", tau, "years\n\n")

# ==============================================================================
# OPTIONAL: TIME-DEPENDENT AUC (timeROC)
# ==============================================================================
if (auc_engine %in% c("timeroc", "both")) {
  cat("=== TIME-DEPENDENT AUC (timeROC) ===\n\n")
  
  get_se_auc <- function(roc, idx) {
    if (!is.null(roc$inference) && !is.null(roc$inference$vect_sd_1)) {
      se <- as.numeric(roc$inference$vect_sd_1[idx]); if (is.finite(se) && se > 0) return(se)
    }
    if (!is.null(roc$inference) && !is.null(roc$inference$vect_iid_AUC)) {
      col <- roc$inference$vect_iid_AUC[, idx]
      if (!all(is.na(col))) {
        se <- sqrt(sum(col^2, na.rm = TRUE)); if (is.finite(se) && se > 0) return(se)
      }
    }
    NA_real_
  }
  
  all_auc_results <- list()
  z_value <- qnorm(0.975)
  
  for (model_name in names(models_list)) {
    fit <- models_list[[model_name]]
    lp <- tryCatch(as.numeric(predict(fit, newdata=DT, type="lp")), error=function(e) NULL)
    if (is.null(lp) || length(unique(lp)) < 2) {
      cat("Skipping model", model_name, ": cannot get LP or no variation\n")
      next
    }
    model_results <- data.table()
    for (t in valid_time_points) {
      at_risk <- sum(DT$TIME >= t)
      events_by_t <- sum(DT$PHENO == 1 & DT$TIME <= t)
      if (at_risk < opt$min_at_risk || events_by_t < opt$min_events) next
      
      t_nudged <- min(t, max_event_time - 1e-8)
      out <- tryCatch({
        roc_obj <- timeROC(T=DT$TIME, delta=DT$PHENO, marker=lp, cause=1, times=t_nudged, iid=!opt$skip_ci)
        auc_idx <- which.min(abs(roc_obj$times - t_nudged))
        auc_t <- roc_obj$AUC[auc_idx]
        if (!opt$skip_ci) {
          se <- get_se_auc(roc_obj, auc_idx)
          lo <- if (is.finite(se)) auc_t - z_value*se else NA_real_
          hi <- if (is.finite(se)) auc_t + z_value*se else NA_real_
        } else { lo <- hi <- NA_real_ }
        list(auc=auc_t, lo=lo, hi=hi)
      }, error=function(e) list(auc=NA_real_, lo=NA_real_, hi=NA_real_))
      
      model_results <- rbind(model_results, data.table(
        model=model_name, time_years=t, n_at_risk=at_risk, n_events=events_by_t,
        AUC=round(out$auc,4),
        AUC_lower=if (!is.na(out$lo)) round(out$lo,4) else NA_real_,
        AUC_upper=if (!is.na(out$hi)) round(out$hi,4) else NA_real_
      ))
    }
    all_auc_results[[model_name]] <- model_results
  }
  
  tdauc_results <- rbindlist(all_auc_results, fill=TRUE)
  if (nrow(tdauc_results) > 0) {
    fwrite(tdauc_results, file.path(opt$outdir, paste0(prefix, "tdauc_results.csv")))
    cat("✓ Saved timeROC AUC(t): tdauc_results.csv\n\n")
  } else {
    cat("No timeROC AUC(t) computed.\n\n")
  }
}

# ==============================================================================
# INTEGRATED METRICS & PREDICTION-ERROR CURVES (riskRegression::Score)
# ==============================================================================
if (opt$compute_integrated || auc_engine %in% c("riskregression", "both")) {
  cat("=== INTEGRATED METRICS & AUC(t) (riskRegression) ===\n\n")
  
  # Use only coxph models so Score() can compute risks for Brier/IBS
  model_objs <- Filter(function(x) inherits(x, "coxph"), models_list)
  if (length(model_objs) == 0) stop("No coxph models available for Score().")

  # Restrict to the manuscript M0-M7 ladder by default. The fitted-model object
  # also carries the per-quantile helper models, and with contrasts=TRUE every
  # extra model adds a full row block to the pairwise contrast tables (16 models
  # -> 240 contrasts, most meaningless). --score_all_models restores the old set.
  if (!isTRUE(opt$score_all_models) && !is.null(results$metadata$model_map)) {
    # CORE ladder only. Clinical_* models are fitted on their own complete-case
    # sample (DT_clinical), so scoring them against the core data loaded here
    # would compare models estimated on different samples.
    mm_dt  <- as.data.table(results$metadata$model_map)
    ladder <- if ("ladder" %in% names(mm_dt)) as.character(mm_dt[ladder == "core", model])
              else as.character(mm_dt$model)
    ladder <- intersect(c("base", ladder), names(model_objs))
    if (length(ladder) >= 2) {
      dropped <- setdiff(names(model_objs), ladder)
      model_objs <- model_objs[ladder]
      if (length(dropped))
        cat("  Scoring the M0-M7 ladder only; excluded", length(dropped),
            "auxiliary model(s):", paste(dropped, collapse = ", "),
            "\n  (use --score_all_models to include them)\n")
    }
  }
  
  # Time grid for Score() = the GUARD-PASSING times only (G-hat, events, at-risk).
  # Do NOT augment with raw follow-up quantiles: those are unvetted and could feed
  # an under-supported horizon into AUC/Brier/iAUC/IBS for a sparse trajectory.
  # Fixed-horizon metrics (incl. the manuscript ΔAUC/ΔBrier contrasts) still run on
  # a single valid time; only the INTEGRATED metrics, which need >=2 points, are
  # skipped when the grid is too sparse.
  times_rr <- sort(unique(valid_time_points[valid_time_points <= tau & valid_time_points > 0]))
  if (length(times_rr) < 1)
    stop("No supported horizon passed the event/at-risk/G-hat guards; cannot compute discrimination.")
  do_integrated <- length(times_rr) >= 2
  if (!do_integrated)
    cat("  ⚠ Only", length(times_rr), "supported horizon; fixed-horizon AUC/Brier + contrasts",
        "reported, integrated iAUC/IBS skipped.\n")
  
  B_val <- if (opt$bootstrap) opt$n_boot else 0
  sc <- riskRegression::Score(
    object       = model_objs,
    formula      = Surv(TIME, PHENO) ~ 1,
    data         = DT,
    times        = times_rr,
    metrics      = c("AUC","Brier"),
    conf.int     = !opt$skip_ci,
    split.method = "none",
    null.model   = TRUE,
    B            = B_val,
    cens.model   = "km",
    # contrasts=TRUE makes Score() return PAIRED model differences with SE, CI and
    # p-value. Subtracting two separately-estimated AUCs is not a valid test; the
    # manuscript needs delta-AUC(M5-M2) and delta-AUC(M7-M6) with inference, and
    # Score() already computes them from the joint influence functions.
    contrasts    = TRUE
  )

  # ------------------- paired contrasts (delta AUC / delta Brier) -------------
  # Score(contrasts=TRUE) returns ALL pairwise model differences. We label only the
  # PRE-SPECIFIED NESTED contrasts (the manuscript reporting hierarchy) with a
  # `report_tier` (main_text / extended_data / supplementary) and a canonical
  # `contrast_label` ("Mlarger:Msmaller"); every other pairwise row is left with
  # report_tier = NA so it is not mistaken for a headline result. `is_key_contrast`
  # stays keyed to --contrast_pairs (default M4:M1,M5:M2,M7:M6) so the completeness banner
  # and run_status remain in sync with the LRT primaries (run_complete = M5vsM2 & M7vsM6).
  mm <- if (!is.null(results$metadata$model_map)) as.data.table(results$metadata$model_map) else NULL
  .n2i <- if (!is.null(mm) && all(c("model","model_id") %in% names(mm)))
    setNames(as.character(mm$model_id), tolower(trimws(as.character(mm$model)))) else character(0)
  nm2id <- function(nm) {                      # model NAME -> model_id, robust to case
    key <- tolower(trimws(nm))
    # `%in%` guard: [[ on a named atomic vector ERRORS for an absent name (unlike a
    # list). Brier contrasts include a "Null model" row not in model_map — return NA.
    if (length(.n2i) && key %in% names(.n2i)) return(unname(.n2i[[key]]))
    if (key == "base") return("M0")            # Score() names M0 "base"; model_map may say "Base"
    NA_character_
  }
  # Fixed reporting hierarchy, as canonical larger-minus-smaller model-id pairs.
  tier_map <- c("M4:M1"="main_text", "M5:M2"="main_text", "M7:M6"="main_text",
                "M7:M5"="extended_data", "M7:M4"="extended_data",
                "M1:M0"="supplementary", "M2:M0"="supplementary", "M3:M0"="supplementary")
  key_labels <- if (!is.null(opt$contrast_pairs) && nzchar(opt$contrast_pairs))
    trimws(strsplit(opt$contrast_pairs, ",")[[1]]) else character(0)

  # Label prespecified nested contrasts, normalizing any swapped orientation to
  # model−reference (negate delta, swap-negate the CI ((lower,upper) -> (-upper,-lower))).
  # Vectorised (no per-row := ), so it is robust across data.table versions.
  label_contrasts <- function(dt, delta_col) {
    dt[, `:=`(is_key_contrast = FALSE, contrast_label = NA_character_, report_tier = NA_character_)]
    if (!nrow(dt)) return(dt)
    mid <- vapply(as.character(dt$model),     nm2id, character(1))
    rid <- vapply(as.character(dt$reference), nm2id, character(1))
    pid <- paste0(mid, ":", rid); rpid <- paste0(rid, ":", mid)
    known <- unique(c(names(tier_map), key_labels))
    canon <- ifelse(pid %in% known, pid, ifelse(rpid %in% known, rpid, NA_character_))
    swap  <- !(pid %in% known) & (rpid %in% known)   # canonical pair is the reversed one
    if (any(swap)) {                                  # normalise to model−reference
      if (delta_col %in% names(dt)) dt[swap, (delta_col) := -get(delta_col)]
      if (all(c("lower","upper") %in% names(dt))) dt[swap, c("lower","upper") := .(-upper, -lower)]
      dt[swap, c("model","reference") := .(reference, model)]
    }
    lab <- !is.na(canon)
    if (any(lab))
      dt[lab, `:=`(contrast_label  = canon[lab],
                   report_tier     = unname(tier_map[canon[lab]]),   # NA if key-only (no tier)
                   is_key_contrast = canon[lab] %in% key_labels)]
    dt
  }
  # Track availability PER key pair PER metric — an aggregate "any found" flag
  # would read TRUE when M5:M2 is present but M7:M6 is missing (or vice versa).
  pair_labels <- key_labels
  found_by_pair <- list(AUC = setNames(rep(FALSE, length(pair_labels)), pair_labels),
                        Brier = setNames(rep(FALSE, length(pair_labels)), pair_labels))
  for (mtc in c("AUC","Brier")) {
    ct <- tryCatch(as.data.table(sc[[mtc]]$contrasts), error = function(e) NULL)
    if (is.null(ct) || !nrow(ct)) {
      cat(sprintf("  Note: %s contrasts not returned by this riskRegression version; skipped.\n", mtc))
      next
    }
    if ("times" %in% names(ct)) setnames(ct, "times", "time_years")
    ct <- label_contrasts(ct, paste0("delta.", mtc))
    for (lb in pair_labels) found_by_pair[[mtc]][[lb]] <- lb %in% ct$contrast_label
    miss <- pair_labels[!found_by_pair[[mtc]]]
    if (length(miss))
      cat(sprintf("  ⚠ Requested %s key contrast(s) not found in Score() output: %s\n",
                  mtc, paste(miss, collapse=", ")))
    n_nested <- sum(!is.na(ct$report_tier))
    fn <- file.path(opt$outdir, paste0(prefix, "rr_", tolower(mtc), "_contrasts.csv"))
    fwrite(ct, fn); cat("  ✓", basename(fn),
                        sprintf("(%d contrasts, %d key, %d prespecified nested)\n",
                                nrow(ct), sum(ct$is_key_contrast), n_nested))
  }
  # Discrimination completeness fragment for central QC — one column per pair×metric,
  # e.g. M5_vs_M2_delta_AUC_available. Column-safe pair id (":" -> "_vs_").
  .cid <- function(lb) gsub(":", "_vs_", lb, fixed = TRUE)
  disc_status <- data.table(key_contrasts_requested =
    if (length(pair_labels)) paste(pair_labels, collapse="; ") else "")
  for (lb in pair_labels) for (mtc in c("AUC","Brier"))
    disc_status[[paste0(.cid(lb), "_delta_", mtc, "_available")]] <- isTRUE(found_by_pair[[mtc]][[lb]])
  # Back-compat aggregate columns (TRUE only if ALL requested pairs found for the metric).
  disc_status[, paired_delta_AUC_available   := length(pair_labels) > 0 && all(unlist(found_by_pair$AUC))]
  disc_status[, paired_delta_Brier_available := length(pair_labels) > 0 && all(unlist(found_by_pair$Brier))]
  fwrite(disc_status, file.path(opt$outdir, paste0(prefix, "run_status_discrimination.csv")))

  # --------------------- AUC(t) by time (riskRegression) ----------------------
  auc_rr <- as.data.table(sc$AUC$score)
  if ("times" %in% names(auc_rr)) setnames(auc_rr, "times", "time_years")
  if (!"AUC" %in% names(auc_rr)) {
    auc_col <- grep("^AUC", names(auc_rr), value = TRUE)[1]
    if (!is.na(auc_col)) setnames(auc_rr, auc_col, "AUC")
  }
  # If Null model AUC not present (some versions), fill as 0.5
  if ("model" %in% names(auc_rr)) {
    null_rows <- auc_rr$model %in% c("Null","Baseline","Reference","CoxNoCovariate")
    if (any(null_rows) && any(is.na(auc_rr$AUC[null_rows]))) {
      auc_rr$AUC[null_rows] <- 0.5
      if ("lower" %in% names(auc_rr)) auc_rr$lower[null_rows] <- NA_real_
      if ("upper" %in% names(auc_rr)) auc_rr$upper[null_rows] <- NA_real_
    }
  }
  fwrite(auc_rr, file.path(opt$outdir, paste0(prefix, "rr_auc_by_time.csv")))
  cat("✓ Saved riskRegression AUC(t): rr_auc_by_time.csv\n")
  
  # --------------------- iAUC over [0, tau] (single number if available) -----
  iauc_rr <- NULL
  if (!is.null(sc$AUC$iauc)) {
    iauc_rr <- as.data.table(sc$AUC$iauc)
  } else if (!is.null(sc$AUC$integrated)) {
    iauc_rr <- as.data.table(sc$AUC$integrated)
  } else if (!is.null(sc$AUC$mean)) {
    tmp <- as.data.table(sc$AUC$mean)
    if (!"model" %in% names(tmp)) {
      rn <- rownames(sc$AUC$mean)
      tmp[, model := if (!is.null(rn)) rn else names(model_objs)]
    }
    numcols <- names(tmp)[sapply(tmp, is.numeric)]
    setnames(tmp, numcols[1], "iAUC")
    iauc_rr <- tmp
  }
  if (!is.null(iauc_rr)) {
    for (nm in c("se","lower","upper","SE","Lower","Upper")) {
      if (nm %in% names(iauc_rr)) setnames(iauc_rr, nm, tolower(nm))
    }
    iauc_rr[, tau := max(times_rr)]
    setcolorder(iauc_rr, unique(c("model","iAUC","lower","upper","se","tau", names(iauc_rr))))
    fwrite(iauc_rr, file.path(opt$outdir, paste0(prefix, "rr_iauc.csv")))
    cat("✓ Saved iAUC (single τ): rr_iauc.csv\n")
  } else {
    cat("⚠ iAUC not returned by Score() in this riskRegression version; skipping rr_iauc.csv\n")
  }
  
  # --------------------- Brier(t) prediction-error curve ----------------------
  brier_rr <- as.data.table(sc$Brier$score)
  if ("times" %in% names(brier_rr)) setnames(brier_rr, "times", "time_years")
  if (!"Brier" %in% names(brier_rr)) {
    bcol <- grep("^Brier", names(brier_rr), value = TRUE)[1]
    if (!is.na(bcol)) setnames(brier_rr, bcol, "Brier")
  }
  fwrite(brier_rr, file.path(opt$outdir, paste0(prefix, "rr_brier_curve.csv")))
  cat("✓ Saved Brier(t): rr_brier_curve.csv\n")
  
  # --------------------- IBS over [0, tau] (single number if available) ------
  ibs_rr <- NULL
  if (!is.null(sc$Brier$IBS)) {
    ibs_rr <- as.data.table(sc$Brier$IBS, keep.rownames = TRUE)
    setnames(ibs_rr, "rn", "model")
  } else if (!is.null(sc$Brier$integrated)) {
    ibs_rr <- as.data.table(sc$Brier$integrated)
  }
  if (!is.null(ibs_rr)) {
    if (!"IBS" %in% names(ibs_rr)) {
      icol <- intersect(c("IBS","integrated","mean","Brier.Int"), names(ibs_rr))
      if (length(icol) > 0) setnames(ibs_rr, icol[1], "IBS")
    }
    for (nm in c("se","lower","upper","SE","Lower","Upper")) {
      if (nm %in% names(ibs_rr)) setnames(ibs_rr, nm, tolower(nm))
    }
    if (!"model" %in% names(ibs_rr)) {
      rn <- rownames(as.data.frame(ibs_rr))
      ibs_rr[, model := if (!is.null(rn)) rn else names(model_objs)]
    }
    ibs_rr[, tau := max(times_rr)]
    setcolorder(ibs_rr, unique(c("model","IBS","lower","upper","se","tau", names(ibs_rr))))
    fwrite(ibs_rr, file.path(opt$outdir, paste0(prefix, "rr_ibs.csv")))
    cat("✓ Saved IBS (single τ): rr_ibs.csv\n")
  } else {
    cat("⚠ IBS not returned by Score() in this riskRegression version; skipping rr_ibs.csv\n")
  }
  
  # --------------------- NEW: per-horizon iAUC(0,t) and IBS(0,t) --------------
  cum_integrate_by_time <- function(curve_dt, value_col, start_value, group_col="model") {
    stopifnot(all(c("time_years", value_col, group_col) %in% names(curve_dt)))
    out_list <- list()
    for (m in unique(curve_dt[[group_col]])) {
      dd <- curve_dt[get(group_col) == m, .(time_years, val = get(value_col))]
      dd <- dd[is.finite(time_years) & is.finite(val)]
      setorder(dd, time_years)
      if (nrow(dd) == 0) next
      if (dd$time_years[1] > 0) {
        dd <- rbind(data.table(time_years = 0, val = start_value), dd)
      } else if (dd$time_years[1] == 0) {
        dd$val[1] <- ifelse(is.finite(dd$val[1]), dd$val[1], start_value)
      }
      area <- numeric(nrow(dd)); avg <- numeric(nrow(dd))
      area[1] <- 0; avg[1] <- start_value
      for (i in 2:nrow(dd)) {
        dt <- dd$time_years[i] - dd$time_years[i-1]
        trap <- 0.5 * dt * (dd$val[i] + dd$val[i-1])
        area[i] <- area[i-1] + trap
        avg[i]  <- if (dd$time_years[i] > 0) area[i] / dd$time_years[i] else start_value
      }
      tmp <- data.table(model = m,
                        time_years = dd$time_years,
                        integrated = avg)
      out_list[[m]] <- tmp[time_years > 0]
    }
    rbindlist(out_list, fill=TRUE)
  }
  
  # iAUC(0,t) / IBS(0,t) — integrated metrics need >=2 supported horizons. Skipped
  # (with the warning above) on a sparse trajectory rather than integrating over an
  # unvetted grid.
  if (do_integrated) {
  iauc_by_time <- cum_integrate_by_time(auc_rr, value_col="AUC", start_value=0.5)
  setnames(iauc_by_time, "integrated", "iAUC_at_tau")
  if (all(c("lower","upper") %in% names(auc_rr))) {
    # These are the integrals of the POINTWISE lower/upper AUC(t) curves. They are
    # a descriptive envelope, NOT a 95% CI for iAUC: integrating pointwise bounds
    # ignores the covariance of AUC(t) across time. For iAUC/IBS inference use the
    # Score()-returned single-tau estimates (rr_iauc.csv / rr_ibs.csv, influence-
    # function CIs) or bootstrap the whole integral. The fixed-horizon paired
    # ΔAUC/ΔBrier contrasts remain the primary inferential result.
    lower_curve <- copy(auc_rr); setnames(lower_curve, "lower", "AUC")
    iauc_lo <- cum_integrate_by_time(lower_curve[, .(model, time_years, AUC)], "AUC", 0.5)
    upper_curve <- copy(auc_rr); setnames(upper_curve, "upper", "AUC")
    iauc_hi <- cum_integrate_by_time(upper_curve[, .(model, time_years, AUC)], "AUC", 0.5)
    setnames(iauc_lo, "integrated", "iAUC_pointwise_env_lower")
    setnames(iauc_hi, "integrated", "iAUC_pointwise_env_upper")
    iauc_by_time <- merge(iauc_by_time, iauc_lo, by=c("model","time_years"), all.x=TRUE)
    iauc_by_time <- merge(iauc_by_time, iauc_hi, by=c("model","time_years"), all.x=TRUE)
  }
  fwrite(iauc_by_time, file.path(opt$outdir, paste0(prefix, "rr_iauc_by_time.csv")))
  cat("✓ Saved per-horizon iAUC(0,t): rr_iauc_by_time.csv\n")
  
  # IBS(0,t)
  ibs_by_time <- cum_integrate_by_time(brier_rr, value_col="Brier", start_value=0.0)
  setnames(ibs_by_time, "integrated", "IBS_at_tau")
  if (all(c("lower","upper") %in% names(brier_rr))) {
    lower_b <- copy(brier_rr); setnames(lower_b, "lower", "Brier")
    ibs_lo <- cum_integrate_by_time(lower_b[, .(model, time_years, Brier)], "Brier", 0.0)
    upper_b <- copy(brier_rr); setnames(upper_b, "upper", "Brier")
    ibs_hi <- cum_integrate_by_time(upper_b[, .(model, time_years, Brier)], "Brier", 0.0)
    # Pointwise-bound envelope, not a CI — see the iAUC note above.
    setnames(ibs_lo, "integrated", "IBS_pointwise_env_lower")
    setnames(ibs_hi, "integrated", "IBS_pointwise_env_upper")
    ibs_by_time <- merge(ibs_by_time, ibs_lo, by=c("model","time_years"), all.x=TRUE)
    ibs_by_time <- merge(ibs_by_time, ibs_hi, by=c("model","time_years"), all.x=TRUE)
  }
  fwrite(ibs_by_time, file.path(opt$outdir, paste0(prefix, "rr_ibs_by_time.csv")))
  cat("✓ Saved per-horizon IBS(0,t): rr_ibs_by_time.csv\n")
  }  # end do_integrated

  # --------------------- Optional plots ---------------------------------------
  if (opt$save_plots) {
    plot_formats_parsed <- trimws(strsplit(opt$plot_formats, ",")[[1]])

    p_pe <- ggplot(brier_rr, aes(x = time_years, y = Brier, color = model)) +
      geom_line(linewidth = 1) +
      { if (!opt$skip_ci && all(c("lower","upper") %in% names(brier_rr)))
        geom_ribbon(aes(ymin = lower, ymax = upper, fill = model), alpha = 0.15, color = NA) } +
      labs(title = sprintf("Prediction-Error Curves (IBS over [0, %.1f])", max(times_rr)),
           x = "Time (years)", y = "Brier score", color = "Model", fill = "Model") +
      theme_bw(base_size = 12) + theme(legend.position = "bottom")
    save_plot_cairo(p_pe, file.path(opt$outdir, paste0(prefix, "pe_curve.png")), width = opt$plot_width, height = opt$plot_height, dpi = 300, formats = plot_formats_parsed)
    cat("✓ Saved plot: pe_curve.png\n")

    p_auc <- ggplot(auc_rr, aes(x = time_years, y = AUC, color = model)) +
      geom_line(linewidth = 1) +
      { if (!opt$skip_ci && all(c("lower","upper") %in% names(auc_rr)))
        geom_ribbon(aes(ymin = lower, ymax = upper, fill = model), alpha = 0.15, color = NA) } +
      labs(title = sprintf("Time-dependent AUC (iAUC over [0, %.1f])", max(times_rr)),
           x = "Time (years)", y = "AUC(t)", color = "Model", fill = "Model") +
      theme_bw(base_size = 12) + theme(legend.position = "bottom")
    save_plot_cairo(p_auc, file.path(opt$outdir, paste0(prefix, "auc_time_curve.png")), width = opt$plot_width, height = opt$plot_height, dpi = 300, formats = plot_formats_parsed)
    cat("✓ Saved plot: auc_time_curve.png\n")
  }
  cat("\n")
}

# ==============================================================================
# MERGED SUMMARY
# ==============================================================================
cat("=== CREATING DISCRIMINATION SUMMARY ===\n\n")
discrimination_summary <- NULL

if (exists("auc_rr") && nrow(auc_rr) > 0) {
  discrimination_summary <- auc_rr[, .(model, time_years, AUC)]
  if (all(c("lower","upper") %in% names(auc_rr))) {
    discrimination_summary[, `:=`(AUC_lower = auc_rr$lower, AUC_upper = auc_rr$upper)]
  }
  cat("✓ Added riskRegression AUC(t)\n")
}

if (exists("tdauc_results") && nrow(tdauc_results) > 0 && auc_engine == "timeroc") {
  discrimination_summary <- tdauc_results[, .(model, time_years, AUC, AUC_lower, AUC_upper)]
  cat("✓ Added timeROC AUC(t)\n")
}

if (exists("brier_rr") && nrow(brier_rr) > 0) {
  br_keep <- brier_rr[, .(model, time_years, Brier)]
  if (all(c("lower","upper") %in% names(brier_rr))) {
    br_keep[, `:=`(Brier_lower = brier_rr$lower, Brier_upper = brier_rr$upper)]
  }
  discrimination_summary <- if (is.null(discrimination_summary)) br_keep
  else merge(discrimination_summary, br_keep, by = c("model","time_years"), all = TRUE)
  cat("✓ Added Brier(t)\n")
}

if (exists("iauc_by_time") && nrow(iauc_by_time) > 0) {
  discrimination_summary <- if (is.null(discrimination_summary)) iauc_by_time
  else merge(discrimination_summary, iauc_by_time, by=c("model","time_years"), all=TRUE)
  cat("✓ Added iAUC(0,t)\n")
}
if (exists("ibs_by_time") && nrow(ibs_by_time) > 0) {
  discrimination_summary <- if (is.null(discrimination_summary)) ibs_by_time
  else merge(discrimination_summary, ibs_by_time, by=c("model","time_years"), all=TRUE)
  cat("✓ Added IBS(0,t)\n")
}

if (!is.null(discrimination_summary) && nrow(discrimination_summary) > 0) {
  ghat_small <- ghat_results[status == "PASS", .(time_years, ghat)]
  discrimination_summary <- merge(discrimination_summary, ghat_small, by = "time_years", all.x = TRUE)
  setorder(discrimination_summary, model, time_years)
  fwrite(discrimination_summary, file.path(opt$outdir, paste0(prefix, "discrimination_summary.csv")))
  cat("✓ Saved combined discrimination_summary.csv\n\n")
} else {
  cat("⚠ Nothing to save in discrimination_summary.csv (no curves computed)\n\n")
}

# ==============================================================================
# SUMMARY & TEARDOWN
# ==============================================================================
cat("==============================================================================\n")
cat("DISCRIMINATION ANALYSIS COMPLETE\n")
cat("==============================================================================\n\n")
cat("AUC Engine:", opt$auc_engine, "\n")
cat("Sample size:", nrow(DT), " | Events:", sum(DT$PHENO, na.rm=TRUE), "\n")
cat("Valid time points (after guards):", paste(valid_time_points, collapse=", "), "\n")
if (opt$compute_integrated) cat("Integration tau:", tau, "years\n")
cat("Outputs saved to:", opt$outdir, "\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("==============================================================================\n")

sink(type="message"); sink(type="output"); close(log_con)
