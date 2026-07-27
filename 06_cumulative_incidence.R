#!/usr/bin/env Rscript
# ==============================================================================
# 06_cumulative_incidence.R - Cumulative incidence analysis
# ==============================================================================
# Purpose:
#   - Calculate cumulative incidence curves by PRS strata
#   - Compute absolute risk differences between strata
#   - Statistical testing (log-rank across PRS strata)
#   - Risk advancement period (RAP)
#   - Population attributable fraction (PAF)         [opt-in]
#   - Bootstrap confidence intervals
#   - Number needed to screen (NNS)                 [opt-in]
#   - Decision curves (net benefit)                  [opt-in]
#   - Max observed cumulative risk                    [opt-in]
#   - Use pre-vetted time points passing G-hat/events/at-risk guards saved by 01*R
#
# ------------------------------------------------------------------------------
# IMPORTANT — THIS IS NOT A COMPETING-RISKS ANALYSIS
# ------------------------------------------------------------------------------
# Despite the file name, this script treats the outcome as a SINGLE event type:
#   - sanitize_surv_data() collapses status to binary (PHENO != 0), so any
#     multi-cause coding is lost;
#   - "cumulative incidence" here is 1 - KM, not a cause-specific CIF, and it
#     therefore ignores competing events (including death);
#   - the stratum comparison is therefore a LOG-RANK test on those 1-KM curves,
#     not Gray's test (which cmprsk::cuminc previously mislabelled it as).
#
# With a single event type and independent censoring this is the correct
# quantity. If competing risks matter for your outcome — in particular if death
# is common before the event of interest — 1 - KM OVERSTATES absolute risk and a
# genuine cause-specific / Fine-Gray module is required instead.
# ------------------------------------------------------------------------------
#
# Inputs:
#   - fitted_models.rds (from 01_fit_models.R)
#   - data_processed.rds (from 01_fit_models.R)
#   - requested_timepoints_ghat.csv (optional; from 01*R)
#
# Author: Ying Wang; updated to consume saved G-hat–passing time points
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

required_packages <- c("optparse", "data.table", "survival", "ggplot2",
                       "boot", "scales", "rlang")  # Cairo omitted (optional; see prs_risk_utils.R)
invisible(lapply(required_packages, load_package))

# Defined before first use. Previously this appeared two lines AFTER its first
# call site and worked only because rlang happens to export the same operator.
`%||%` <- function(a, b) if (is.null(a)) b else a

# Shared risk-prediction helper, located next to this script.
# Prefer PRS_SCRIPT_DIR (exported by the driver) so we never parse --file, whose spaces
# some macOS/OneDrive setups encode as ~+~ (which breaks normalizePath()). Fall back to
# --file (with ~+~ decoded), then to getwd().
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

# ==============================================================================
# CLI
# ==============================================================================
option_list <- list(
  make_option("--models_file", type="character",
              help="Path to fitted_models.rds from 01_fit_models.R [REQUIRED]"),
  make_option("--data_file", type="character", default=NULL,
              help="Path to data_processed.rds [default: auto-detect beside models_file]"),
  make_option("--score_role", type="character", default="all",
              help=paste("Which PRS drives the risk strata: 'all' (default), or a comma-separated",
                         "subset of onset,outcome,progression. Legacy behaviour was the primary",
                         "(ONSET) score only.")),
  make_option("--role_quantiles", type="character", default="Q5,Q10",
              help=paste("Quantile schemes emitted per role [default: Q5,Q10]; 'all' for every",
                         "scheme. Extreme-tail schemes stay on the primary score only.")),
  make_option("--enable_clinical_utility", action="store_true", default=FALSE,
              help=paste("Emit the clinical-utility outputs: net_benefit, clinical_utility_metrics,",
                         "population_attributable_fraction, nns_*, lifetime_risks. These are",
                         "APPARENT (fit and evaluated in the same sample) and, for PAF/NNS, assume",
                         "an intervention that is not defined here. Descriptive risk outputs",
                         "(absolute_risks, ci_curves, risk_table, risk_differences, median_times,",
                         "risk_advancement_period) are always written. [default: FALSE]")),
  make_option("--risk_model", type="character", default=NULL,
              help=paste("Model used for decision curves / predicted risk. Accepts a manuscript id",
                         "(M0-M7, resolved via metadata$model_map) or a raw fitted-model name.",
                         "Default is the manuscript primary comparison model (M5) when available,",
                         "otherwise 'minimal'. NOTE: 'minimal' carries only the PRIMARY PRS, which",
                         "is the onset score - rarely what you want here.")),
  make_option("--outdir", type="character", default=NULL,
              help="Output directory [default: same as models_file]"),
  make_option("--prefix", type="character", default="",
              help="Prefix for output files [default: inferred from models_file/metadata]"),
  
  # Time points behavior. default=FALSE so the driver's --use_supported_times flag
  # (emitted only when --use-supported-times=1) fully controls it: setting 0 in the
  # driver now actually disables supported-time selection here. Stages 03/05 match.
  make_option("--use_supported_times", action="store_true", default=FALSE,
              help="Use time points saved by 01*R (requested_timepoints_ghat.csv or metadata aliases). The driver passes this when --use-supported-times=1 (its default)."),
  make_option("--time_points", type="character", default="5,10,15,20",
              help="Fallback comma-separated time points (years) if no saved times are found"),
  make_option("--max_time", type="numeric", default=NA_real_,
              help="Max follow-up to display; default: 95th percentile of TIME"),
  
  # Analysis settings
  make_option("--reference_quantile", type="character", default="Q1",
              help="Reference stratum label for risk differences (e.g., Q1)"),
  make_option("--high_risk_threshold", type="numeric", default=0.20,
              help="Threshold for defining high risk (decision curves)"),
  
  # Bootstrap
  make_option("--bootstrap_risk_differences", action="store_true", default=FALSE,
              help=paste("Bootstrap CIs for between-stratum risk differences. OFF by default: it",
                         "resamples survfit per role x quantile x group x horizon and dominates",
                         "runtime. Point estimates are always produced. [default: FALSE]")),
  make_option("--min_cell_count", type="integer", default=10,
              help=paste("Suppress risk/curve rows whose at-risk count is below this, for",
                         "disclosure control in small strata; a minimum_cell_count_pass flag",
                         "records the outcome. [default: 10]")),
  make_option("--n_bootstrap", type="integer", default=1000,
              help="Bootstrap iterations for CIs"),
  make_option("--bootstrap_seed", type="integer", default=123,
              help="Seed for bootstrap"),
  
  # RAP thresholds
  make_option("--rap_thresholds", type="character", default="0.30,0.50,0.70",
              help="Comma-separated risk levels for RAP (proportions or percents)"),
  
  # Decision curves
  make_option("--decision_curves_all_times", action="store_true", default=FALSE,
              help="Generate decision curves for all time horizons (default: only first)"),
  
  # Plots
  make_option("--generate_plots", action="store_true", default=FALSE,
              help="Generate plots (PNG)"),
  make_option("--plot_formats", type="character", default="png",
              help="Comma-separated formats to save (pdf,png) [default: png]"),
  make_option("--figure_width", type="numeric", default=10,
              help="Figure width (in)"),
  make_option("--figure_height", type="numeric", default=8,
              help="Figure height (in)"),
  make_option("--dpi", type="numeric", default=300,
              help="Figure resolution (DPI)"),
  
  make_option("--verbose", action="store_true", default=TRUE,
              help="Print detailed output")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$models_file)) stop("--models_file is required")

# ==============================================================================
# SETUP & LOAD
# ==============================================================================
cat("==============================================================================\n")
cat("ENHANCED CUMULATIVE INCIDENCE ANALYSIS\n")
cat("==============================================================================\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("==============================================================================\n\n")

if (!file.exists(opt$models_file)) stop("Models file not found: ", opt$models_file)
results <- readRDS(opt$models_file)
cat("✓ Loaded fitted models\n")

fitted_models <- results$models
metadata      <- results$metadata
# See 02_kaplan_meier.R for the rationale: metadata$quantiles is the PRIMARY
# (onset) score only; quantiles_by_role carries all three, keyed so the existing
# "PRS_<key>" lookup resolves to the PRS_<role>_<scheme> columns.
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

# data file
if (is.null(opt$data_file)) {
  models_dir <- dirname(opt$models_file)
  base <- basename(opt$models_file)
  data_file <- file.path(models_dir, sub("fitted_models", "data_processed", base))
} else {
  data_file <- opt$data_file
}
if (!file.exists(data_file)) stop("Data file not found: ", data_file)
DT <- readRDS(data_file)
cat("✓ Loaded processed data (N=", nrow(DT), ")\n\n", sep="")

# columns
timeCol   <- metadata$columns$time %||% "TIME"
statusCol <- metadata$columns$status %||% "PHENO"

# outdir/prefix
if (is.null(opt$outdir)) opt$outdir <- dirname(opt$models_file)
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

if (nchar(opt$prefix) > 0) {
  prefix <- paste0(opt$prefix, "_")
} else {
  parts <- c()
  if (!is.null(metadata$options$cohort)     && metadata$options$cohort    != "") parts <- c(parts, metadata$options$cohort)
  if (!is.null(metadata$options$ancestry)   && metadata$options$ancestry  != "") parts <- c(parts, metadata$options$ancestry)
  if (!is.null(metadata$options$pheno_name) && metadata$options$pheno_name!= "") parts <- c(parts, metadata$options$pheno_name)
  prefix <- if (length(parts)>0) paste0(paste(parts, collapse="_"), "_06_cumulative_incidence_") else "06_cumulative_incidence_"
}
cat("Output prefix:", prefix, "\n")
cat("Generate plots:", opt$generate_plots, "\n\n")

# Logging
log_file <- file.path(opt$outdir, paste0(prefix, "log.txt"))
log_con  <- file(log_file, open="wt")
sink(log_con, type="output", split=opt$verbose)
sink(log_con, type="message")
cat("Log file:", log_file, "\n\n")

# ==============================================================================
# TIME POINTS: LOAD SAVED G-HAT–PASSING TIMES FROM 01*R (preferred)
# ==============================================================================
# Sources we try, in order:
# 1) requested_timepoints_ghat.csv (if present)
# 2) metadata$uno_supported_times / metadata$time_points_supported / metadata$time_support$t_pass
# 3) CLI --time_points fallback

read_saved_requested_times <- function(models_file) {
  dd <- dirname(models_file)
  # Try to detect the base prefix used in 01*R (replace "fitted_models.rds")
  base <- basename(models_file)
  file_stub <- sub("fitted_models\\.rds$", "", base)
  # Common filename
  cand <- c(
    file.path(dd, paste0(file_stub, "requested_timepoints_ghat.csv")),
    file.path(dd, "requested_timepoints_ghat.csv")
  )
  cand <- cand[file.exists(cand)]
  if (!length(cand)) return(NULL)
  out <- tryCatch(fread(cand[1]), error=function(e) NULL)
  if (is.null(out)) return(NULL)
  # Expect columns: time_years, pass (logical/0-1). If no pass column, assume all pass.
  if (!"time_years" %in% names(out)) {
    # try "time" or first col
    if ("time" %in% names(out)) setnames(out, "time", "time_years")
    else setnames(out, names(out)[1], "time_years")
  }
  if (!"pass" %in% names(out)) out[, pass := TRUE]
  out[pass == TRUE & is.finite(time_years), sort(unique(time_years))]
}

find_supported_times_from_metadata <- function(res, DT) {
  cand <- list()
  if (!is.null(res$metadata$uno_supported_times)) cand[[length(cand)+1]] <- res$metadata$uno_supported_times
  if (!is.null(res$metadata$time_points_supported)) cand[[length(cand)+1]] <- res$metadata$time_points_supported
  if (!is.null(res$metadata$time_support) && !is.null(res$metadata$time_support$t_pass))
    cand[[length(cand)+1]] <- res$metadata$time_support$t_pass
  if (!is.null(res$stats) && length(res$stats)) {
    per_model <- unique(unlist(lapply(res$stats, function(x) {
      if (!is.null(x$uno_supported_times)) return(x$uno_supported_times)
      if (!is.null(x$time_points_supported)) return(x$time_points_supported)
      NULL
    })))
    if (length(per_model)) cand[[length(cand)+1]] <- per_model
  }
  if (!length(cand)) return(NULL)
  merged <- unique(sort(na.omit(as.numeric(unlist(cand)))))
  if (!length(merged)) return(NULL)
  merged <- merged[merged > 0 & merged <= max(DT[[timeCol]], na.rm=TRUE)]
  if (!length(merged)) return(NULL)
  merged
}

if (isTRUE(opt$use_supported_times)) {
  saved_times <- read_saved_requested_times(opt$models_file)
  if (is.null(saved_times)) saved_times <- find_supported_times_from_metadata(results, DT)
} else saved_times <- NULL

if (!is.null(saved_times) && length(saved_times) > 0) {
  time_points <- saved_times
  cat("Using saved time points passing guards:", paste(time_points, collapse=", "), "years\n\n")
} else {
  # Fallback to CLI
  tp <- as.numeric(trimws(strsplit(opt$time_points, ",")[[1]]))
  time_points <- tp[is.finite(tp) & tp > 0]
  if (!length(time_points)) stop("No valid time points available (neither saved nor CLI).")
  cat("Using user-specified time points (--time_points):", paste(time_points, collapse=", "), "years\n\n")
}

# Max time
if (!is.finite(opt$max_time) || is.na(opt$max_time)) {
  opt$max_time <- as.numeric(stats::quantile(DT[[timeCol]], 0.95, na.rm=TRUE))
}
cat("Max follow-up for display:", round(opt$max_time,2), "years\n\n")

# RAP thresholds parse (accept 0.3 or 30)
parse_rap_thresholds <- function(s) {
  raw <- as.numeric(trimws(strsplit(s, ",")[[1]]))
  raw <- raw[is.finite(raw) & raw > 0]
  if (!length(raw)) return(c(0.3, 0.5, 0.7))
  out <- ifelse(raw > 1, raw/100, raw)
  pmin(pmax(out, 1e-6), 0.999)
}
rap_thresholds <- parse_rap_thresholds(opt$rap_thresholds)

cat("Analysis settings:\n")
cat("  Time points:", paste(time_points, collapse=", "), "years\n")
cat("  Reference label:", opt$reference_quantile, "\n")
cat("  Bootstrap iterations:", opt$n_bootstrap, "\n")
cat("  RAP thresholds:", paste(sprintf("%.0f%%", rap_thresholds*100), collapse=", "), "\n")
cat("  Decision curves all times:", opt$decision_curves_all_times, "\n\n")

set.seed(opt$bootstrap_seed)

# ==============================================================================
# HELPERS
# ==============================================================================
resolve_strata_column <- function(DT, quantile_list, key) {
  cand <- NULL
  if (!is.null(quantile_list) && length(quantile_list) > 0 && key %in% names(quantile_list)) {
    val <- quantile_list[[key]]
    if (is.character(val) && length(val)==1) cand <- val
  }
  hits <- character(0)
  if (!is.null(cand) && cand %in% names(DT)) hits <- c(hits, cand)
  patterns <- unique(c(
    paste0("^", key, "$"),
    paste0("^", toupper(key), "$"),
    paste0("^PRS_", key, "$"),
    paste0("^PRS_", toupper(key), "$"),
    paste0("^", key, "_group$"),
    paste0("^PRS_", key, "_group$"),
    paste0("^PRS_", key, "_QUANT$"),
    paste0("^", sub("^quantile_", "PRS_", key), "$"),
    paste0("^", sub("^quantile_", "PRS_", key), "_group$")
  ))
  for (p in patterns) hits <- unique(c(hits, grep(p, names(DT), value=TRUE)))
  if (!length(hits)) {
    likely <- names(DT)[sapply(DT, function(x) is.factor(x) || is.character(x))]
    likely <- grep("PRS|Q[0-9]+|quantile|extreme|group", likely, value=TRUE, ignore.case=TRUE)
    likely <- likely[sapply(likely, function(cn) length(unique(DT[[cn]]))>=2)]
    hits <- unique(c(hits, likely))
  }
  if (!length(hits)) return(NA_character_)
  hits_key <- hits[grepl(key, hits, ignore.case=TRUE)]
  if (length(hits_key)>0) return(hits_key[1])
  hits[1]
}

sanitize_surv_data <- function(DT, strata_col, time_col = "TIME", status_col = "PHENO") {
  stopifnot(strata_col %in% names(DT), time_col %in% names(DT), status_col %in% names(DT))
  D <- data.table::copy(DT[, .SD, .SDcols = c(strata_col, time_col, status_col)])
  setnames(D, c(strata_col, time_col, status_col), c("STRATA","TIME","PHENO"))
  if (!is.numeric(D$PHENO)) D[, PHENO := as.integer(PHENO)]
  D[is.na(PHENO), PHENO := 0L]
  D[, PHENO := as.integer(PHENO != 0)]
  D <- D[!is.na(STRATA) & !is.na(TIME)]
  if (!is.factor(D$STRATA)) D[, STRATA := factor(STRATA)]
  D[]
}

get_strata_survival <- function(DT, strata_col, time_col="TIME", status_col="PHENO") {
  if (!is.factor(DT[[strata_col]])) DT[[strata_col]] <- factor(DT[[strata_col]])
  lv <- levels(DT[[strata_col]])
  out <- vector("list", length(lv))
  for (i in seq_along(lv)) {
    st <- lv[i]
    df <- DT[get(strata_col)==st]
    if (nrow(df) < 2) next
    fit <- survfit(Surv(get(time_col), get(status_col)) ~ 1, data=df)
    out[[i]] <- data.table(
      stratum = st,
      time = fit$time,
      surv = fit$surv,
      surv_lower = fit$lower,
      surv_upper = fit$upper,
      n_risk = fit$n.risk,
      n_event = fit$n.event
    )
  }
  rbindlist(out, fill=TRUE)
}

# FIXED-HORIZON LOOKUP ON A STEP FUNCTION.
# A Kaplan-Meier curve is a step function: S is constant on [t_i, t_{i+1}). The
# value at horizon tp is therefore the estimate at the LAST event time <= tp.
# The previous nearest-neighbour match, which.min(abs(sd$time - tp)), could
# return an estimate from AFTER the horizon -- ask for 5 years, get the 5.1-year
# risk -- silently inflating every published absolute risk. findInterval() gives
# the correct left-continuous lookup, and a horizon beyond the last observed
# event time returns NA rather than the carried-forward final value, matching the
# extrapolate=FALSE semantics in prs_risk_utils.R.
.km_at_horizon <- function(sd, tp) {
  if (nrow(sd) == 0) return(NULL)
  tmax <- max(sd$time, na.rm = TRUE)
  if (is.finite(tp) && tp > tmax)
    return(list(ci = NA_real_, cil = NA_real_, ciu = NA_real_,
                n_risk = NA_integer_, supported = FALSE))
  idx <- findInterval(tp, sd$time)
  if (idx == 0L)      # before the first event time: S = 1, risk = 0
    return(list(ci = 0, cil = 0, ciu = 0,
                n_risk = sd$n_risk[1], supported = TRUE))
  list(ci      = 1 - sd$surv[idx],
       cil     = 1 - sd$surv_upper[idx],
       ciu     = 1 - sd$surv_lower[idx],
       n_risk  = sd$n_risk[idx],
       supported = TRUE)
}

calculate_absolute_risks <- function(surv_data, time_points, min_cell_count = 10L) {
  res <- list(); n_unsupported <- 0L; n_suppressed <- 0L
  for (tp in time_points) {
    for (st in unique(surv_data$stratum)) {
      sd <- surv_data[stratum == st]
      v  <- .km_at_horizon(sd, tp)
      if (is.null(v)) next
      if (!isTRUE(v$supported)) n_unsupported <- n_unsupported + 1L
      # Disclosure control: at-risk below the threshold -> suppress the estimate
      # and flag it, so exact small-group risks are never published.
      cell_ok <- is.finite(v$n_risk) && v$n_risk >= min_cell_count
      if (!cell_ok) n_suppressed <- n_suppressed + 1L
      res[[length(res)+1]] <- data.table(
        stratum = st, time = tp,
        ci = if (cell_ok) v$ci else NA_real_,
        ci_lower = if (cell_ok) v$cil else NA_real_,
        ci_upper = if (cell_ok) v$ciu else NA_real_,
        # Suppress the count itself on failing rows: printing n_risk = 7 would
        # disclose exactly what the threshold is meant to hide.
        n_risk = if (cell_ok) v$n_risk else NA_integer_,
        supported = v$supported,
        minimum_cell_count_pass = cell_ok
      )
    }
  }
  if (n_unsupported > 0)
    cat(sprintf("  ⚠ %d stratum/horizon combination(s) exceed observed follow-up; risk is NA (not 0).\n",
                n_unsupported))
  if (n_suppressed > 0)
    cat(sprintf("  ⚠ %d stratum/horizon combination(s) suppressed (< %d at risk) for disclosure control.\n",
                n_suppressed, min_cell_count))
  rbindlist(res, fill=TRUE)
}

calculate_risk_differences <- function(absolute_risks, ref_label) {
  if (nrow(absolute_risks)==0) return(data.table())
  # Horizons beyond observed follow-up now carry ci = NA (see .km_at_horizon).
  # Drop them before differencing so an unsupported horizon is reported as absent
  # rather than silently becoming an NA "difference" that looks like a result.
  n_drop <- sum(!is.finite(absolute_risks$ci))
  if (n_drop > 0) {
    cat(sprintf("    %d unsupported stratum/horizon row(s) excluded from risk differences.\n", n_drop))
    absolute_risks <- absolute_risks[is.finite(ci)]
    if (nrow(absolute_risks)==0) return(data.table())
  }
  ref <- absolute_risks[stratum==ref_label]
  if (nrow(ref)==0) return(data.table())
  setnames(ref, c("ci","ci_lower","ci_upper"), c("ref_ci","ref_ci_lower","ref_ci_upper"))
  res <- merge(absolute_risks, ref[, .(time, ref_ci, ref_ci_lower, ref_ci_upper)], by="time", allow.cartesian=TRUE)
  res <- res[stratum != ref_label]
  # These bounds subtract two MARGINAL KM confidence limits — they are conservative
  # marginal bounds, NOT a paired 95% CI for the risk difference (which ignores the
  # covariance of the two group estimates). Named accordingly; a valid paired CI
  # comes only from --bootstrap_risk_differences.
  res[, `:=`(
    risk_diff = ci - ref_ci,
    risk_diff_marginal_lower = ci_lower - ref_ci_upper,
    risk_diff_marginal_upper = ci_upper - ref_ci_lower
  )]
  res
}

calculate_nns <- function(risk_differences) {
  if (nrow(risk_differences)==0) return(data.table())
  nd <- copy(risk_differences)
  nd <- nd[is.finite(risk_diff) & risk_diff > 0]
  if (nrow(nd)==0) return(data.table())
  nd[, nns := ceiling(1 / risk_diff)]
  nd[, nns_lower := ifelse(is.finite(risk_diff_marginal_upper) & risk_diff_marginal_upper>0, ceiling(1/risk_diff_marginal_upper), NA)]
  nd[, nns_upper := ifelse(is.finite(risk_diff_marginal_lower) & risk_diff_marginal_lower>0, ceiling(1/risk_diff_marginal_lower), NA)]
  nd[, .(stratum, time, risk_diff, nns, nns_lower, nns_upper)]
}

# STRATUM COMPARISON — LOG-RANK, NOT GRAY'S TEST.
# Gray's test compares cumulative INCIDENCE functions across competing causes.
# This pipeline models a single event type (sanitize_surv_data() collapses status
# to binary and death is not a competing event), so feeding it to cmprsk::cuminc
# produced a test that was neither labelled nor interpretable as advertised.
# survdiff() is the correct comparison of these 1-KM curves across PRS strata.
perform_logrank_test <- function(D, strata_col="STRATA", time_col="TIME", status_col="PHENO") {
  tryCatch({
    d <- data.frame(TIME_ = as.numeric(D[[time_col]]),
                    EV_   = as.integer(D[[status_col]] != 0),
                    G_    = factor(D[[strata_col]]))
    d <- d[is.finite(d$TIME_) & !is.na(d$G_), ]
    if (nlevels(droplevels(d$G_)) < 2) return(NULL)
    sd_ <- survival::survdiff(survival::Surv(TIME_, EV_) ~ G_, data = d)
    df_ <- length(sd_$n) - 1L
    tt <- data.table(
      comparison = "overall",
      test       = "log-rank (survdiff) on 1-KM curves across strata",
      statistic  = as.numeric(sd_$chisq),
      df         = df_,
      pvalue     = stats::pchisq(as.numeric(sd_$chisq), df = df_, lower.tail = FALSE),
      n_groups   = length(sd_$n)
    )
    list(tests = tt, pvalue_overall = tt$pvalue)
  }, error = function(e) {
    cat("  Log-rank stratum comparison failed:", e$message, "\n"); NULL
  })
}

bootstrap_risk_difference_safe <- function(data, indices, time_point,
                                           stratum_high, stratum_ref,
                                           time_col="TIME", status_col="PHENO", strata_col="STRATA") {
  d <- data[indices, ]
  dh <- d[get(strata_col)==stratum_high]
  dr <- d[get(strata_col)==stratum_ref]
  if (nrow(dh) < 2 || nrow(dr) < 2) return(NA_real_)
  out <- tryCatch({
    fh <- survfit(Surv(get(time_col), get(status_col)) ~ 1, data=dh)
    fr <- survfit(Surv(get(time_col), get(status_col)) ~ 1, data=dr)
    sh <- summary(fh, times=time_point, extend=TRUE)
    sr <- summary(fr, times=time_point, extend=TRUE)
    (1 - sh$surv) - (1 - sr$surv)
  }, error=function(e) NA_real_)
  out
}

calculate_rap <- function(surv_curves, stratum_high, stratum_ref, risk_level) {
  r  <- surv_curves[stratum==stratum_ref]
  h  <- surv_curves[stratum==stratum_high]
  if (nrow(r)==0 || nrow(h)==0) return(data.table(
    risk_level=risk_level, stratum_high=stratum_high, stratum_ref=stratum_ref,
    time_high=NA_real_, time_ref=NA_real_, rap_years=NA_real_))
  r[, ci := 1 - surv]; h[, ci := 1 - surv]
  tr <- r[ci >= risk_level]
  th <- h[ci >= risk_level]
  time_ref  <- if (nrow(tr)==0) NA_real_ else min(tr$time)
  time_high <- if (nrow(th)==0) NA_real_ else min(th$time)
  data.table(
    risk_level=risk_level, stratum_high=stratum_high, stratum_ref=stratum_ref,
    time_high=time_high, time_ref=time_ref,
    rap_years = if (is.na(time_ref) || is.na(time_high)) NA_real_ else (time_ref - time_high)
  )
}

calculate_paf <- function(absolute_risks, reference_stratum) {
  if (nrow(absolute_risks)==0) return(data.table())
  res <- absolute_risks[, .(
    pop_risk = weighted.mean(ci, w = pmax(n_risk, 1), na.rm=TRUE),
    ref_risk = ci[stratum==reference_stratum][1]
  ), by=time]
  res[, paf := ifelse(is.finite(pop_risk) & pop_risk>0 & is.finite(ref_risk),
                      (pop_risk - ref_risk)/pop_risk, NA_real_)]
  res
}

calculate_median_times <- function(DT, strata_col, time_col="TIME", status_col="PHENO") {
  if (!is.factor(DT[[strata_col]])) DT[[strata_col]] <- factor(DT[[strata_col]])
  lv <- levels(DT[[strata_col]])
  out <- vector("list", length(lv))
  for (i in seq_along(lv)) {
    st <- lv[i]
    df <- DT[get(strata_col)==st]
    if (nrow(df) < 2) next
    fit <- survfit(Surv(get(time_col), get(status_col)) ~ 1, data=df)
    sm  <- summary(fit)
    tab <- suppressWarnings(as.list(sm$table))
    out[[i]] <- data.table(
      stratum = st,
      median_time = suppressWarnings(as.numeric(tab$median)),
      median_lower = suppressWarnings(as.numeric(tab$`0.95LCL`)),
      median_upper = suppressWarnings(as.numeric(tab$`0.95UCL`))
    )
  }
  rbindlist(out, fill=TRUE)
}


# Predicted absolute risk now comes from the shared helper (prs_risk_utils.R).
#
# The previous implementation had two defects:
#   (1) survfit(fit, newdata=)$surv is a times x individuals matrix, but the
#       code took S[, idx] with idx a TIME index — returning one individual's
#       whole survival curve instead of every individual's risk at t. It raised
#       no error whenever n_times happened to equal n_individuals.
#   (2) it used a nearest-neighbour time lookup, which can select a time AFTER
#       the requested horizon, and indexed a start-stop model's baseline (age
#       scale) with a follow-up-scale horizon.
# This function feeds calculate_decision_curve(), so both defects propagated
# into net_benefit.csv and clinical_utility_metrics.csv.
predict_risk_at_time <- function(fit, newdata, t) {
  predict_risk_at_horizon(fit, newdata, t, start_col = "START")
}

# SURVIVAL DECISION CURVE (Vickers). The event probability among the
# threshold-positive group is estimated by Kaplan-Meier AT THE HORIZON, not by
# counting eventual events:
#
#   NB(p) = [1-S(t | risk>=p)] * P(risk>=p) - (p/(1-p)) * S(t | risk>=p) * P(risk>=p)
#
# Counting eventual `PHENO` instead (the previous implementation) treats an event
# occurring after t as a true positive at t, and someone censored before t as a
# true negative — so net benefit, sensitivity, specificity, PPV and NPV were not
# horizon-specific quantities at all.
.km_risk_at <- function(time, status, t) {
  if (!length(time)) return(NA_real_)
  s <- tryCatch(survival::survfit(survival::Surv(time, status) ~ 1), error=function(e) NULL)
  if (is.null(s) || !length(s$time)) return(0)
  idx <- findInterval(t, s$time)
  1 - if (idx == 0L) 1 else s$surv[idx]
}

calculate_decision_curve <- function(DT, t, fit, status_col="PHENO", time_col="TIME",
                                     high_risk_threshold=0.20) {
  pred_risk <- predict_risk_at_time(fit, DT, t)
  tt <- DT[[time_col]]; ss <- DT[[status_col]]
  # Unsupported horizons now yield NA risk (see prs_risk_utils.R); those rows
  # cannot be classified and must not silently count as low risk.
  usable   <- is.finite(pred_risk)
  n_unusable <- sum(!usable)
  if (n_unusable > 0)
    cat(sprintf("  ⚠ %d/%d individuals have NA predicted risk at t=%.2f (outside baseline support); excluded.\n",
                n_unusable, length(pred_risk), t))
  if (!any(usable)) return(NULL)
  pr <- pred_risk[usable]; tt <- tt[usable]; ss <- ss[usable]
  n  <- length(pr)

  # Event probability in the whole (usable) sample, by KM at the horizon.
  risk_all <- .km_risk_at(tt, ss, t)

  thr_seq <- seq(0.01, 0.99, by=0.01)
  nb_list <- vector("list", length(thr_seq))
  for (i in seq_along(thr_seq)) {
    thr <- thr_seq[i]
    w   <- thr/(1-thr)
    treat_all_nb <- risk_all - (1 - risk_all) * w
    hi  <- pr >= thr
    p_hi <- mean(hi)
    if (sum(hi) > 0) {
      risk_hi  <- .km_risk_at(tt[hi], ss[hi], t)
      nb_model <- risk_hi * p_hi - (1 - risk_hi) * p_hi * w
    } else { risk_hi <- NA_real_; nb_model <- 0 }
    nb_list[[i]] <- data.table(
      time_horizon=t, threshold=thr,
      nb_model=nb_model, nb_all=treat_all_nb, nb_none=0,
      km_risk_high_risk=risk_hi, km_risk_all=risk_all,
      n_high_risk=sum(hi), prop_high_risk=p_hi
    )
  }
  nb_dt <- rbindlist(nb_list)

  # Horizon-specific operating characteristics, all from KM at t.
  thr <- high_risk_threshold
  hi  <- pr >= thr
  p_hi <- mean(hi)
  if (sum(hi) == 0)
    cat(sprintf(paste0("  ⚠ No individual reaches the %.0f%% risk threshold at t=%.2f ",
                       "(max predicted risk %.4f); operating characteristics are NA.\n",
                       "    Lower --high_risk_threshold or use a longer horizon.\n"),
                100*thr, t, max(pr, na.rm=TRUE)))
  risk_hi <- if (sum(hi)  > 0) .km_risk_at(tt[hi],  ss[hi],  t) else NA_real_
  risk_lo <- if (sum(!hi) > 0) .km_risk_at(tt[!hi], ss[!hi], t) else NA_real_
  # Expected event / non-event mass, from KM rather than observed counts.
  tp <- risk_hi * p_hi
  fn <- risk_lo * (1 - p_hi)
  fp <- (1 - risk_hi) * p_hi
  tn <- (1 - risk_lo) * (1 - p_hi)
  sens <- if (is.finite(tp + fn) && (tp + fn) > 0) tp/(tp + fn) else NA_real_
  spec <- if (is.finite(tn + fp) && (tn + fp) > 0) tn/(tn + fp) else NA_real_
  ppv  <- if (is.finite(risk_hi)) risk_hi else NA_real_
  npv  <- if (is.finite(risk_lo)) 1 - risk_lo else NA_real_
  row  <- nb_dt[threshold==thr]
  metrics <- data.table(
    time_horizon=t, threshold=thr, sensitivity=sens, specificity=spec,
    ppv=ppv, npv=npv,
    km_risk_high_risk=risk_hi, km_risk_low_risk=risk_lo, km_risk_all=risk_all,
    n_excluded_unsupported=n_unusable,
    estimator="KM-based survival decision curve (Vickers); apparent, same-sample",
    nb_model=if (nrow(row)) row$nb_model else NA_real_,
    nb_all=if (nrow(row)) row$nb_all else NA_real_,
    prop_high_risk=if (nrow(row)) row$prop_high_risk else p_hi
  )
  list(net_benefit=nb_dt, metrics=metrics, predicted_risk=pred_risk)
}

# ==============================================================================
# MAIN
# ==============================================================================
cat("==============================================================================\n")
cat("CUMULATIVE INCIDENCE ANALYSIS\n")
cat("==============================================================================\n\n")

all_results <- list()

# Discover keys from quantile_list or DT
keys <- names(quantile_list)
if (is.null(keys) || length(keys)==0) {
  cand_cols <- grep("PRS|Q[0-9]+|quantile|extreme|group", names(DT), value=TRUE, ignore.case=TRUE)
  keys <- unique(c("Q5","Q10","quantile_extreme_1","quantile_extreme_5", cand_cols))
}

for (q_name in keys) {
  cat("Analyzing key:", q_name, "\n")
  cat("----------------------------------------------------------------------\n")
  
  strata_col_raw <- resolve_strata_column(DT, quantile_list, q_name)
  if (is.na(strata_col_raw) || !strata_col_raw %in% names(DT)) {
    cat("  Error: No valid strata column found for key '", q_name, "'.\n\n", sep="")
    all_results[[q_name]] <- NULL
    next
  }
  cat("  Using strata column:", strata_col_raw, "\n")
  
  if (!is.factor(DT[[strata_col_raw]])) DT[[strata_col_raw]] <- factor(DT[[strata_col_raw]])
  strata_levels <- levels(DT[[strata_col_raw]])
  cat("  Strata levels:", paste(strata_levels, collapse=", "), "\n")
  
  # 1) Survival curves by stratum
  cat("  1. Calculating cumulative incidence curves...\n")
  surv_curves <- get_strata_survival(DT, strata_col_raw, timeCol, statusCol)
  if (nrow(surv_curves)==0) {
    cat("  Warning: no survival data for this key; skipping.\n\n")
    all_results[[q_name]] <- NULL
    next
  }
  # max_time is a DISPLAY bound (95th pct by default). It must NOT truncate the
  # curve used for estimation, or a valid horizon past the 95th percentile is
  # wrongly reported unsupported. Keep the full curve for every numeric
  # calculation and the (private) dense export; truncate only for plotting.
  surv_curves[, `:=`(
    ci = 1 - surv,
    ci_lower = 1 - surv_upper,
    ci_upper = 1 - surv_lower
  )]
  surv_curves_plot <- surv_curves[time <= opt$max_time]
  
  # 2) Absolute risks at target times (pre-filtered times from 01*R)
  cat("  2. Estimating absolute risks at target times...\n")
  absolute_risks <- calculate_absolute_risks(surv_curves, time_points, min_cell_count = opt$min_cell_count)
  
  if (nrow(absolute_risks)>0) {
    cat("    Absolute risk ranges:\n")
    for (tp in time_points) {
      rr <- absolute_risks[time==tp]
      if (nrow(rr)>0) {
        r <- range(rr$ci*100, na.rm=TRUE)
        cat(sprintf("      %g years: %.1f%% - %.1f%%\n", tp, r[1], r[2]))
      }
    }
  }
  
  # 3) Risk differences vs reference
  ref_label <- opt$reference_quantile
  if (!ref_label %in% strata_levels) {
    cat("  Note: reference label '", ref_label, "' not present; using '", strata_levels[1], "'.\n", sep="")
    ref_label <- strata_levels[1]
  }
  cat("  Reference stratum:", ref_label, "\n")
  risk_differences <- calculate_risk_differences(absolute_risks, ref_label)
  
  # 4) NNS
  cat("  3. Calculating number needed to screen (NNS)...\n")
  # NNS is a clinical-utility output: computed only when opted in. Inverting a
  # between-quantile risk difference is not a number needed to screen unless a
  # screening action, comparator and causal effect are defined.
  nns <- if (isTRUE(opt$enable_clinical_utility)) calculate_nns(risk_differences) else data.table()
  
  # 5) Plots (optional)
  if (opt$generate_plots) {
    cat("  4. Generating cumulative incidence plots...\n")
    p_ci <- ggplot(surv_curves_plot, aes(x=time, y=ci, color=stratum, fill=stratum)) +
      geom_line(linewidth=1) +
      geom_ribbon(aes(ymin=ci_lower, ymax=ci_upper), alpha=0.2, color=NA) +
      scale_y_continuous(labels=scales::percent_format(accuracy=1), limits=c(0,1)) +
      scale_x_continuous(breaks=seq(0, max(surv_curves$time, na.rm=TRUE), by=5)) +
      labs(title=paste0("Cumulative risk (1-KM) by ", q_name),
           x="Time (years)", y="Cumulative risk (1-KM)",
           color="Stratum", fill="Stratum") +
      theme_bw(base_size=12) +
      theme(legend.position="right", panel.grid.minor=element_blank(),
            plot.title=element_text(hjust=0.5, face="bold"))
    out_ci <- file.path(opt$outdir, paste0(prefix, "cumulative_incidence_", q_name, ".png"))
    save_plot_cairo(p_ci, out_ci, width=opt$figure_width, height=opt$figure_height, dpi=opt$dpi, formats=trimws(strsplit(opt$plot_formats, ",")[[1]]))
    cat("    ✓ Saved:", basename(out_ci), "\n")
  } else {
    cat("  4. Skipping plots (--generate_plots not set)\n")
  }
  
  # 6) Risk table (wide)
  cat("  5. Building risk table...\n")
  risk_table <- dcast(absolute_risks, stratum ~ time, value.var="ci")
  risk_table[, N := sapply(stratum, function(s) sum(DT[[strata_col_raw]]==s))]
  # Disclosure: blank N and every risk cell for strata below the threshold.
  small <- risk_table$N < opt$min_cell_count
  risk_table[, minimum_cell_count_pass := !small]
  if (any(small)) {
    time_cols <- setdiff(names(risk_table), c("stratum","N","minimum_cell_count_pass"))
    for (col in time_cols) risk_table[small, (col) := NA_real_]
    risk_table[small, N := NA_integer_]
  }
  setcolorder(risk_table, c("stratum","N","minimum_cell_count_pass",
                            setdiff(names(risk_table), c("stratum","N","minimum_cell_count_pass"))))
  for (col in names(risk_table)) if (!col %in% c("stratum","N","minimum_cell_count_pass"))
    risk_table[[col]] <- ifelse(is.na(risk_table[[col]]), NA_character_, sprintf("%.1f%%", 100*risk_table[[col]]))
  
  # Store
  all_results[[q_name]] <- list(
    surv_curves = surv_curves,
    absolute_risks = absolute_risks,
    risk_differences = risk_differences,
    nns = nns,
    risk_table = risk_table,
    strata_col = strata_col_raw,
    reference = ref_label
  )
  cat("\n")
}

# ==============================================================================
# STRATUM COMPARISON (log-rank on 1-KM curves; NOT Gray's test — single event type)
# ==============================================================================
cat("==============================================================================\n")
cat("STRATUM COMPARISON (log-rank on 1-KM curves)\n")
cat("==============================================================================\n\n")

gray_results <- list()
for (q_name in names(all_results)) {
  if (is.null(all_results[[q_name]])) next
  strata_col_raw <- all_results[[q_name]]$strata_col
  cat("Testing", q_name, "...\n")
  Dclean <- sanitize_surv_data(DT, strata_col=strata_col_raw, time_col=timeCol, status_col=statusCol)
  gt <- perform_logrank_test(Dclean, strata_col="STRATA", time_col="TIME", status_col="PHENO")
  if (!is.null(gt)) {
    cat("  Overall p-value:", if (!is.null(gt$pvalue_overall)) format.pval(gt$pvalue_overall) else "NA", "\n\n")
    gray_results[[q_name]] <- gt
  } else {
    cat("  Gray's test not available.\n\n")
  }
}

# ==============================================================================
# BOOTSTRAP CIs FOR RISK DIFFERENCES
# ==============================================================================
cat("==============================================================================\n")
cat("BOOTSTRAP CONFIDENCE INTERVALS\n")
cat("==============================================================================\n\n")

# OFF BY DEFAULT. This resamples survfit for every
# (PRS role x quantile scheme x non-reference group x horizon) combination at
# n_bootstrap iterations — hundreds of thousands of survfit() calls per trajectory
# — and it ran regardless of --enable_clinical_utility. The risk-difference POINT
# estimates and their KM-derived bounds are always produced; only the resampling
# CIs are gated. Enable with --bootstrap_risk_differences.
bootstrap_ci_results <- list()
bootstrap_ci_dt <- data.table()
if (isTRUE(opt$bootstrap_risk_differences) && opt$n_bootstrap > 0) {
cat("Computing bootstrap CIs for risk differences (", opt$n_bootstrap, " iterations)...\n\n", sep="")

for (q_name in names(all_results)) {
  if (is.null(all_results[[q_name]])) next
  strata_col_raw <- all_results[[q_name]]$strata_col
  ref_label      <- all_results[[q_name]]$reference
  
  Dclean <- sanitize_surv_data(DT, strata_col=strata_col_raw, time_col=timeCol, status_col=statusCol)
  lv <- levels(Dclean$STRATA)
  if (!ref_label %in% lv) {
    cat("  Skipping", q_name, "- reference label not present in cleaned data.\n")
    next
  }
  
  for (st in lv) {
    if (st == ref_label) next
    for (tp in time_points) {
      cat(sprintf("  Bootstrapping %s : %s vs %s @ %g years...\r", q_name, st, ref_label, tp))
      br <- boot(
        data = Dclean[, .(STRATA, TIME, PHENO)],
        statistic = function(dat, ind) bootstrap_risk_difference_safe(
          data = dat, indices = ind, time_point = tp,
          stratum_high = st, stratum_ref = ref_label,
          time_col = "TIME", status_col = "PHENO", strata_col = "STRATA"
        ),
        R = opt$n_bootstrap
      )
      tvals <- br$t[is.finite(br$t)]
      if (length(tvals) < max(50, 0.2*opt$n_bootstrap)) {
        ci_low <- NA_real_; ci_up <- NA_real_; est <- mean(tvals)
      } else {
        est <- mean(tvals)
        ci  <- tryCatch(boot.ci(br, type="perc"), error=function(e) NULL)
        if (!is.null(ci) && !is.null(ci$percent)) {
          ci_low <- ci$percent[4]; ci_up <- ci$percent[5]
        } else { ci_low <- NA_real_; ci_up <- NA_real_ }
      }
      bootstrap_ci_results[[length(bootstrap_ci_results)+1]] <- data.table(
        quantile_type = q_name, stratum = st, time = tp,
        risk_diff = est, ci_lower = ci_low, ci_upper = ci_up,
        n_boot = length(tvals)
      )
    }
  }
}
cat("\n")
bootstrap_ci_dt <- if (length(bootstrap_ci_results)>0) rbindlist(bootstrap_ci_results, fill=TRUE) else data.table()
cat("✓ Bootstrap confidence intervals computed\n\n")
} else {
  cat("Skipping risk-difference bootstrap (default). Use --bootstrap_risk_differences to enable.\n\n")
}

# ==============================================================================
# RAP
# ==============================================================================
cat("==============================================================================\n")
cat("RISK ADVANCEMENT PERIOD (RAP)\n")
cat("==============================================================================\n\n")

rap_results <- list()
for (q_name in names(all_results)) {
  if (is.null(all_results[[q_name]])) next
  surv_curves <- all_results[[q_name]]$surv_curves
  lv <- sort(unique(surv_curves$stratum))
  if (length(lv) < 2) next
  high <- lv[length(lv)]
  ref  <- all_results[[q_name]]$reference
  for (rl in rap_thresholds) {
    rap <- calculate_rap(surv_curves, stratum_high=high, stratum_ref=ref, risk_level=rl)
    rap[, quantile_type := q_name]
    rap_results[[length(rap_results)+1]] <- rap
    if (!is.na(rap$rap_years)) {
      cat(sprintf("  %s @ %.0f%%: RAP = %.1f years (high vs %s)\n",
                  q_name, rl*100, rap$rap_years, ref))
    }
  }
}
rap_dt <- if (length(rap_results)>0) rbindlist(rap_results, fill=TRUE) else data.table()
cat("\n✓ Risk advancement period calculated\n\n")

# ==============================================================================
# PAF
# ==============================================================================
cat("==============================================================================\n")
cat("POPULATION ATTRIBUTABLE FRACTION (PAF)\n")
cat("==============================================================================\n\n")

paf_results <- list()
# Population attributable fraction is clinical-utility only: it presumes an
# intervention that removes the excess risk, which this pipeline does not define.
if (isTRUE(opt$enable_clinical_utility)) for (q_name in names(all_results)) {
  if (is.null(all_results[[q_name]])) next
  abs_r <- all_results[[q_name]]$absolute_risks
  ref   <- all_results[[q_name]]$reference
  paf   <- calculate_paf(abs_r, ref)
  paf[, quantile_type := q_name]
  paf_results[[length(paf_results)+1]] <- paf
}
paf_dt <- if (length(paf_results)>0) rbindlist(paf_results, fill=TRUE) else data.table()
if (nrow(paf_dt)>0) {
  cat("  Example PAF at last time:\n")
  lastt <- max(paf_dt$time)
  ex <- paf_dt[time==lastt][1]
  if (nrow(ex)) cat(sprintf("  %g years: PAF = %.1f%% (Pop %.1f%%, Ref %.1f%%)\n",
                            lastt, ex$paf*100, ex$pop_risk*100, ex$ref_risk*100))
}
cat("\n✓ Population attributable fraction calculated\n\n")

# ==============================================================================
# MEDIAN TIME TO EVENT
# ==============================================================================
cat("==============================================================================\n")
cat("MEDIAN TIME TO EVENT\n")
cat("==============================================================================\n\n")

median_results <- list()
for (q_name in names(all_results)) {
  if (is.null(all_results[[q_name]])) next
  sc <- all_results[[q_name]]$strata_col
  med <- calculate_median_times(DT, sc, timeCol, statusCol)
  med[, quantile_type := q_name]
  median_results[[length(median_results)+1]] <- med
}
median_dt <- if (length(median_results)>0) rbindlist(median_results, fill=TRUE) else data.table()
cat("✓ Median time to event calculated\n\n")

# ==============================================================================
# RISK PER SD OF PRS - REMOVED
# ==============================================================================
# The former risk_per_sd output merged a GROUP-level KM risk onto individuals and
# regressed it on individual PRS (lm(ci ~ PRS)). That is pseudo-replication: every
# member of a quantile shares one outcome value, so the SE and p-value were
# meaningless, and it always used PRS_std (the ONSET score) even when the strata
# came from another PRS role. The correctly estimated per-SD effect is the
# continuous-PRS Cox coefficient, already reported in *_all_coefficients.csv.


# ==============================================================================
# DECISION CURVES
# ==============================================================================
cat("==============================================================================\n")
cat("CLINICAL UTILITY ANALYSIS (Decision Curves)\n")
cat("==============================================================================\n\n")

# Model selection. Previously hardcoded to "minimal", which contains only the
# PRIMARY PRS (the onset score as the driver supplies it) — so the decision
# curves described M1, not the manuscript's M5/M2 comparison.
.resolve_risk_model <- function(want) {
  mm <- results$metadata$model_map
  if (!is.null(want) && nzchar(want)) {
    if (want %in% names(fitted_models)) return(want)
    if (!is.null(mm)) {
      hit <- as.data.table(mm)[model_id == want, model]
      if (length(hit) && hit[1] %in% names(fitted_models)) return(as.character(hit[1]))
    }
    stop("--risk_model '", want, "' not found. Available: ",
         paste(names(fitted_models), collapse = ", "))
  }
  # Default: manuscript M5 (outcome + progression) if present, else legacy behaviour.
  if (!is.null(mm)) {
    hit <- as.data.table(mm)[model_id == "M5", model]
    if (length(hit) && hit[1] %in% names(fitted_models)) return(as.character(hit[1]))
  }
  if ("minimal" %in% names(fitted_models)) "minimal" else names(fitted_models)[1]
}
minimal_model_name <- .resolve_risk_model(opt$risk_model)
fit <- fitted_models[[minimal_model_name]]
cat("Risk model for decision curves / predicted risk:", minimal_model_name, "\n\n")

# Per-stratum baseline support at each evaluated horizon. Beyond the last observed
# baseline event time the conditional survival ratio collapses to 1, which would
# report an UNSUPPORTED horizon as exactly zero risk; the helper now returns NA
# there, and this table says which strata/horizons are affected. In a stratified
# prospective model one T1_STATUS stratum can be unsupported while the pooled
# sample looks fine, so support is reported per stratum.
hs_all <- tryCatch({
  rbindlist(lapply(time_points, function(tp) {
    h <- horizon_support(fit, DT, tp)
    h$horizon <- tp
    h$model   <- minimal_model_name
    h
  }), fill = TRUE)
}, error = function(e) { cat("  Note: horizon support unavailable:", conditionMessage(e), "\n"); NULL })
if (!is.null(hs_all) && nrow(hs_all)) {
  # Disclosure: per-stratum support counts can be small in a thin stratum.
  suppress_small_cells(hs_all,
                       count_cols = intersect(c("n","n_at_risk_at_end","n_events_by_end"), names(hs_all)),
                       threshold = opt$min_cell_count,
                       estimate_cols = intersect(c("censoring_survival"), names(hs_all)))
  fwrite(hs_all, file.path(opt$outdir, paste0(prefix, "horizon_support.csv")))
  cat("✓ Saved horizon_support.csv\n")
  bad <- hs_all[supported == FALSE]
  if (nrow(bad))
    cat("  ⚠ Unsupported horizon/stratum combinations (risks are NA, not 0):\n",
        paste0("      ", bad$stratum, " @ ", bad$horizon, "y\n", collapse = ""), sep = "")
  cat("\n")
}

# Decision curves are clinical-utility only. horizon_support above stays
# unconditional: it is a QC diagnostic for the descriptive risks, not a utility metric.
clinical_utility <- NULL
if (isTRUE(opt$enable_clinical_utility)) {
if (opt$decision_curves_all_times) {
  cat("Calculating decision curves at all time points:", paste(time_points, collapse=", "), "\n\n")
  nb_all <- list(); met_all <- list()
  for (tp in time_points) {
    dc <- calculate_decision_curve(DT, tp, fit, status_col=statusCol, high_risk_threshold=opt$high_risk_threshold)
    nb_all[[length(nb_all)+1]] <- dc$net_benefit
    met_all[[length(met_all)+1]] <- dc$metrics
  }
  net_benefit_all_times <- rbindlist(nb_all, fill=TRUE)
  metrics_all_times     <- rbindlist(met_all, fill=TRUE)
  if (opt$generate_plots) {
    for (tp in time_points) {
      tp_data <- net_benefit_all_times[time_horizon==tp]
      p <- ggplot(tp_data, aes(x=threshold)) +
        geom_line(aes(y=nb_model, color="PRS Model"), linewidth=1) +
        geom_line(aes(y=nb_all,   color="Treat All"), linewidth=1, linetype="dashed") +
        geom_line(aes(y=nb_none,  color="Treat None"), linewidth=1, linetype="dotted") +
        scale_color_manual(values=c("PRS Model"="#d7191c", "Treat All"="#2c7bb6", "Treat None"="black")) +
        scale_x_continuous(labels=scales::percent_format(accuracy=1)) +
        labs(title=paste0("Decision Curve (", tp, "-year)"), x="Risk Threshold", y="Net Benefit", color="Strategy") +
        theme_bw(base_size=12) + theme(legend.position="right", panel.grid.minor=element_blank(),
                                       plot.title=element_text(hjust=0.5, face="bold"))
      f <- file.path(opt$outdir, paste0(prefix, "decision_curve_", tp, "yr.png"))
      save_plot_cairo(p, f, width=opt$figure_width, height=opt$figure_height*0.75, dpi=opt$dpi, formats=trimws(strsplit(opt$plot_formats, ",")[[1]]))
      cat("  ✓", basename(f), "\n")
    }
  }
  clinical_utility <- list(net_benefit=net_benefit_all_times, metrics_at_threshold=metrics_all_times, time_points=time_points)
} else {
  # Use the LONGEST available horizon, not the shortest. At the shortest horizon
  # predicted risks are far below any clinically plausible threshold, so the
  # threshold-positive group is empty and every operating characteristic is NA.
  # (The old eventual-PHENO counting masked this by returning numbers that were
  # not horizon-specific.)
  target_time <- time_points[length(time_points)]
  cat("Target time horizon:", target_time, "years (longest available; use --decision_curves_all_times for all)\n\n")
  dc <- calculate_decision_curve(DT, target_time, fit, status_col=statusCol, high_risk_threshold=opt$high_risk_threshold)
  net_benefit_dt <- dc$net_benefit
  if (opt$generate_plots) {
    p <- ggplot(net_benefit_dt, aes(x=threshold)) +
      geom_line(aes(y=nb_model, color="PRS Model"), linewidth=1) +
      geom_line(aes(y=nb_all,   color="Treat All"), linewidth=1, linetype="dashed") +
      geom_line(aes(y=nb_none,  color="Treat None"), linewidth=1, linetype="dotted") +
      scale_color_manual(values=c("PRS Model"="#d7191c", "Treat All"="#2c7bb6", "Treat None"="black")) +
      scale_x_continuous(labels=scales::percent_format(accuracy=1)) +
      labs(title=paste0("Decision Curve (", target_time, "-year)"), x="Risk Threshold", y="Net Benefit", color="Strategy") +
      theme_bw(base_size=12) + theme(legend.position="right", panel.grid.minor=element_blank(),
                                     plot.title=element_text(hjust=0.5, face="bold"))
    f <- file.path(opt$outdir, paste0(prefix, "decision_curve.png"))
    save_plot_cairo(p, f, width=opt$figure_width, height=opt$figure_height*0.75, dpi=opt$dpi, formats=trimws(strsplit(opt$plot_formats, ",")[[1]]))
    cat("  ✓", basename(f), "\n")
  }
  clinical_utility <- list(net_benefit=net_benefit_dt, metrics_at_threshold=dc$metrics, time_points=target_time)
}
}  # end enable_clinical_utility (decision curves)

# ==============================================================================
# LIFETIME RISK
# ==============================================================================
cat("==============================================================================\n")
cat("LIFETIME RISK ESTIMATION\n")
cat("==============================================================================\n\n")

lifetime_results <- list()
# Maximum observed cumulative risk (NOT lifetime risk) - clinical-utility only.
if (isTRUE(opt$enable_clinical_utility)) for (q_name in names(all_results)) {
  if (is.null(all_results[[q_name]])) next
  sc <- all_results[[q_name]]$surv_curves
  lr <- sc[, .(max_observed_cumulative_risk = max(ci, na.rm=TRUE), max_time = max(time, na.rm=TRUE)), by=stratum]
  lr[, quantile := q_name]
  lifetime_results[[q_name]] <- lr
  print(lr[, .(stratum, max_observed_cumulative_risk=sprintf("%.1f%%", max_observed_cumulative_risk*100),
               max_time=sprintf("%.1f years", max_time))])
  cat("\n")
}
lifetime_dt <- rbindlist(lifetime_results, fill=TRUE)

# ==============================================================================
# SAVE RESULTS
# ==============================================================================
cat("==============================================================================\n")
cat("SAVING RESULTS\n")
cat("==============================================================================\n\n")

for (q_name in names(all_results)) {
  if (is.null(all_results[[q_name]])) next
  fwrite(all_results[[q_name]]$surv_curves,     file.path(opt$outdir, paste0(prefix, "ci_curves_", q_name, ".csv")))
  fwrite(all_results[[q_name]]$absolute_risks,  file.path(opt$outdir, paste0(prefix, "absolute_risks_", q_name, ".csv")))
  fwrite(all_results[[q_name]]$risk_differences,file.path(opt$outdir, paste0(prefix, "risk_differences_", q_name, ".csv")))
  if (opt$enable_clinical_utility &&
      !is.null(all_results[[q_name]]$nns) && nrow(all_results[[q_name]]$nns)>0) {
    fwrite(all_results[[q_name]]$nns,           file.path(opt$outdir, paste0(prefix, "nns_", q_name, ".csv")))
  }
  fwrite(all_results[[q_name]]$risk_table,      file.path(opt$outdir, paste0(prefix, "risk_table_", q_name, ".csv")))
}
cat("✓ Saved cumulative incidence results (1-KM cumulative risk under independent censoring)\n")

if (length(gray_results)>0) {
  # NOT Gray's test: Gray's test compares cumulative INCIDENCE functions under
  # competing risks, and this pipeline models a single event type with death not
  # treated as a competing event. What is computed is a log-rank-family comparison
  # of 1-KM curves across PRS strata, so the file is named for what it is.
  for (q_name in names(gray_results)) if (!is.null(gray_results[[q_name]])) {
    fwrite(gray_results[[q_name]]$tests,        file.path(opt$outdir, paste0(prefix, "strata_comparison_test_", q_name, ".csv")))
  }
  cat("✓ Saved stratum comparison tests (log-rank family on 1-KM curves)\n")
}

if (nrow(bootstrap_ci_dt)>0) {
  fwrite(bootstrap_ci_dt,                       file.path(opt$outdir, paste0(prefix, "bootstrap_ci.csv")))
  cat("✓ Saved bootstrap confidence intervals\n")
}

if (nrow(rap_dt)>0) {
  fwrite(rap_dt,                                file.path(opt$outdir, paste0(prefix, "risk_advancement_period.csv")))
  cat("✓ Saved risk advancement period\n")
}

if (opt$enable_clinical_utility && nrow(paf_dt)>0) {
  fwrite(paf_dt,                                file.path(opt$outdir, paste0(prefix, "population_attributable_fraction.csv")))
  cat("✓ Saved population attributable fraction (opt-in)\n")
}

if (nrow(median_dt)>0) {
  fwrite(median_dt,                             file.path(opt$outdir, paste0(prefix, "median_times.csv")))
  cat("✓ Saved median time to event\n")
}

if (opt$enable_clinical_utility && !is.null(clinical_utility)) {
  fwrite(clinical_utility$net_benefit,         file.path(opt$outdir, paste0(prefix, "net_benefit.csv")))
  fwrite(clinical_utility$metrics_at_threshold,file.path(opt$outdir, paste0(prefix, "clinical_utility_metrics.csv")))
  cat("✓ Saved clinical utility results (opt-in; KM-based survival DCA, APPARENT performance)\n")
}

if (opt$enable_clinical_utility) {
  fwrite(lifetime_dt,                            file.path(opt$outdir, paste0(prefix, "max_observed_cumulative_risk.csv")))
  cat("✓ Saved max_observed_cumulative_risk.csv (opt-in; max observed 1-KM risk at longest follow-up, NOT a lifetime projection)\n")
}

# Save a comprehensive RDS
all_outputs <- list(
  cumulative_incidence = all_results,
  statistical_tests = list(gray_test = gray_results, bootstrap_ci = bootstrap_ci_dt),
  advanced_metrics = list(risk_advancement_period = rap_dt, population_attributable_fraction = paf_dt,
                          median_times = median_dt),
  clinical_utility = clinical_utility,
  lifetime_risks = lifetime_dt,
  settings = list(
    time_points = time_points,
    max_time = opt$max_time,
    reference_quantile = opt$reference_quantile,
    high_risk_threshold = opt$high_risk_threshold,
    n_bootstrap = opt$n_bootstrap,
    bootstrap_seed = opt$bootstrap_seed,
    rap_thresholds = rap_thresholds,
    generate_plots = opt$generate_plots,
    decision_curves_all_times = opt$decision_curves_all_times,
    used_saved_timepoints = isTRUE(opt$use_supported_times)
  ),
  # also surface where times came from
  consumed_timepoints_source = list(
    from_requested_timepoints_ghat_csv = !is.null(read_saved_requested_times(opt$models_file)),
    from_metadata_alias = !is.null(find_supported_times_from_metadata(results, DT))
  ),
  metadata = metadata
)
saveRDS(all_outputs, file.path(opt$outdir, paste0(prefix, "results.rds")))
cat("✓ Saved comprehensive results object\n\n")

# ==============================================================================
# SUMMARY
# ==============================================================================
cat("==============================================================================\n")
cat("SUMMARY\n")
cat("==============================================================================\n\n")
cat("✓ Analysis completed!\n\n")
cat("Outputs saved to:", opt$outdir, "\n")
cat("Log file:", log_file, "\n")

sink(type="message"); sink(type="output"); close(log_con)
