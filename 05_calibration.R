#!/usr/bin/env Rscript
# ==============================================================================
# 05_calibration.R - Calibration Analysis for Survival Models
#   - All Cox models from 01_fit_models.R (base, minimal, full, quantiles, extremes, ...)
#   - Metrics at chosen times: slope, intercept, R², ICI, E-statistic
#   - Point metrics: CIL, O/E, obs_minus_pred_q50, obs_minus_pred_q90
#   - Grouped calibration (Nam–D’Agostino–Greenwood) by risk quantiles
#   - Recalibration modes (applied per model × time):
#       raw, recal_int, recal_slope, recal_int_slope, platt, isotonic
#   - Optional bootstrap ribbons for grouped curves
#   - Optional plots; PNG default
#
# Inputs:
#   --models_file  : fitted_models.rds from 01_fit_models.R
#   --data_file    : data_processed.rds from 01_fit_models.R
#
# Outputs (key):
#   <prefix>05_calibration_calibration_metrics.csv
#   <prefix>05_calibration_calibration_groups.csv
#   <prefix>05_calibration_calibration_plot_points.csv
#   <prefix>05_calibration_calibration_plot_bands.csv
#   <prefix>05_calibration_calibration_by_time.csv
#   <prefix>05_calibration_cal_<model>_t<time>_<mode>.(png|pdf)
#   <prefix>05_calibration_slope_trajectory.(png|pdf)
#   <prefix>05_calibration_ici_trajectory.(png|pdf)
#   <prefix>05_calibration_log.txt
#   <prefix>05_calibration_cal_times_ghat_audit.csv   # (NEW) audit of calibration times and guards
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

required_packages <- c("optparse", "data.table", "survival", "ggplot2", "rms", "pec")
invisible(lapply(required_packages, load_package))  # Cairo intentionally omitted (optional; see prs_risk_utils.R)

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
# COMMAND-LINE ARGUMENTS
# ==============================================================================
option_list <- list(
  # Input files
  make_option("--models_file", type="character", default=NULL,
              help="Path to fitted_models.rds from 01_fit_models.R [REQUIRED]"),
  make_option("--data_file", type="character", default=NULL,
              help="Path to data_processed.rds from 01_fit_models.R [REQUIRED]"),
  
  # (NEW) Use time points saved by 01* (metadata$time_support$eval_time_points)
  make_option("--use_supported_times", action="store_true", default=FALSE,
              help="Use eval_time_points saved by 01* (falls back to requested_valid_time_points, then t_pass_auto)."),
  
  # Calibration time points (fallback if not using supported times)
  make_option("--cal_times", type="character", default=NULL,
              help="Comma-separated calibration time points in years (e.g., '5,10,15') [default: auto/IQR]"),
  make_option("--n_cal_times", type="integer", default=3,
              help="Number of calibration time points if auto-selected [default: 3]"),
  
  # (NEW) Guard thresholds (applied to chosen times)
  make_option("--min_events", type="integer", default=10,
              help="Minimum cumulative events by t to keep t [default: 10]"),
  make_option("--min_at_risk", type="integer", default=50,
              help="Minimum at-risk at t to keep t [default: 50]"),
  make_option("--ghat_threshold", type="numeric", default=0.02,
              help="Minimum G-hat (censoring survival) at t to keep t [default: 0.02]"),
  
  # Grouping / smoothing
  make_option("--calibration_modes", type="character", default="raw",
              help=paste("Comma-separated recalibration modes [default: raw]. Options:",
                         "raw,recal_int,recal_slope,recal_int_slope,platt,isotonic. The non-raw",
                         "modes are trained and evaluated in the same sample (optimistic).")),
  make_option("--calibration_models", type="character", default=NULL,
              help=paste("Comma-separated models to calibrate (manuscript ids M0-M7 or raw names).",
                         "Default: the core M0-M7 ladder. Clinical models are never calibrated here",
                         "(they live on a different sample).")),
  make_option("--min_cell_count", type="integer", default=10,
              help="Disclosure: drop calibration groups smaller than this [default: 10]"),
  make_option("--n_groups", type="integer", default=10,
              help="Number of quantile risk groups for grouped calibration [default: 10]"),
  make_option("--smooth_calibration", action="store_true", default=TRUE,
              help="Add LOESS smoothed curve on plots [default: TRUE]"),
  
  # Stratified calibration (optional; must exist in data)
  make_option("--stratify_by", type="character", default=NULL,
              help="Variable for stratified calibration (e.g., 'sex', 'age_strata') [default: NULL]"),
  
  # Bootstrap bands
  make_option("--bootstrap_bands", action="store_true", default=FALSE,
              help="Compute bootstrap ribbons for grouped curves [default: FALSE]"),
  make_option("--n_boot", type="integer", default=500,
              help="Bootstrap replicates for ribbons [default: 500]"),
  make_option("--boot_seed", type="integer", default=123,
              help="Seed for bootstrap resampling [default: 123]"),
  make_option("--ci_level", type="numeric", default=0.95,
              help="Confidence level for ribbons [default: 0.95]"),
  
  # Plot options
  make_option("--make_plots", action="store_true", default=FALSE,
              help="Generate plots [default: FALSE]"),
  make_option("--format", type="character", default="png",
              help="Output formats: 'png', 'pdf', or 'pdf,png' [default: png]"),
  make_option("--plot_width", type="numeric", default=8,
              help="Plot width [inches] [default: 8]"),
  make_option("--plot_height", type="numeric", default=8,
              help="Plot height [inches] [default: 8]"),
  make_option("--dpi", type="numeric", default=300,
              help="DPI for PNG output [default: 300]"),
  
  # Output & logging
  make_option("--outdir", type="character", default="output",
              help="Output directory [default: output]"),
  make_option("--prefix", type="character", default="",
              help="Prefix for output files [default: auto]"),
  make_option("--verbose", action="store_true", default=TRUE,
              help="Print detailed output [default: TRUE]")
)

opt <- parse_args(OptionParser(
  option_list=option_list,
  description="\nCalibration of survival models from 01_fit_models.R with raw + multiple recalibration modes.\n",
  epilogue="
Examples:
  Rscript 05_calibration.R \\
    --models_file results/fitted_models.rds \\
    --data_file   results/data_processed.rds \\
    --use_supported_times \\
    --make_plots --bootstrap_bands --n_boot 1000
"
))

# Validate
if (is.null(opt$models_file)) stop("--models_file is required")
if (is.null(opt$data_file))   stop("--data_file is required")
if (!file.exists(opt$models_file)) stop("Models file not found: ", opt$models_file)
if (!file.exists(opt$data_file))   stop("Data file not found: ", opt$data_file)

# Create output
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

cat("==============================================================================\n")
cat("CALIBRATION ANALYSIS\n")
cat("==============================================================================\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("==============================================================================\n\n")

# ==============================================================================
# LOAD DATA
# ==============================================================================
cat("LOADING DATA\n")
cat("----------------------------------------------------------------------\n")
results <- readRDS(opt$models_file)
DT      <- readRDS(opt$data_file)

cat("✓ Loaded fitted models\n")
cat("✓ Loaded processed data (N=", nrow(DT), ")\n\n", sep="")

metadata  <- results$metadata
timeCol   <- "TIME"
statusCol <- "PHENO"

# Construct prefix
if (nchar(opt$prefix) > 0) {
  prefix <- paste0(opt$prefix, "_05_calibration_")
} else {
  base_prefix_parts <- c()
  if (!is.null(metadata$options$cohort)    && metadata$options$cohort    != "") base_prefix_parts <- c(base_prefix_parts, metadata$options$cohort)
  if (!is.null(metadata$options$ancestry)  && metadata$options$ancestry  != "") base_prefix_parts <- c(base_prefix_parts, metadata$options$ancestry)
  if (!is.null(metadata$options$pheno_name)&& metadata$options$pheno_name!= "") base_prefix_parts <- c(base_prefix_parts, metadata$options$pheno_name)
  prefix <- if (length(base_prefix_parts) > 0) paste0(paste(base_prefix_parts, collapse="_"), "_05_calibration_") else "05_calibration_"
}
cat("Output prefix:", prefix, "\n\n")

# Logging
log_file <- file.path(opt$outdir, paste0(prefix, "log.txt"))
log_con  <- file(log_file, open="wt")
sink(log_con, type="output", split=opt$verbose)
sink(log_con, type="message")

cat("Configuration:\n")
cat("  Outcome:", metadata$options$pheno_name, "\n")
cat("  N observations:", format(nrow(DT), big.mark=","), "\n")
cat("  N events:",      format(sum(DT[[statusCol]], na.rm=TRUE), big.mark=","), "\n")
cat("  Models available:", paste(names(results$models), collapse=", "), "\n\n")

# ==============================================================================
# CALIBRATION TIME POINTS  (with --use_supported_times + guards)
# ==============================================================================
cat("CALIBRATION TIME POINTS\n")
cat("----------------------------------------------------------------------\n")

parse_times_cli <- function(x) {
  if (is.null(x) || nchar(x)==0) return(numeric(0))
  t <- suppressWarnings(as.numeric(trimws(strsplit(x, ",")[[1]])))
  sort(unique(t[is.finite(t) & t > 0]))
}

# choose candidate times
candidate_times <- numeric(0)
if (isTRUE(opt$use_supported_times)) {
  ts <- metadata$time_support
  if (!is.null(ts$eval_time_points) && length(ts$eval_time_points)>0) {
    candidate_times <- as.numeric(ts$eval_time_points)
    cat("Using eval_time_points from 01*: ", paste(candidate_times, collapse=", "), "\n", sep="")
  } else if (!is.null(ts$requested_valid_time_points) && length(ts$requested_valid_time_points)>0) {
    candidate_times <- as.numeric(ts$requested_valid_time_points)
    cat("Using requested_valid_time_points from 01*: ", paste(candidate_times, collapse=", "), "\n", sep="")
  } else if (!is.null(ts$t_pass_auto) && length(ts$t_pass_auto)>0) {
    candidate_times <- as.numeric(ts$t_pass_auto)
    cat("Using t_pass_auto from 01*: ", paste(candidate_times, collapse=", "), "\n", sep="")
  } else {
    cat("No supported times present in 01* metadata; falling back to --cal_times / auto.\n")
  }
}

# fallback to CLI or IQR-based auto if none from metadata
if (length(candidate_times) == 0) {
  if (!is.null(opt$cal_times)) {
    candidate_times <- parse_times_cli(opt$cal_times)
    cat("Using --cal_times: ", paste(candidate_times, collapse=", "), "\n", sep="")
  } else {
    max_time <- max(DT[[timeCol]], na.rm=TRUE)
    iqr      <- stats::quantile(DT[[timeCol]], c(0.25, 0.75), na.rm=TRUE)
    candidate_times <- seq(iqr[1], iqr[2], length.out=opt$n_cal_times)
    candidate_times <- round(candidate_times, 1)
    candidate_times <- sort(unique(candidate_times))
    candidate_times <- candidate_times[candidate_times > 0 & candidate_times <= max_time]
    cat("Auto (IQR) cal_times: ", paste(candidate_times, collapse=", "), "\n", sep="")
  }
}
stopifnot(length(candidate_times) > 0)

# hard cap to observed follow-up
candidate_times <- candidate_times[candidate_times <= max(DT[[timeCol]], na.rm=TRUE)]
stopifnot(length(candidate_times) > 0)

# ---------- Apply guards (G-hat, min events, min at risk) & save audit ----------
cens_surv <- survfit(Surv(DT[[timeCol]], 1 - DT[[statusCol]]) ~ 1)
audit <- data.table(time_years = candidate_times)
audit[, `:=`(
  ghat          = vapply(time_years, function(t) {
    out <- tryCatch(summary(cens_surv, times=t)$surv, error=function(e) NA_real_)
    as.numeric(out)
  }, numeric(1)),
  n_events_by_t = vapply(time_years, function(t) sum(DT[[statusCol]] == 1 & DT[[timeCol]] <= t), integer(1)),
  n_at_risk     = vapply(time_years, function(t) sum(DT[[timeCol]] >= t), integer(1))
)]
audit[, status := fifelse(!is.finite(ghat), "No estimate",
                          fifelse(ghat < opt$ghat_threshold, sprintf("FAIL (<%g)", opt$ghat_threshold),
                                  fifelse(n_at_risk < opt$min_at_risk, sprintf("FAIL (at_risk<%d)", opt$min_at_risk),
                                          fifelse(n_events_by_t < opt$min_events, sprintf("FAIL (events<%d)", opt$min_events), "PASS"))))]
cal_times <- audit[status=="PASS", time_years]
cat("Cal times after guards: ", if (length(cal_times)) paste(cal_times, collapse=", ") else "NONE", "\n\n", sep="")
if (length(cal_times)==0) stop("No calibration time points passed the guards.")

audit_file <- file.path(opt$outdir, paste0(prefix, "cal_times_ghat_audit.csv"))
fwrite(audit, audit_file)
cat("✓ Saved calibration time audit: ", basename(audit_file), "\n\n", sep="")

# ==============================================================================
# HELPERS
# ==============================================================================
# Predicted survival now comes from the shared helper (prs_risk_utils.R), which
# indexes the baseline on the model's own time scale and conditions on entry.
# The previous local version took the baseline from survfit(fit) — which for a
# start-stop model is on the START/STOP (age) scale — and indexed it with a
# follow-up-scale horizon, silently returning the wrong survival whenever
# --time_entry was used.
get_predicted_survival <- function(fit, newdata, time_point) {
  s <- predict_survival_at_horizon(fit, newdata, time_point, start_col = "START")
  if (!length(s) || all(!is.finite(s))) return(NULL)
  s
}

km_at_time <- function(df, time_col, status_col, t) {
  sf <- tryCatch(survfit(Surv(get(time_col), get(status_col)) ~ 1, data=df), error=function(e) NULL)
  if (is.null(sf)) return(list(S=NA_real_, SE=NA_real_))
  ssum <- tryCatch(summary(sf, times=t), error=function(e) NULL)
  if (is.null(ssum) || length(ssum$surv)==0) {
    if (t <= min(sf$time, na.rm=TRUE)) return(list(S=1, SE=0))
    return(list(S=utils::tail(sf$surv,1), SE=if (!is.null(sf$std.err)) utils::tail(sf$std.err,1) else NA_real_))
  }
  S  <- as.numeric(ssum$surv[1])
  SE <- if (!is.null(ssum$std.err)) as.numeric(ssum$std.err[1]) else NA_real_
  list(S=S, SE=SE)
}

grouped_calibration_from_risk <- function(DT, risk, t, n_groups=10, timeCol="TIME", statusCol="PHENO", alpha=0.05,
                                          min_cell_count=10L) {
  r <- as.numeric(risk)
  if (length(unique(r[is.finite(r)])) < 3) return(NULL)
  probs <- seq(0, 1, length.out=n_groups+1)
  brks  <- unique(stats::quantile(r, probs=probs, na.rm=TRUE))
  if (length(brks) < 3) return(NULL)
  grp   <- cut(r, breaks=brks, include.lowest=TRUE, labels=FALSE)

  tab <- data.table(group=integer(), n=integer(), pred_risk=double(),
                    obs_risk=double(), obs_se=double())
  for (g in sort(unique(grp))) {
    idx <- which(grp==g)
    # Disclosure + stability: drop calibration groups below the threshold, so no
    # small group size is published (was a hardcoded 5).
    if (length(idx) < min_cell_count) next
    df_g <- DT[idx, ]
    k    <- km_at_time(df_g, timeCol, statusCol, t)
    tab  <- rbind(tab, data.table(group=g, n=length(idx),
                                  pred_risk=mean(r[idx], na.rm=TRUE),
                                  obs_risk=1 - k$S, obs_se=k$SE))
  }
  if (nrow(tab) < 3) return(NULL)
  tab[, var_O := obs_se^2]
  tab <- tab[is.finite(var_O) & var_O > 0]
  chisq <- sum((tab$obs_risk - tab$pred_risk)^2 / tab$var_O, na.rm=TRUE)
  df    <- max(1, nrow(tab) - 1)
  pval  <- stats::pchisq(chisq, df=df, lower.tail=FALSE)
  list(table=tab, chisq=chisq, df=df, pval=pval, breaks=brks)
}

point_calibration_metrics_from_risk <- function(DT, risk, t, timeCol="TIME", statusCol="PHENO") {
  r <- as.numeric(risk)
  k_all <- km_at_time(DT, timeCol, statusCol, t)
  O_all <- 1 - k_all$S
  E_all <- mean(r, na.rm=TRUE)
  
  eps <- 1e-6
  logit <- function(x) log((pmax(pmin(x, 1-eps), eps)) / (1 - pmax(pmin(x, 1-eps), eps)))
  CIL <- logit(O_all) - logit(E_all)
  OE  <- if (E_all > 0) O_all / E_all else NA_real_
  
  q50 <- as.numeric(stats::quantile(r, 0.50, na.rm=TRUE))
  q90 <- as.numeric(stats::quantile(r, 0.90, na.rm=TRUE))
  k_nn <- max(30, floor(0.03 * sum(is.finite(r))))
  nn_idx <- function(q) order(abs(r - q))[1:k_nn]
  k50 <- km_at_time(DT[nn_idx(q50), ], timeCol, statusCol, t)
  k90 <- km_at_time(DT[nn_idx(q90), ], timeCol, statusCol, t)
  O50 <- 1 - k50$S; O90 <- 1 - k90$S
  # Signed observed-minus-predicted gap among individuals near the 50th and 90th
  # percentiles of PREDICTED risk. These are NOT "expected" risks, and NOT Austin's
  # E50/E90 (percentiles of |obs-pred| over the whole range) — named for exactly what
  # they are so the published metric is not misread as an expected value.
  obs_minus_pred_q50 <- O50 - q50;  obs_minus_pred_q90 <- O90 - q90

  list(CIL=CIL, OE=OE, O=O_all, E=E_all,
       obs_minus_pred_q50=obs_minus_pred_q50, obs_minus_pred_q90=obs_minus_pred_q90)
}

fit_recalibration_params <- function(tab) {
  out <- list()
  lin <- tryCatch(stats::lm(obs_risk ~ pred_risk, data=tab, weights=tab$n), error=function(e) NULL)
  if (!is.null(lin)) {
    co <- stats::coef(lin)
    out$lin_intercept <- unname(ifelse(is.na(co[1]), 0, co[1]))
    out$lin_slope     <- unname(ifelse(is.na(co[2]), 1, co[2]))
  } else { out$lin_intercept <- 0; out$lin_slope <- 1 }
  
  eps <- 1e-6
  safe_qlogis <- function(p) stats::qlogis(pmax(pmin(p, 1-eps), eps))
  tab$z <- safe_qlogis(tab$pred_risk)
  pl <- tryCatch(stats::glm(obs_risk ~ z, weights=n, family=binomial(), data=tab), error=function(e) NULL)
  if (!is.null(pl)) {
    cop <- stats::coef(pl)
    out$platt_a <- unname(ifelse(is.na(cop[1]), 0, cop[1]))
    out$platt_b <- unname(ifelse(is.na(cop[2]), 1, cop[2]))
  } else { out$platt_a <- 0; out$platt_b <- 1 }
  
  iso <- tryCatch(stats::isoreg(x = tab$pred_risk, y = tab$obs_risk, wt = tab$n), error=function(e) NULL)
  out$isotonic <- iso
  out
}

apply_recalibration <- function(risk, mode = c("raw","recal_int","recal_slope","recal_int_slope","platt","isotonic"),
                                params = NULL) {
  mode <- match.arg(mode)
  r <- as.numeric(risk); clamp01 <- function(x) pmin(pmax(x, 0), 1)
  eps <- 1e-6
  if (mode == "raw" || is.null(params)) return(clamp01(r))
  if (mode == "recal_int")       return(clamp01(params$lin_intercept + r))
  if (mode == "recal_slope")     return(clamp01(params$lin_slope * r))
  if (mode == "recal_int_slope") return(clamp01(params$lin_intercept + params$lin_slope * r))
  if (mode == "platt")           return(stats::plogis(params$platt_a + params$platt_b * stats::qlogis(pmax(pmin(r, 1-eps), eps))))
  if (mode == "isotonic") {
    iso <- params$isotonic
    if (is.null(iso) || is.null(iso$yf) || is.null(iso$x)) return(clamp01(r))
    ord <- order(iso$x); x0 <- iso$x[ord]; y0 <- iso$yf[ord]
    idx <- findInterval(r, x0, all.inside = TRUE)
    return(clamp01(y0[idx]))
  }
  clamp01(r)
}

bootstrap_bands_from_risk <- function(DT, risk, t, n_groups=10, B=500, seed=123,
                                      timeCol="TIME", statusCol="PHENO", level=0.95) {
  set.seed(seed); r0 <- as.numeric(risk); probs <- seq(0, 1, length.out=n_groups+1)
  keep <- vector("list", B)
  for (b in seq_len(B)) {
    idx <- sample(seq_len(nrow(DT)), replace=TRUE)
    Dtb <- DT[idx, ]; rb <- r0[idx]
    br  <- unique(stats::quantile(rb, probs=probs, na.rm=TRUE))
    if (length(br) < 3) next
    gb  <- cut(rb, breaks=br, include.lowest=TRUE, labels=FALSE)
    tab_b <- data.table(group=integer(), pred_mean=double(), obs_risk=double(), n=integer())
    for (g in sort(unique(gb))) {
      ids <- which(gb==g); if (length(ids) < 5) next
      k   <- km_at_time(Dtb[ids, ], timeCol, statusCol, t)
      tab_b <- rbind(tab_b, data.table(group=g,
                                       pred_mean=mean(rb[ids], na.rm=TRUE),
                                       obs_risk=1 - k$S,
                                       n=length(ids)))
    }
    if (nrow(tab_b) > 0) keep[[b]] <- tab_b
  }
  keep <- keep[!vapply(keep, is.null, logical(1))]
  if (length(keep) == 0) return(NULL)
  big <- data.table::rbindlist(keep, idcol="boot")
  alpha <- 1 - level
  bands <- big[, .(
    pred_mean = stats::median(pred_mean, na.rm=TRUE),
    obs_med   = stats::median(obs_risk,  na.rm=TRUE),
    obs_lo    = stats::quantile(obs_risk, probs=alpha/2, na.rm=TRUE),
    obs_hi    = stats::quantile(obs_risk, probs=1 - alpha/2, na.rm=TRUE),
    n_med     = stats::median(n, na.rm=TRUE)
  ), by=.(group)]
  bands[]
}

# ==============================================================================
# MAIN CALIBRATION
# ==============================================================================
cat("==============================================================================\n")
cat("CALIBRATION: ALL COX MODELS\n")
cat("==============================================================================\n\n")

# Collectors
all_cal_metrics  <- list()
all_group_tables <- list()
all_band_tables  <- list()
plot_points      <- list()
plot_bands       <- list()

# Choose Cox models only
models_to_assess <- names(results$models)[vapply(results$models, function(x) inherits(x, "coxph"), logical(1))]

# Restrict to the CORE M0-M7 ladder by default (parallel with stage 03). Clinical
# models are fitted on DT_clinical, a different sample, so calibrating them against
# this DT would mix samples. --calibration_models narrows further (manuscript ids
# via metadata$model_map, or raw model names).
core_models <- if (!is.null(results$metadata$model_map)) {
  mm <- as.data.table(results$metadata$model_map)
  if ("ladder" %in% names(mm)) as.character(mm[ladder == "core", model]) else as.character(mm$model)
} else NULL
if (!is.null(core_models)) {
  keep_core <- intersect(unique(c("base", core_models)), models_to_assess)
  if (length(keep_core) >= 1) models_to_assess <- keep_core
}
if (!is.null(opt$calibration_models) && nzchar(opt$calibration_models)) {
  want <- trimws(strsplit(opt$calibration_models, ",")[[1]])
  id2name <- function(id) {
    if (is.null(results$metadata$model_map)) return(id)
    mm <- as.data.table(results$metadata$model_map); hit <- mm[model_id == id, model]
    if (length(hit)) as.character(hit[1]) else id
  }
  want <- vapply(want, id2name, character(1))
  sel <- intersect(want, models_to_assess)
  if (length(sel)) models_to_assess <- sel
  else cat("⚠ --calibration_models matched none; keeping the core ladder.\n")
}
cat("Models to assess:", paste(models_to_assess, collapse=", "), "\n\n")

# Recalibration modes. DEFAULT = raw only. The recalibration modes (intercept,
# slope, Platt, isotonic) are fitted AND evaluated in the same cohort, so they are
# optimistic and must not be read as evidence that calibration improved; they are
# opt-in via --calibration_modes.
all_recal_modes <- c("raw","recal_int","recal_slope","recal_int_slope","platt","isotonic")
recal_modes <- if (!is.null(opt$calibration_modes) && nzchar(opt$calibration_modes)) {
  req <- trimws(strsplit(opt$calibration_modes, ",")[[1]])
  bad <- setdiff(req, all_recal_modes)
  if (length(bad)) stop("Unknown --calibration_modes: ", paste(bad, collapse=", "),
                        ". Valid: ", paste(all_recal_modes, collapse=", "))
  req
} else "raw"
cat("Recalibration modes:", paste(recal_modes, collapse=", "), "\n\n")

for (model_name in models_to_assess) {
  cat("Model:", model_name, "\n")
  cat("----------------------------------------------------------------------\n")
  fit <- results$models[[model_name]]
  
  for (t in cal_times) {
    cat("  Time:", t, "years\n")
    
    # RAW predictions → event risk r_i(t)
    S_i <- get_predicted_survival(fit, DT, t)
    if (is.null(S_i)) { cat("    ⚠ Could not compute predictions; skipping\n"); next }
    r_raw <- 1 - S_i

    # RESTRICT TO FINITE PREDICTIONS, AND USE THE SAME ROWS FOR OBSERVED RISK.
    # predict_risk_at_horizon() now returns NA where the horizon lies outside the
    # baseline support. Left unfiltered, the predicted mean would be taken over
    # the finite rows while the observed KM used everyone — so O/E,
    # calibration-in-the-large and the grouped table would compare two different
    # samples. DT_t is the sample every metric at this horizon is computed on.
    finite_pred <- is.finite(r_raw)
    n_unsupported <- sum(!finite_pred)
    if (n_unsupported > 0)
      cat(sprintf("    %d/%d individual(s) have NA predicted risk at t=%g (outside baseline support); excluded from BOTH predicted and observed.\n",
                  n_unsupported, length(r_raw), t))
    if (!any(finite_pred)) { cat("    ⚠ No finite predictions at this horizon; skipping\n"); next }
    DT_t  <- DT[finite_pred]
    r_raw <- r_raw[finite_pred]

    if (stats::sd(r_raw, na.rm=TRUE) < 1e-8) { cat("    ⚠ Near-constant predicted risk; skipping\n"); next }

    # Base grouped calibration for raw predictions (to fit recalibration params)
    gc_raw <- grouped_calibration_from_risk(DT_t, risk = r_raw, t = t,
                                            n_groups = opt$n_groups, timeCol=timeCol, statusCol=statusCol,
                                            min_cell_count = opt$min_cell_count)
    if (is.null(gc_raw)) { cat("    ⚠ Grouped calibration (raw) not available; skipping time\n"); next }
    
    # Fit recalibration params (from raw grouped table)
    params <- fit_recalibration_params(gc_raw$table)
    
    # Iterate all modes
    for (mode_label in recal_modes) {
      # Recalibrated predictions for this mode
      r_mode <- apply_recalibration(r_raw, mode = mode_label, params = params)
      
      # Recompute grouped calibration using recalibrated risk (new quantile groups)
      gc <- grouped_calibration_from_risk(DT_t, risk = r_mode, t = t,
                                          n_groups = opt$n_groups, timeCol=timeCol, statusCol=statusCol,
                                          min_cell_count = opt$min_cell_count)
      if (is.null(gc)) { cat("    ⚠ Grouped calibration not available for mode:", mode_label, "\n"); next }
      
      # Metrics on grouped table
      tbl <- data.table(gc$table)  # has: group, n, pred_risk, obs_risk, obs_se, var_O
      ICI_val <- mean(abs(tbl$pred_risk - tbl$obs_risk), na.rm=TRUE)
      E_stat  <- sum((tbl$pred_risk - tbl$obs_risk)^2, na.rm=TRUE)
      slope <- intercept <- r2 <- NA_real_
      lmfit <- tryCatch(stats::lm(obs_risk ~ pred_risk, data = tbl, weights = tbl$n), error=function(e) NULL)
      if (!is.null(lmfit)) {
        co  <- stats::coef(lmfit)
        slope <- unname(ifelse(is.na(co[2]), NA, co[2]))
        intercept <- unname(ifelse(is.na(co[1]), NA, co[1]))
        r2 <- summary(lmfit)$r.squared
      }
      
      # Point metrics (CIL, O/E, obs_minus_pred_q50/q90) for this mode
      pm <- point_calibration_metrics_from_risk(DT_t, risk = r_mode, t = t, timeCol=timeCol, statusCol=statusCol)
      
      cat(sprintf("    [%s] χ²=%.3f (df=%d), p=%.3g | slope=%.3f, int=%.3f | ICI=%.4f | CIL=%.3f, O/E=%.3f\n",
                  mode_label, gc$chisq, gc$df, gc$pval, slope, intercept, ICI_val, pm$CIL, pm$OE))
      
      # Save grouped table
      tbl_out <- data.table(model = model_name, time = t, mode = mode_label, tbl)
      all_group_tables[[paste(model_name, t, mode_label, sep="__")]] <- tbl_out
      
      # Save metrics
      all_cal_metrics[[paste(model_name, t, mode_label, sep="__")]] <- data.table(
        model = model_name, time = t, mode = mode_label,
        # NOT the standard individual-level calibration slope / calibration-in-the-
        # large. These come from a weighted regression of GROUPED observed risk on
        # GROUPED predicted risk (obs_risk ~ pred_risk over the risk deciles), so
        # they are named for what they are and must not be reported as the usual
        # survival calibration slope/intercept.
        grouped_calibration_slope = slope,
        grouped_calibration_intercept = intercept,
        r_squared = r2,
        ICI = ICI_val, E_statistic = E_stat,
        CIL = pm$CIL, OE = pm$OE, O_overall = pm$O, E_overall = pm$E,
        obs_minus_pred_q50 = pm$obs_minus_pred_q50, obs_minus_pred_q90 = pm$obs_minus_pred_q90,
        grouped_chisq = gc$chisq, grouped_df = gc$df, grouped_p = gc$pval,
        n_groups = nrow(tbl),
        # Sample actually used at this horizon: predicted and observed are computed
        # on the SAME rows (see the finite-prediction filter above).
        n_analysed = nrow(DT_t), n_dropped_unsupported = n_unsupported
      )
      
      # Optional bootstrap ribbons
      bands <- NULL
      if (isTRUE(opt$bootstrap_bands)) {
        bands <- bootstrap_bands_from_risk(DT_t, risk = r_mode, t = t,
                                           n_groups = opt$n_groups,
                                           B = opt$n_boot, seed = opt$boot_seed,
                                           timeCol = timeCol, statusCol = statusCol,
                                           level = opt$ci_level)
        if (!is.null(bands)) {
          bands[, `:=`(model = model_name, time = t, mode = mode_label)]
          data.table::setcolorder(bands, c("model","time","mode","group","pred_mean","obs_lo","obs_med","obs_hi","n_med"))
          all_band_tables[[paste(model_name, t, mode_label, sep="__")]] <- bands
        }
      }
      
      # Plot collectors
      pp <- data.table(tbl_out)[, .(model, time, mode, group, n, pred_risk, obs_risk)]
      plot_points[[length(plot_points) + 1]] <- pp
      if (!is.null(bands)) plot_bands[[length(plot_bands) + 1]] <- data.table(bands)
      
      # Plots (now opt-in via --make_plots)
      if (isTRUE(opt$make_plots)) {
        p <- ggplot(tbl_out, aes(x=pred_risk, y=obs_risk)) +
          geom_abline(slope=1, intercept=0, linetype="dashed", color="gray50") +
          geom_point(aes(size=n), alpha=0.7, color="#2c7fb8") +
          { if (opt$smooth_calibration && nrow(tbl_out) >= 5)
            geom_smooth(method="loess", se=TRUE, color="#d95f02", fill="#d95f02", alpha=0.15) } +
          labs(
            title = sprintf("Calibration at t=%.1f years (%s | %s)", t, model_name, mode_label),
            x = "Predicted risk by t",
            y = "Observed risk by t (KM)",
            size = "N"
          ) +
          coord_equal(xlim=c(0,1), ylim=c(0,1)) +
          theme_minimal(base_size=12) +
          theme(legend.position="bottom", plot.title=element_text(face="bold", hjust=0.5))
        
        if (isTRUE(opt$bootstrap_bands) && !is.null(bands)) {
          p <- p + geom_ribbon(data=bands, aes(x=pred_mean, ymin=obs_lo, ymax=obs_hi),
                               inherit.aes=FALSE, alpha=0.15, fill="#2c7fb8")
        }
        
        formats <- trimws(strsplit(opt$format, ",")[[1]])
        save_plot_cairo(p, file.path(opt$outdir, paste0(prefix, "cal_", model_name, "_t", t, "_", mode_label, ".png")),
                        width=opt$plot_width, height=opt$plot_height, dpi=opt$dpi, formats=formats)
        cat(sprintf("    ✓ Saved plot(s) for mode: %s\n", mode_label))
      }
    } # modes
  }   # times
  cat("\n")
}

# ==============================================================================
# SAVE TABLES
# ==============================================================================
cal_groups_file  <- file.path(opt$outdir, paste0(prefix, "calibration_groups.csv"))
cal_metrics_file <- file.path(opt$outdir, paste0(prefix, "calibration_metrics.csv"))
cal_bands_file   <- file.path(opt$outdir, paste0(prefix, "calibration_bands.csv"))
by_time_file     <- file.path(opt$outdir, paste0(prefix, "calibration_by_time.csv"))

if (length(all_group_tables) > 0) {
  groups_dt <- data.table::rbindlist(all_group_tables, fill=TRUE)
  data.table::setorder(groups_dt, model, time, mode, group)
  data.table::fwrite(groups_dt, cal_groups_file)
  cat("✓ Saved grouped calibration points:", basename(cal_groups_file), "\n")
} else {
  cat("ℹ No grouped calibration points to save\n")
}

if (length(all_cal_metrics) > 0) {
  metrics_dt <- data.table::rbindlist(all_cal_metrics, fill=TRUE)
  data.table::setorder(metrics_dt, model, time, mode)
  data.table::fwrite(metrics_dt, cal_metrics_file)
  cat("✓ Saved calibration metrics:", basename(cal_metrics_file), "\n")
} else {
  cat("ℹ No calibration metrics to save\n")
}

if (length(all_band_tables) > 0) {
  bands_dt <- data.table::rbindlist(all_band_tables, fill=TRUE)
  data.table::setorder(bands_dt, model, time, mode, group)
  data.table::fwrite(bands_dt, cal_bands_file)
  cat("✓ Saved bootstrap calibration bands:", basename(cal_bands_file), "\n")
} else if (isTRUE(opt$bootstrap_bands)) {
  cat("ℹ Bands requested but none saved (insufficient data during bootstraps)\n")
}

if (length(all_group_tables) > 0) {
  by_time <- data.table::rbindlist(all_group_tables, fill=TRUE)[,
                                                                .(n_groups = .N,
                                                                  mean_abs_error = mean(abs(pred_risk - obs_risk), na.rm=TRUE),
                                                                  mean_sq_error  = mean((pred_risk - obs_risk)^2, na.rm=TRUE)),
                                                                by=.(model, time, mode)]
  data.table::setorder(by_time, model, time, mode)
  data.table::fwrite(by_time, by_time_file)
  cat("✓ Saved calibration-by-time summary:", basename(by_time_file), "\n")
}

# -------- Save combined tidy plot data --------
plot_points_file <- file.path(opt$outdir, paste0(prefix, "calibration_plot_points.csv"))
plot_bands_file  <- file.path(opt$outdir, paste0(prefix, "calibration_plot_bands.csv"))

if (length(plot_points)) {
  pp <- data.table::rbindlist(plot_points, fill = TRUE)
  for (nm in c("time","group","n","pred_risk","obs_risk")) if (nm %in% names(pp)) pp[[nm]] <- as.numeric(pp[[nm]])
  data.table::setorder(pp, model, time, mode, group)
  data.table::fwrite(pp, plot_points_file)
  cat("✓ Saved plot points:", basename(plot_points_file), "\n")
} else {
  cat("ℹ No plot points to save (none computed)\n")
}

if (length(plot_bands)) {
  pb <- data.table::rbindlist(plot_bands, fill = TRUE)
  for (nm in c("time","group","pred_mean","obs_lo","obs_med","obs_hi","n_med")) if (nm %in% names(pb)) pb[[nm]] <- as.numeric(pb[[nm]])
  data.table::setorder(pb, model, time, mode, group)
  data.table::fwrite(pb, plot_bands_file)
  cat("✓ Saved plot ribbons:", basename(plot_bands_file), "\n")
} else if (isTRUE(opt$bootstrap_bands)) {
  cat("ℹ Bands requested but none saved (insufficient data)\n")
}

# ==============================================================================
# OPTIONAL TRAJECTORY PLOTS
# ==============================================================================
if (isTRUE(opt$make_plots) && length(all_cal_metrics) > 0) {
  md <- data.table::rbindlist(all_cal_metrics, fill=TRUE)
  
  p_slope <- ggplot(md, aes(x=time, y=grouped_calibration_slope, color=mode, group=interaction(model, mode))) +
    geom_hline(yintercept=1, linetype="dashed", color="gray50") +
    geom_line(linewidth=1) + geom_point(size=2) +
    labs(title="Grouped Calibration Slope Over Time", x="Time (years)",
         y="Grouped calibration slope", color="Mode") +
    theme_minimal(base_size=12) + theme(legend.position="bottom")
    save_plot_cairo(p_slope, file.path(opt$outdir, paste0(prefix, "slope_trajectory.", "png")), # Use .png as a base filename, save_plot_cairo handles formats
                    width=8, height=6, dpi=opt$dpi, formats=trimws(strsplit(opt$format, ",")[[1]]))
  
  p_ici <- ggplot(md, aes(x=time, y=ICI, color=mode, group=interaction(model, mode))) +
    geom_hline(yintercept=0, linetype="dashed", color="gray50") +
    geom_line(linewidth=1) + geom_point(size=2) +
    labs(title="Integrated Calibration Index Over Time", x="Time (years)", y="ICI", color="Mode") +
    theme_minimal(base_size=12) + theme(legend.position="bottom")
    save_plot_cairo(p_ici, file.path(opt$outdir, paste0(prefix, "ici_trajectory.", "png")), # Use .png as a base filename, save_plot_cairo handles formats
                    width=8, height=6, dpi=opt$dpi, formats=trimws(strsplit(opt$format, ",")[[1]]))
  
  cat("✓ Saved trajectory plots\n")
}

# ==============================================================================
# SUMMARY
# ==============================================================================
cat("==============================================================================\n")
cat("SUMMARY\n")
cat("==============================================================================\n\n")
cat("✓ Calibration analysis completed (modes: ", paste(recal_modes, collapse=", "), ")\n\n", sep="")
if (length(all_cal_metrics) > 0) {
  md <- data.table::rbindlist(all_cal_metrics)
  cat("Models assessed:", paste(unique(md$model), collapse=", "), "\n")
  cat("Times assessed:",  paste(sort(unique(md$time)), collapse=", "), "years\n")
  cat("Modes applied:",   paste(unique(md$mode), collapse=", "), "\n\n")
}
cat("Outputs saved to:", opt$outdir, "\n")
cat("Log:", log_file, "\n")

sink(type="message"); sink(type="output"); close(log_con)
