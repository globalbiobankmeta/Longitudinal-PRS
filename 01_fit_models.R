#!/usr/bin/env Rscript
# ==============================================================================
# 01_fit_models.R - Centralized Model Fitting & Intermediate Output Storage
# ==============================================================================
# Purpose:
#   - Load and process PRS, phenotype, and covariate data
#   - Fit Cox models ONCE and store all intermediate outputs
#   - Provide reusable objects for downstream evaluation metric scripts
#   - Avoid redundant model fitting across multiple analyses
#
# Evaluation metrics provided:
#   C-index, HRs, LRT, PH tests, Royston & Sauerbrei's R²
#
# NEW (this version):
#   - Default survival is start–stop: Surv(START, STOP, PHENO)
#   - START defaults to 0 unless --time_entry is provided
#   - STOP is provided by --time_exit (default: time_to_event)
#   - Reverse-KM and follow-up QC always use FOLLOWUP = STOP - START
#   - Back-compat aliases: --time_col -> --time_exit; --index_time_col/--entry_time_col -> --time_entry
#
# Outputs:
# ├── prefix_fitted_models.rds
# ├── prefix_data_processed.rds
# ├── prefix_model_metadata.rds
# ├── prefix_model_comparison.csv
# ├── prefix_all_coefficients.csv
# ├── prefix_ph_tests.csv
# ├── prefix_time_summary.csv
# ├── prefix_followup_reverseKM.csv
# ├── prefix_time_support.csv
# ├── prefix_requested_timepoints_audit.csv
# ├── prefix_requested_timepoints_ghat.csv
# ├── prefix_harrells_ci_details.csv
# ├── prefix_uno_cindex_details.csv           # if --calculate_uno
# ├── prefix_fit_models.log
# └── prefix_load_models.R
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
required_packages <- c("optparse", "data.table", "survival", "splines")
invisible(lapply(required_packages, load_package))

# ==============================================================================
# COMMAND-LINE ARGUMENTS
# ==============================================================================
option_list <- list(
  # Input files
  make_option("--prs_file", type="character",
              help="Path to PRS score file (.profile or .sscore) [REQUIRED]"),
  make_option("--pheno_file", type="character",
              help="Path to phenotype file (with time and status) [REQUIRED]"),
  make_option("--cov_file", type="character", default=NULL,
              help="Path to covariate file [default: NULL]"),
  make_option("--pop_file", type="character", default=NULL,
              help="Population subset file (FID IID format) [default: NULL]"),
  
  # Column names
  make_option("--id_col", type="character", default="IID",
              help="ID column name [default: IID]"),
  make_option("--prs_col", type="character", default="SCORESUM",
              help="PRS column name [default: SCORESUM]"),
  
  # >>> Age-derived time scale (preferred). The script derives the T1->T2 time
  # from two AGE columns, so biobanks need not precompute a time column:
  #     time_since_T1 = age_at_exit - age_at_T1  ;  START = 0 ; STOP = time_since_T1
  # age_at_T1 also serves as the ns() age covariate (one input, two uses).
  make_option("--age_t1_col",   type="character", default="diagAge",
              help="Age at T1 diagnosis (index age AND spline covariate) [default: diagAge]"),
  make_option("--age_exit_col", type="character", default=NULL,
              help="Age at T2 (events) or last follow-up (non-events). If set with --age_t1_col, STOP is derived as age_exit - age_t1."),
  make_option("--age_spline_df", type="integer", default=3L,
              help="df for ns(age_at_T1) in every model; 0 = enter age linearly [default: 3]"),

  # >>> start–stop columns (used directly if age columns are not supplied)
  make_option("--time_exit",  type="character", default="time_to_event",
              help="EXIT/STOP time column (used when --age_exit_col is not given) [default: time_to_event]"),
  make_option("--time_entry", type="character", default=NULL,
              help="ENTRY/START time column (e.g., age at entry/index). If omitted, START=0 [default: NULL]"),

  # Back-compat (aliases)
  make_option("--time_col",       type="character", default=NULL,
              help="DEPRECATED alias for --time_exit"),
  make_option("--index_time_col", type="character", default=NULL,
              help="DEPRECATED alias for --time_entry"),
  make_option("--entry_time_col", type="character", default=NULL,
              help="DEPRECATED alias for --time_entry"),
  
  make_option("--allow_missing_covariates", action="store_true", default=FALSE,
              help="Downgrade a missing requested covariate from an error to a warning [default: FALSE = strict]"),

  make_option("--status_col", type="character", default="PHENO",
              help="Event status column name [default: PHENO]"),
  make_option("--covariates", type="character", default="",
              help="Comma-separated basic covariate names (e.g., 'age,sex,PC1,PC2,...,PC10')"),
  make_option("--clinical_covariates", type="character", default=NULL,
              help="Comma-separated clinical risk factor names [default: NULL]"),
  make_option("--sex_col", type="character", default=NULL,
              help="Sex/gender column name for summaries [default: NULL]"),
  make_option("--sex_labels", type="character", default="0=Female,1=Male",
              help="Sex coding: value=label pairs [default: 0=Female,1=Male]"),
  
  # PRS processing
  make_option("--prs_signs", type="character", default=NULL,
              help="Prespecified per-PRS direction multipliers +1/-1, comma-separated in --prs_labels (onset,progression,outcome) order. -1 = flip that score, +1 = keep. E.g. '-1,-1,-1' flips all three; '1,-1,1' flips only progression. Use for a score with a known reversed allele/sign convention."),
  make_option("--auto_flip_prs", action="store_true", default=FALSE,
              help="DISCOURAGED: auto-flip a PRS if it is negatively correlated with the OUTCOME in the target cohort. This is target-outcome dependent; direction should come from discovery-stage allele harmonisation. [default: FALSE]"),
  make_option("--zscore_by_controls", action="store_true", default=FALSE,
              help="DISCOURAGED: z-score PRS using control (non-event) mean/SD. Affects only the HR-per-SD scale (not LRT/AUC) but is weakly outcome-dependent; whole-sample scaling is the default [default: FALSE]"),
  
  # Analysis options & guards
  make_option("--prs_quantiles", type="character", default="5",
              help="Comma-separated quantile specs [default: 5]"),
  make_option("--prs_quantile_reference", type="character", default="lowest",
              help="Reference group: 'lowest', 'middle', or e.g. 'Q3' [default: lowest]"),
  make_option("--prs_extremes", type="character", default=NULL,
              help="Extreme percentiles (e.g., '1,5') [default: NULL]"),
  make_option("--age_col", type="character", default=NULL,
              help="Age column name (for metadata only) [default: NULL]"),
  make_option("--min_events",     type="integer", default=10,
              help="Minimum cumulative events by t [default: 10]"),
  make_option("--min_at_risk",    type="integer", default=50,
              help="Minimum at-risk at t [default: 50]"),
  make_option("--ghat_threshold", type="numeric", default=0.02,
              help="Minimum G-hat(t) [default: 0.02]"),
  make_option("--min_cell_count", type="integer", default=10,
              help="Disclosure: blank subgroup counts (timing/T1-status) below this [default: 10]"),
  make_option("--min_events_total", type="integer", default=50,
              help="Hard floor on total T2 events; stop below this (full model not estimable) [default: 50]"),

  # >>> Analysis mode & prospective (Layer 2)
  make_option("--analysis_mode", type="character", default="gwas_aligned",
              help="'gwas_aligned' (Layer 1, START=0, STOP=time since T1) or 'prospective' (Layer 2, index at recruitment/T1) [default: gwas_aligned]"),
  make_option("--recruit_file", type="character", default=NULL,
              help="Prospective mode: file with IID + age at recruitment [default: NULL]"),
  make_option("--recruit_age_col", type="character", default="age_at_recruitment",
              help="Column in --recruit_file holding age at recruitment [default: age_at_recruitment]"),
  make_option("--incident_lag_days", type="numeric", default=0,
              help="Prospective mode: add this lag (days) to the T1 index for incident T1 [default: 0]"),
  # Which T1 states enter the prospective analysis. 'both' (default, unchanged
  # behaviour) keeps prevalent T1 alongside incident T1, indexing prevalent people at
  # recruitment and separating the two with strata(T1_STATUS) + a T1_DURATION term.
  # 'incident' restricts to T1 occurring after recruitment, giving one clean
  # time-since-T1 clock at the cost of the prevalent sample.
  make_option("--prospective_t1", type="character", default="both",
              help="Prospective mode: which T1 states to include, 'both' or 'incident' [default: both]"),
  # >>> Lag sensitivity (Layer 3): drop anyone with time_since_T1 <= lag/365.25
  make_option("--lag_days", type="numeric", default=0,
              help="Exclude transitions with time_since_T1 <= lag_days (a single value; the driver sweeps 0/30/90/365) [default: 0]"),
  
  # Requested evaluation horizons (years)
  make_option("--time_points", type="character", default="",
              help="Comma-separated time points (years) to audit & save (e.g., '1,2,5,10')."),
  
  # Uno's C-index (optional)
  make_option("--calculate_uno", action="store_true", default=FALSE,
              help="Calculate Uno's C-index for all models [default: FALSE]"),
  make_option("--uno_bootstrap", action="store_true", default=FALSE,
              help="Bootstrap CI for Uno's C (slower) [default: FALSE]"),
  make_option("--uno_bootstrap_n", type="integer", default=500,
              help="Bootstrap iterations for Uno's CI [default: 500]"),
  make_option("--tau",            type="numeric", default=NA,
              help="Integration horizon for Uno's C (years). NA -> auto via guards."),
  
  # Stratification
  make_option("--strata", type="character", default=NULL,
              help="Stratification variable (must exist in data) [default: NULL]"),
  make_option("--strata_var", type="character", default=NULL,
              help="Variable to create strata from [default: NULL]"),
  make_option("--strata_breaks", type="character", default=NULL,
              help="Breakpoints for strata [default: NULL]"),
  make_option("--strata_labels", type="character", default=NULL,
              help="Labels for strata [default: auto]"),
  make_option("--test_strata_interaction", action="store_true", default=FALSE,
              help="Test PRS × strata interaction [default: FALSE]"),
  
  # Output
  make_option("--outdir", type="character", default="output",
              help="Output directory [default: output]"),
  make_option("--cohort", type="character", default=NULL,
              help="Cohort for prefix [default: NULL]"),
  make_option("--ancestry", type="character", default=NULL,
              help="Ancestry for prefix [default: NULL]"),
  make_option("--pheno_name", type="character", default=NULL,
              help="Phenotype name for prefix [default: NULL]"),
  make_option("--prs_method", type="character", default=NULL,
              help="PRS method, inserted in the auto-prefix between ancestry and trait [default: NULL]"),
  make_option("--prefix", type="character", default="",
              help="Manual prefix [default: auto-generated]"),
  make_option("--verbose", action="store_true", default=TRUE,
              help="Print detailed output [default: TRUE]"),
  
  # Time summaries (descriptive only)
  make_option("--time_summary_cols", type="character", default="",
              help="Comma-separated column names to add to descriptive time summaries"),
  
  # Follow-up threshold. Default 0 = keep every valid ordered transition (STOP>START),
  # matching the progression GWAS, which used no washout. Same-age T1/T2 (STOP=0) are
  # still excluded by STOP>START. Raise this only for a deliberate washout analysis.
  make_option("--min_followup_threshold", type="numeric", default=0,
              help="Minimum follow-up time to keep (STOP − START) [default: 0 = no washout]"),

  # ===========================================================================
  # MULTI-PRS OPTIONS (optional; backward-compatible)
  # ===========================================================================
  # Supply additional PRS files beyond the primary --prs_file.
  # When --prs_labels is provided the script enters multi-PRS mode and builds
  # every model combination defined by --models (or auto-built from --prs_types).
  # The primary --prs_file / --prs_col is treated as the FIRST entry; its label
  # should be the first element of --prs_labels.
  make_option("--prs_files", type="character", default=NULL,
              help="Comma-separated ADDITIONAL PRS file paths (beyond --prs_file). Order must match --prs_labels[2..N]."),
  make_option("--prs_labels", type="character", default=NULL,
              help="Comma-separated labels for ALL PRS (primary first). E.g. 'Onset,Progression,Outcome'. Activates multi-PRS mode."),
  make_option("--prs_types", type="character", default=NULL,
              help="Comma-separated PRS roles, matching --prs_labels order. Each must be: onset | progression | outcome."),
  make_option("--prs_col_list", type="character", default=NULL,
              help="Comma-separated score column names, one per PRS file. Falls back to --prs_col for unlisted entries."),
  make_option("--models", type="character", default=NULL,
              help="Semicolon-separated explicit model specs: 'ModelName=Label1+Label2'. If NULL and multi-PRS mode is active, auto-builds all combinations from --prs_types.")
  # NOTE: --include_base was removed. It was declared but never referenced, and
  # the covariate-only Base model (M0) is fitted unconditionally, so the flag
  # could never change behaviour.
)

opt <- parse_args(OptionParser(
  option_list=option_list,
  description="\nFit Cox models with default start–stop Surv(START, STOP, PHENO). START=0 unless --time_entry is provided; reverse-KM uses FOLLOWUP=STOP−START.",
  epilogue="
Examples:
  Rscript 01_fit_models.R --prs_file PRS.sscore --pheno_file pheno.csv --time_exit age_event --status_col PHENO --time_points '1,2,5,10'
  Rscript 01_fit_models.R --prs_file PRS.sscore --pheno_file pheno.csv --time_entry age_enroll --time_exit age_event
"
))

# Back-compat: promote deprecated aliases
if (!is.null(opt$time_col) && nzchar(opt$time_col) && (is.null(opt$time_exit) || !nzchar(opt$time_exit))) {
  warning("--time_col is deprecated; use --time_exit. Promoting value to --time_exit.")
  opt$time_exit <- opt$time_col
}
if (!is.null(opt$index_time_col) && nzchar(opt$index_time_col) && (is.null(opt$time_entry) || !nzchar(opt$time_entry))) {
  warning("--index_time_col is deprecated; use --time_entry. Promoting value to --time_entry.")
  opt$time_entry <- opt$index_time_col
}
if (!is.null(opt$entry_time_col) && nzchar(opt$entry_time_col) && (is.null(opt$time_entry) || !nzchar(opt$time_entry))) {
  warning("--entry_time_col is deprecated; use --time_entry. Promoting value to --time_entry.")
  opt$time_entry <- opt$entry_time_col
}

# Required
if (is.null(opt$prs_file)) stop("--prs_file is required")
if (is.null(opt$pheno_file)) stop("--pheno_file is required")
# Validate analysis_mode explicitly (an invalid value must NOT silently behave
# like gwas_aligned).
opt$analysis_mode <- tryCatch(match.arg(opt$analysis_mode, c("gwas_aligned", "prospective")),
                              error = function(e) stop("--analysis_mode must be 'gwas_aligned' or 'prospective' (got '", opt$analysis_mode, "')"))

# PRS direction options must not conflict (a score could otherwise be flipped twice).
if (!is.null(opt$prs_signs) && nzchar(opt$prs_signs) && isTRUE(opt$auto_flip_prs))
  stop("--prs_signs cannot be combined with --auto_flip_prs (double-flip risk). ",
       "Use exactly one direction mechanism.")

# ==============================================================================
# OUTPUT PREFIX & LOG
# ==============================================================================
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)
if (nchar(opt$prefix) > 0) {
  prefix <- paste0(opt$prefix, "_")
} else {
  prefix_parts <- c()
  if (!is.null(opt$cohort)   && opt$cohort   != "") prefix_parts <- c(prefix_parts, opt$cohort)
  if (!is.null(opt$ancestry) && opt$ancestry != "") prefix_parts <- c(prefix_parts, opt$ancestry)
  # PRS method between ancestry and trait, matching the driver stem
  # ${COHORT}_${ANCESTRY}_${PRS_METHOD}_${TRAIT}, so FM_MAIN resolves and a P+T run
  # cannot overwrite a PRS-CS run.
  if (!is.null(opt$prs_method) && opt$prs_method != "") prefix_parts <- c(prefix_parts, opt$prs_method)
  if (!is.null(opt$pheno_name) && opt$pheno_name != "") prefix_parts <- c(prefix_parts, opt$pheno_name)
  if (!is.null(opt$strata_var) && opt$strata_var != "") prefix_parts <- c(prefix_parts, paste0("strata_", opt$strata_var))
  else if (!is.null(opt$strata) && opt$strata != "")  prefix_parts <- c(prefix_parts, paste0("strata_", opt$strata))
  # Include analysis mode / lag in the auto-prefix so gwas_aligned, prospective and
  # lagged runs into the same directory do not overwrite each other. NON-DEFAULT
  # only, so the default gwas_aligned + lag-0 prefix is unchanged (the driver's
  # FM_MAIN path still resolves).
  if (!identical(opt$analysis_mode, "gwas_aligned")) prefix_parts <- c(prefix_parts, opt$analysis_mode)
  if (is.finite(opt$lag_days) && opt$lag_days > 0) prefix_parts <- c(prefix_parts, paste0("lag", opt$lag_days, "d"))
  if (is.finite(opt$incident_lag_days) && opt$incident_lag_days > 0) prefix_parts <- c(prefix_parts, paste0("incidentLag", opt$incident_lag_days, "d"))
  prefix <- if (length(prefix_parts) > 0) paste0(paste(prefix_parts, collapse="_"), "_01_fit_models_") else "01_fit_models_"
}
cat("Output prefix:", prefix, "\n\n")

log_file <- file.path(opt$outdir, paste0(prefix, "fit_models.log"))
log_con <- file(log_file, open="wt")
sink(log_con, type="output", split=opt$verbose)
sink(log_con, type="message")

cat("==============================================================================\n")
cat("CENTRALIZED MODEL FITTING - PRS EVALUATION\n")
cat("==============================================================================\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("==============================================================================\n\n")

# ==============================================================================
# HELPERS
# ==============================================================================
# Duplicate IIDs cause many-to-many blow-up in merges. Error out early with a
# clear message (a harmonised analysis should have one row per individual).
assert_unique_iid <- function(dt, id = "IID", what = "file") {
  if (!id %in% names(dt)) return(invisible(NULL))
  dups <- sum(duplicated(dt[[id]]))
  if (dups > 0) stop(what, " has ", dups, " duplicated ", id,
                     " value(s); de-duplicate before running (one row per individual).")
  invisible(NULL)
}

# Cohort-flow (attrition) tracking: append (stage, N, events) as the pipeline
# filters. status_col name varies (PHENO after recode, else the raw status col).
attrition_log <- list()
record_stage <- function(dt, label, status_col = "PHENO") {
  n  <- if (is.null(dt)) 0L else nrow(dt)
  ev <- if (!is.null(dt) && status_col %in% names(dt)) sum(dt[[status_col]] == 1, na.rm = TRUE) else NA_integer_
  attrition_log[[length(attrition_log) + 1L]] <<- data.table::data.table(stage = label, N_remaining = n, events_remaining = ev)
  invisible(NULL)
}

# ------------------------------------------------------------------------------
# .aeqsurv_keep() — which rows survive survival's own time-resolution rounding.
#
# survival::aeqSurv() (called by coxph via timefix=TRUE) rounds event times to the
# data's effective resolution and HARD-ERRORS — "aeqSurv exception, an interval has
# effective length 0" — if any interval collapses. A strict STOP > START test does
# NOT prevent this: two ages derived from the SAME date can differ by ~1e-13 through
# float cancellation, pass `> START`, and then round to zero length. Those rows are
# same-age T1/T2 transitions, which this pipeline already excludes by policy (see
# --min_followup_threshold), so apply aeqSurv's own rule and make the exclusion
# float-safe instead of leaving coxph to abort mid-fit.
#
# The rule is aeqSurv's verbatim: consecutive sorted unique times merge when the gap
# is <= tolerance OR gap/mean(abs(times)) <= tolerance. Returns a logical vector,
# TRUE = the interval still spans at least one cut boundary.
# ------------------------------------------------------------------------------
.aeqsurv_keep <- function(start, stop, tolerance = sqrt(.Machine$double.eps)) {
  y <- sort(unique(c(start, stop))); y <- y[is.finite(y)]
  if (length(y) < 2L) return(rep(TRUE, length(start)))
  dy   <- diff(y)
  tied <- (dy <= tolerance) | (dy / mean(abs(y)) <= tolerance)
  if (!any(tied)) return(rep(TRUE, length(start)))
  cuts <- y[c(TRUE, !tied)]
  findInterval(start, cuts) != findInterval(stop, cuts)
}

# Drop rows whose Surv interval collapses under aeqSurv, asking aeqSurv itself
# (the arbiter coxph uses) rather than trusting one pass of the rule above: the
# cut points and mean(abs(y)) both shift as rows are removed. Returns the subset.
.drop_degenerate_intervals <- function(dt, what = "sample") {
  for (.i in 1:5) {
    ok <- tryCatch({ survival::aeqSurv(survival::Surv(dt$START, dt$STOP, dt$PHENO)); TRUE },
                   error = function(e) FALSE)
    if (ok) return(dt)
    k <- .aeqsurv_keep(dt$START, dt$STOP)
    if (all(k))
      stop("Surv intervals in the ", what, " still collapse under aeqSurv, but no row ",
           "is identifiable as degenerate; inspect START/STOP.")
    cat("  ⚠ Removing ", sum(!k), " row(s) from the ", what,
        " whose interval collapses at survival's time resolution\n", sep = "")
    dt <- dt[k]
  }
  stop("Could not resolve degenerate Surv intervals in the ", what, " after 5 passes.")
}

# Disclosure helper (identical to prs_risk_utils::suppress_small_cells; inlined
# because this script sources nothing). Suppresses a row when ANY count column is a
# small non-zero integer (1..threshold-1), NA-ing counts + estimates.
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

# Partition table where the rows sum to a total that IS published elsewhere (e.g.
# t1_status_final / timing_categories vs the overall N in final_analysis_summary).
# Per-row suppression is insufficient — a suppressed row is reconstructable by
# subtraction — so suppress ALL rows when ANY row is small.
suppress_partition_all_or_none <- function(x, count_cols, threshold = 10L) {
  count_cols <- intersect(count_cols, names(x))
  if (!length(count_cols) || !nrow(x)) { x[, minimum_cell_count_pass := TRUE]; return(x[]) }
  any_small <- any(unlist(lapply(count_cols, function(cc) {
    v <- suppressWarnings(as.numeric(x[[cc]])); any(is.finite(v) & v > 0 & v < threshold)
  })))
  x[, minimum_cell_count_pass := !any_small]
  if (any_small) for (cc in count_cols) x[, (cc) := NA]
  x[]
}

parse_sex_labels <- function(sex_labels_str) {
  if (is.null(sex_labels_str) || sex_labels_str == "") return(NULL)
  pairs <- strsplit(sex_labels_str, ",")[[1]]
  sex_map <- list()
  for (pair in pairs) {
    parts <- strsplit(trimws(pair), "=")[[1]]
    if (length(parts) == 2) sex_map[[trimws(parts[1])]] <- trimws(parts[2])
  }
  sex_map
}

compute_supported_times_from_data <- function(DT, time_col="TIME", status_col="PHENO",
                                              min_events=10L, min_at_risk=50L,
                                              ghat_threshold=0.02, cap_quantile=0.95) {
  time   <- as.numeric(DT[[time_col]])
  status <- as.integer(DT[[status_col]])
  kmC <- survival::survfit(survival::Surv(time, 1L - status) ~ 1)
  getG <- function(t) {
    idx <- max(which(kmC$time <= t), 0L)
    if (idx == 0L) 1 else kmC$surv[idx]
  }
  ev_times <- sort(unique(time[status == 1L]))
  if (!length(ev_times)) return(numeric(0))
  cap <- stats::quantile(ev_times, cap_quantile, type = 2, na.rm = TRUE)
  grid <- ev_times[ev_times <= cap]
  if (!length(grid)) grid <- ev_times
  ok <- vapply(grid, function(t) {
    ghat_t   <- getG(t)
    n_ev_t   <- sum(status == 1L & time <= t)
    n_risk_t <- sum(time >= t)
    is.finite(ghat_t) && (ghat_t >= ghat_threshold) &&
      (n_ev_t >= min_events) && (n_risk_t >= min_at_risk)
  }, logical(1))
  unique(grid[ok])
}

build_time_support_table <- function(DT, candidate_times, time_col="TIME", status_col="PHENO",
                                     min_events=10L, min_at_risk=50L, ghat_threshold=0.02) {
  time   <- as.numeric(DT[[time_col]])
  status <- as.integer(DT[[status_col]])
  kmC <- survival::survfit(survival::Surv(time, 1L - status) ~ 1)
  getG <- function(t) {
    idx <- max(which(kmC$time <= t), 0L)
    if (idx == 0L) 1 else kmC$surv[idx]
  }
  cand <- sort(unique(candidate_times))
  if (!length(cand)) return(data.table())
  out <- lapply(cand, function(t) {
    ghat_t   <- suppressWarnings(getG(t))
    n_ev_t   <- sum(status == 1L & time <= t)
    n_risk_t <- sum(time >= t)
    pass <- is.finite(ghat_t) && (ghat_t >= ghat_threshold) &&
      (n_ev_t >= min_events) && (n_risk_t >= min_at_risk)
    data.table(
      time_years = as.numeric(t),
      ghat = as.numeric(ghat_t),
      n_events_by_t = n_ev_t,
      n_at_risk = n_risk_t,
      ghat_threshold = ghat_threshold,
      min_events = min_events,
      min_at_risk = min_at_risk,
      pass = pass
    )
  })
  rbindlist(out)
}

audit_requested_time_points <- function(DT, req_times,
                                        time_col="TIME", status_col="PHENO",
                                        min_events=10L, min_at_risk=50L, ghat_threshold=0.02) {
  if (length(req_times) == 0) return(list(ghat_table=data.table(), valid_times=numeric(0)))
  time   <- as.numeric(DT[[time_col]])
  status <- as.integer(DT[[status_col]])
  kmC <- survival::survfit(survival::Surv(time, 1L - status) ~ 1)
  getG <- function(t) {
    idx <- max(which(kmC$time <= t), 0L)
    if (idx == 0L) 1 else kmC$surv[idx]
  }
  out <- lapply(req_times, function(t) {
    ghat_t   <- suppressWarnings(getG(t))
    n_ev_t   <- sum(status == 1L & time <= t)
    n_risk_t <- sum(time >= t)
    # Status mirrors the FULL horizon guard (G-hat AND at-risk AND events), so a
    # horizon that will be dropped for too few events/at-risk is not shown as PASS
    # once the count columns are stripped for disclosure in the published copy.
    status_flag <- if (!is.finite(ghat_t)) "No estimate"
    else if (ghat_t < ghat_threshold) sprintf("FAIL (G-hat<%g)", ghat_threshold)
    else if (n_risk_t < min_at_risk)   sprintf("FAIL (at_risk<%d)", min_at_risk)
    else if (n_ev_t < min_events)      sprintf("FAIL (events<%d)", min_events)
    else "PASS"
    data.table(time_years=t, ghat=as.numeric(ghat_t), n_events_by_t=n_ev_t,
               n_at_risk=n_risk_t, status=status_flag)
  })
  gtab <- rbindlist(out)
  valid <- gtab[status=="PASS" & is.finite(ghat) & ghat >= ghat_threshold &
                  n_events_by_t >= min_events & n_at_risk >= min_at_risk, time_years]
  list(ghat_table=gtab, valid_times=as.numeric(valid))
}

safe_formula <- function(response, predictors, strata=NULL) {
  pred_str <- if (length(predictors) > 1) paste(predictors, collapse=" + ") else predictors
  if (!is.null(strata) && nchar(strata) > 0) as.formula(sprintf("%s ~ %s + strata(%s)", response, pred_str, strata))
  else as.formula(sprintf("%s ~ %s", response, pred_str))
}

calculate_harrells_ci <- function(model, data, n_boot = 1000, seed = 42, max_pairs_for_bootstrap = 20000) {
  # n_case_control_pairs = events x non-events: an UPPER BOUND on comparable
  # survival pairs (censoring makes many pairs non-comparable). Used only to
  # decide bootstrap-vs-SE CI, not as a concordance-pairs count.
  n <- nrow(data); n_events <- sum(data$PHENO == 1); n_case_control_pairs <- n_events * (n - n_events)
  n_pairs <- n_case_control_pairs  # legacy alias kept in the returned list
  smry <- summary(model)
  c_index <- smry$concordance[1]; c_se <- smry$concordance[2]
  cat(sprintf("    N=%s, Events=%s, Case-control pairs=%s\n",
              format(n, big.mark=","), format(n_events, big.mark=","), format(n_case_control_pairs, big.mark=",")))
  if (n_pairs <= max_pairs_for_bootstrap) {
    cat("    Using bootstrap CI (comparison pairs ≤ 20,000)...\n")
    set.seed(seed); boot_cs <- rep(NA_real_, n_boot); n_success <- 0
    for (i in seq_len(n_boot)) {
      boot_idx <- sample.int(n, n, replace = TRUE); boot_data <- data[boot_idx, ]
      boot_cs[i] <- tryCatch({ summary(update(model, data=boot_data))$concordance[1] }, error=function(e) NA_real_)
      if (is.finite(boot_cs[i])) n_success <- n_success + 1
    }
    boot_cs <- boot_cs[is.finite(boot_cs)]
    if (length(boot_cs) < 10) {
      ci_lower <- c_index - 1.96*c_se; ci_upper <- c_index + 1.96*c_se; method <- "SE (bootstrap failed)"
    } else {
      ci_lower <- quantile(boot_cs, 0.025); ci_upper <- quantile(boot_cs, 0.975); method <- "bootstrap"
      cat(sprintf("    Bootstrap success: %d/%d iterations\n", n_success, n_boot))
    }
  } else {
    cat("    Using SE-based CI (comparison pairs > 20,000)...\n")
    ci_lower <- c_index - 1.96*c_se; ci_upper <- c_index + 1.96*c_se; method <- "SE"
  }
  ci_lower <- max(0, ci_lower); ci_upper <- min(1, ci_upper)
  cat(sprintf("    Harrell's C: %.4f (95%% CI: %.4f-%.4f) [%s]\n", c_index, ci_lower, ci_upper, method))
  list(c_index=c_index, ci_lower=ci_lower, ci_upper=ci_upper, se=c_se,
       n_case_control_pairs=n_case_control_pairs, n_pairs=n_pairs, method=method)
}

calculate_uno_cindex <- function(model, data, time_col="TIME", status_col="PHENO",
                                 min_events=10L, min_at_risk=50L, ghat_threshold=0.02,
                                 skip_ci=FALSE, n_boot=500L, seed=42L, tau=NA_real_) {
  if (!requireNamespace("survAUC", quietly = TRUE)) {
    stop("Package 'survAUC' is required for Uno's C (install.packages('survAUC')).")
  }
  time <- as.numeric(data[[time_col]]); status <- as.integer(data[[status_col]])
  if (length(unique(time)) < 2L || sum(status == 1L) < 5L) {
    return(list(cindex=NA_real_, se=NA_real_, ci_lower=NA_real_, ci_upper=NA_real_, method="insufficient_support", tau=NA_real_))
  }
  lp <- tryCatch(as.numeric(predict(model, newdata=data, type="lp")), error=function(e) NA_real_)
  if (!is.finite(sd(lp)) || all(!is.finite(lp))) {
    return(list(cindex=NA_real_, se=NA_real_, ci_lower=NA_real_, ci_upper=NA_real_, method="no_variation", tau=NA_real_))
  }
  set.seed(1L); lp <- lp + rnorm(length(lp), sd=1e-12)
  kmC <- survival::survfit(survival::Surv(time, 1L - status) ~ 1)
  get_Ghat <- function(t){ idx <- max(which(kmC$time <= t), 0L); if (idx==0L) 1 else kmC$surv[idx] }
  pick_tau <- function() {
    ev_times <- sort(unique(time[status == 1L])); if (!length(ev_times)) return(NA_real_)
    cap95 <- stats::quantile(ev_times, 0.95, type=2, na.rm=TRUE); grid <- ev_times[ev_times <= cap95]; if (!length(grid)) grid <- ev_times
    ok <- vapply(grid, function(t){ ghat_t <- get_Ghat(t); n_ev <- sum(status==1 & time<=t); n_risk <- sum(time>=t);
    is.finite(ghat_t) && ghat_t >= ghat_threshold && n_ev >= min_events && n_risk >= min_at_risk }, logical(1))
    if (!any(ok)) return(NA_real_); max(grid[ok])
  }
  tau_val <- if (is.finite(tau)) tau else pick_tau()
  if (!is.finite(tau_val) || tau_val <= 0) {
    ev_times <- time[status==1]; if (!length(ev_times)) return(list(cindex=NA_real_, se=NA_real_, ci_lower=NA_real_, ci_upper=NA_real_, method="failed_tau", tau=NA_real_))
    tau_val <- stats::median(ev_times, na.rm=TRUE)
  }
  uno_c <- tryCatch({
    survAUC::UnoC(Surv.rsp=survival::Surv(time, status),
                  Surv.rsp.new=survival::Surv(time, status),
                  lpnew=lp, time=tau_val)
  }, error=function(e) NA_real_)
  if (!is.finite(uno_c) || uno_c <= 0 || uno_c > 1) return(list(cindex=NA_real_, se=NA_real_, ci_lower=NA_real_, ci_upper=NA_real_, method="failed_estimation", tau=tau_val))
  se <- ci_lower <- ci_upper <- NA_real_; method <- sprintf("UnoC_tau=%.4f", tau_val)
  if (!isTRUE(skip_ci)) {
    set.seed(seed); B <- as.integer(n_boot); boot_vals <- rep(NA_real_, B)
    idx_e <- which(status==1L); idx_c <- which(status==0L)
    for (b in seq_len(B)) {
      sb <- c(sample(idx_e, length(idx_e), replace=TRUE), sample(idx_c, length(idx_c), replace=TRUE))
      d_b <- data[sb, ]; m_b <- tryCatch(update(model, data=d_b), error=function(e) NULL); if (is.null(m_b)) next
      lp_b <- tryCatch(as.numeric(predict(m_b, newdata=d_b, type="lp")), error=function(e) NULL); if (is.null(lp_b)) next
      lp_b <- lp_b + rnorm(length(lp_b), sd=1e-12)
      t_b <- min(tau_val, max(d_b[[time_col]][d_b[[status_col]]==1L], na.rm=TRUE))
      uno_b <- tryCatch({
        survAUC::UnoC(Surv.rsp=survival::Surv(d_b[[time_col]], d_b[[status_col]]),
                      Surv.rsp.new=survival::Surv(d_b[[time_col]], d_b[[status_col]]),
                      lpnew=lp_b, time=t_b)
      }, error=function(e) NA_real_)
      if (is.finite(uno_b) && uno_b>0 && uno_b<=1) boot_vals[b] <- uno_b
    }
    boot_vals <- boot_vals[is.finite(boot_vals)]
    if (length(boot_vals) >= 50L) {
      se <- stats::sd(boot_vals); ci_lower <- stats::quantile(boot_vals, 0.025); ci_upper <- stats::quantile(boot_vals, 0.975)
      uno_c <- mean(boot_vals); method <- paste0(method, "_boot")
    }
  }
  list(cindex=uno_c, se=se, ci_lower=ci_lower, ci_upper=ci_upper, method=method, tau=tau_val)
}

extract_comprehensive_stats <- function(fit, model_name) {
  smry <- summary(fit)
  coef_table <- data.table(
    variable = rownames(smry$coefficients),
    coef = smry$coefficients[, "coef"],
    exp_coef = smry$coefficients[, "exp(coef)"],
    se_coef = smry$coefficients[, "se(coef)"],
    z = smry$coefficients[, "z"],
    p = smry$coefficients[, "Pr(>|z|)"],
    HR = smry$coefficients[, "exp(coef)"],
    CI_lower = exp(smry$coefficients[, "coef"] - 1.96 * smry$coefficients[, "se(coef)"]),
    CI_upper = exp(smry$coefficients[, "coef"] + 1.96 * smry$coefficients[, "se(coef)"])
  )
  coef_table[, model := model_name]
  fit_stats <- data.table(
    model = model_name,
    n = fit$n,
    n_events = fit$nevent,
    loglik_null = fit$loglik[1],
    loglik_model = fit$loglik[2],
    AIC = AIC(fit),
    BIC = BIC(fit),
    concordance = smry$concordance[1],
    concordance_se = smry$concordance[2],
    rsquare = smry$rsq[1],
    max_rsquare = smry$rsq[2],
    wald_test = smry$waldtest[1],
    wald_p = smry$waldtest[3],
    logrank_test = smry$sctest[1],
    logrank_p = smry$sctest[3]
  )
  list(coefficients=coef_table, fit_stats=fit_stats, summary=smry, model_name=model_name)
}

# ==============================================================================
# LOADING DATA
# ==============================================================================
cat("==============================================================================\nLOADING DATA\n==============================================================================\n\n")
cat("1. Loading PRS scores...\n")
if (!file.exists(opt$prs_file)) stop("PRS file not found: ", opt$prs_file)
prs_all <- fread(opt$prs_file)
cat("   Loaded", format(nrow(prs_all), big.mark=","), "PRS scores\n")

# Population subset
if (!is.null(opt$pop_file) && file.exists(opt$pop_file)) {
  ids <- fread(opt$pop_file, header=FALSE)
  prs_all <- prs_all[get(opt$id_col) %in% ids$V2]
}

# Validate PRS column
if (!opt$prs_col %in% names(prs_all)) {
  score_cols <- grep("SCORE.*SUM|SCORE1_SUM|SCORESUM", names(prs_all), value=TRUE, ignore.case=TRUE)
  if (length(score_cols) > 0) opt$prs_col <- score_cols[1] else stop("Cannot find PRS score column")
}
assert_unique_iid(prs_all, opt$id_col, "primary PRS file")
prs <- prs_all[, c(opt$id_col, opt$prs_col), with=FALSE]
setnames(prs, opt$prs_col, "SCORESUM")
# Standardise the ID column to "IID" so the later merge-by-IID works for a
# non-default --id_col (previously the column kept its original name and the
# merge would fail).
if (opt$id_col != "IID" && opt$id_col %in% names(prs)) setnames(prs, opt$id_col, "IID")
record_stage(prs, "primary_PRS_loaded")

# Load phenotype
cat("2. Loading phenotype...\n")
phen <- fread(opt$pheno_file)
if (!"IID" %in% names(phen)) {
  if (opt$id_col %in% names(phen)) setnames(phen, opt$id_col, "IID") else setnames(phen, names(phen)[1], "IID")
}
assert_unique_iid(phen, "IID", "phenotype file")

# ------------------------------------------------------------------------------
# Decide how STOP is obtained.
#   Preferred: two AGE columns -> STOP = age_exit - age_t1 (age_t1 also the covariate)
#   Fallback : an explicit --time_exit column (back-compat)
# ------------------------------------------------------------------------------
# If an age-exit column is EXPLICITLY requested but absent, stop rather than
# silently falling back to --time_exit (a different time scale — dangerous for a
# harmonised multi-biobank analysis).
if (!is.null(opt$age_exit_col) && nzchar(opt$age_exit_col) && !(opt$age_exit_col %in% names(phen)))
  stop("--age_exit_col='", opt$age_exit_col, "' not found in the phenotype file.\n",
       "  Provide the column, or clear --age_exit_col to use the --time_exit time scale.")
use_age_derivation <- !is.null(opt$age_exit_col) && nzchar(opt$age_exit_col) &&
                      !is.null(opt$age_t1_col)   && nzchar(opt$age_t1_col)   &&
                      (opt$age_t1_col %in% names(phen)) && (opt$age_exit_col %in% names(phen))
age_t1_present <- !is.null(opt$age_t1_col) && nzchar(opt$age_t1_col) && (opt$age_t1_col %in% names(phen))

# Columns to pull from phenotype (status + time + age + descriptive)
ph_cols <- c("IID", opt$status_col)
if (use_age_derivation) {
  ph_cols <- unique(c(ph_cols, opt$age_t1_col, opt$age_exit_col))
} else {
  ph_cols <- unique(c(ph_cols, opt$time_exit))
  if (age_t1_present) ph_cols <- unique(c(ph_cols, opt$age_t1_col))  # still want age covariate
}
extra_sum_cols <- character(0)
if (!is.null(opt$time_summary_cols) && nzchar(opt$time_summary_cols)) {
  extra_sum_cols <- trimws(strsplit(opt$time_summary_cols, ",")[[1]])
  extra_sum_cols <- extra_sum_cols[nzchar(extra_sum_cols)]
  ph_cols <- unique(c(ph_cols, extra_sum_cols))
}
if (!is.null(opt$time_entry) && nzchar(opt$time_entry)) ph_cols <- unique(c(ph_cols, opt$time_entry))

prs <- merge(prs, phen[, intersect(ph_cols, names(phen)), with=FALSE], by="IID")

# Standardize PHENO
setnames(prs, opt$status_col, "PHENO")

# Keep the age-at-T1 column under a stable name AGE_T1 (covariate + index age)
AGE_T1_AVAILABLE <- age_t1_present
if (AGE_T1_AVAILABLE) prs[, AGE_T1 := as.numeric(get(opt$age_t1_col))]

# Derive or take START/STOP
if (use_age_derivation) {
  prs[, AGE_EXIT := as.numeric(get(opt$age_exit_col))]
  prs[, START := 0]
  prs[, STOP  := AGE_EXIT - AGE_T1]
  cat(sprintf("\n  Time scale derived from ages: STOP = %s - %s (START = 0)\n",
              opt$age_exit_col, opt$age_t1_col))
} else {
  setnames(prs, opt$time_exit, "STOP")
  if (!is.null(opt$time_entry) && nzchar(opt$time_entry)) setnames(prs, opt$time_entry, "START") else prs[, START := 0]
  # Reconstruct age-at-exit for the descriptive summaries (age-at-T2, interval) when
  # only a time column is supplied but age-at-T1 is present: AGE_EXIT = AGE_T1 + STOP.
  if (AGE_T1_AVAILABLE) prs[, AGE_EXIT := AGE_T1 + STOP]
  cat("\n  Time scale taken from --time_exit column:", opt$time_exit, "\n")
}

# Drop incomplete PHENO/STOP
prs <- prs[!is.na(PHENO) & !is.na(STOP)]
# Recode 1/2 -> 0/1, but ONLY when BOTH codes are present. The old
# all(unique %in% c(1,2)) test was TRUE for an all-events stratum (only 1s),
# silently converting every event to a censor. setequal requires both values.
.u <- sort(unique(prs$PHENO))
if (setequal(.u, c(1, 2))) {
  prs[, PHENO := PHENO - 1L]
} else if (!all(.u %in% c(0, 1))) {
  stop("Event status must be coded 0/1 or 1/2; got: ", paste(.u, collapse=", "))
}
record_stage(prs, "phenotype_matched")

# ==============================================================================
# START–STOP ENFORCEMENT & FOLLOW-UP
# ==============================================================================
cat("\nUSING START–STOP: Surv(START, STOP, PHENO)\n")

# TIME_SINCE_T1 = time from T1 to T2/censor, preserved as its own column so lag
# and prospective re-indexing (which overwrite STOP) can still refer to it.
if (exists("AGE_T1_AVAILABLE") && isTRUE(AGE_T1_AVAILABLE) && "AGE_EXIT" %in% names(prs)) {
  prs[, TIME_SINCE_T1 := AGE_EXIT - AGE_T1]
} else {
  prs[, TIME_SINCE_T1 := STOP]   # gwas-aligned fallback: STOP already = time since T1 (START=0)
}

# Load the recruitment-age file if supplied, in ANY mode (prospective needs it for
# re-indexing; gwas_aligned uses it only for the descriptive recruitment-timing QC).
if (!is.null(opt$recruit_file) && file.exists(opt$recruit_file)) {
  .rec <- fread(opt$recruit_file)
  if (!"IID" %in% names(.rec)) setnames(.rec, names(.rec)[1], "IID")
  assert_unique_iid(.rec, "IID", "recruitment file")
  if (!opt$recruit_age_col %in% names(.rec))
    stop("Column '", opt$recruit_age_col, "' not found in --recruit_file")
  prs <- merge(prs, .rec[, c("IID", opt$recruit_age_col), with=FALSE], by="IID", all.x=TRUE)
  prs[, AGE_RECRUIT := as.numeric(get(opt$recruit_age_col))]
  # COVERAGE. This is a LEFT join: IIDs absent from the recruitment file get
  # AGE_RECRUIT = NA. In prospective mode every NA row is unusable, and it used to
  # vanish silently (NA propagated through timing_category, and data.table drops NA
  # rows on a `!=` filter). A partial recruitment file could therefore redefine the
  # analysis cohort without saying so. Report it loudly in EVERY mode.
  recruit_coverage_total <- nrow(prs)
  recruit_coverage_n     <- sum(!is.na(prs$AGE_RECRUIT))
  recruit_coverage_pct   <- if (recruit_coverage_total > 0)
                              100 * recruit_coverage_n / recruit_coverage_total else NA_real_
  cat(sprintf("  Recruitment-age coverage: %d of %d (%.1f%%)\n",
              recruit_coverage_n, recruit_coverage_total, recruit_coverage_pct))
  if (recruit_coverage_n < recruit_coverage_total) {
    cat("  ⚠ ", recruit_coverage_total - recruit_coverage_n,
        " individual(s) have NO recruitment age.\n", sep="")
    if (identical(opt$analysis_mode, "prospective"))
      cat("    In prospective mode these CANNOT be indexed and are excluded below.\n",
          "    Check that --recruit_file covers the same individuals as the phenotype file.\n", sep="")
    else
      cat("    gwas_aligned: recruitment age is descriptive only; no rows are excluded for this.\n")
  }
}

# Snapshot for the temporal-order QC, taken BEFORE any strict-order/lag/prospective
# exclusion so same-day (TIME_SINCE_T1==0) transitions are counted.
timing_qc_snapshot <- prs[, intersect(c("IID","PHENO","TIME_SINCE_T1","AGE_T1","AGE_EXIT","AGE_RECRUIT"), names(prs)), with=FALSE]

# Prospective mode and a disease-duration lag define different time origins;
# combining them is ambiguous. Disallow.
if (identical(opt$analysis_mode, "prospective") && is.finite(opt$lag_days) && opt$lag_days > 0)
  stop("--analysis_mode=prospective cannot be combined with --lag_days>0 (ambiguous time origin).\n",
       "  Run the prospective analysis and the lag sensitivity as separate jobs.")

# ------------------------------------------------------------------------------
# PROSPECTIVE RE-INDEXING (Layer 2) — before the STOP>START filter
# ------------------------------------------------------------------------------
# Requires the age-derived scale + a recruitment-age file. Index differs by
# whether T1 is prevalent (index = recruitment) or incident (index = T1[+lag]).
timing_category <- NULL
if (identical(opt$analysis_mode, "prospective")) {
  if (!AGE_T1_AVAILABLE) stop("prospective mode requires --age_t1_col present in the phenotype file")
  if (!"AGE_RECRUIT" %in% names(prs))   # loaded above from --recruit_file
    stop("prospective mode requires --recruit_file (IID + ", opt$recruit_age_col, ")")
  cat("  ANALYSIS MODE: prospective (Layer 2)\n")
  if (!"AGE_EXIT" %in% names(prs)) prs[, AGE_EXIT := AGE_T1 + STOP]  # reconstruct if fallback scale

  lag_yr <- opt$incident_lag_days / 365.25

  # Drop individuals with no recruitment age EXPLICITLY and FIRST. They cannot be
  # indexed, so they have no place in a prospective analysis — but until this was
  # added they were removed silently: is_incident became NA, timing_category became
  # NA, and `prs[timing_category != "Historical"]` below drops NA rows without a
  # word (data.table treats NA as FALSE in i). A recruitment file covering only
  # part of the cohort therefore shrank the sample invisibly.
  n_no_recruit <- sum(is.na(prs$AGE_RECRUIT))
  if (n_no_recruit > 0) {
    cat("  Excluding ", n_no_recruit, " of ", nrow(prs),
        " individual(s) with no recruitment age (cannot be indexed)\n", sep="")
    prs <- prs[!is.na(AGE_RECRUIT)]
    record_stage(prs, "recruit_age_available")
    if (!nrow(prs))
      stop("No individual has a recruitment age; prospective mode cannot proceed.\n",
           "  Check that --recruit_file uses the same IIDs as the phenotype file.")
  }

  prs[, is_incident := AGE_T1 > AGE_RECRUIT]
  # timing category
  prs[, timing_category := fifelse(!is_incident & PHENO == 1 & AGE_EXIT <= AGE_RECRUIT, "Historical",
                             fifelse(!is_incident, "Prevalent_future",
                                     "Incident"))]
  # Backstop: every row must now classify. An NA here means a new path produced a
  # missing AGE_RECRUIT/AGE_EXIT after the exclusion above, and would again be
  # dropped silently by the `!= "Historical"` filter.
  if (any(is.na(prs$timing_category)))
    stop(sum(is.na(prs$timing_category)), " row(s) have an NA timing_category after ",
         "recruitment-age filtering; refusing to drop them silently.")
  timing_category <- prs$timing_category
  # Record the full classification (before excluding Historical) for the manuscript.
  timing_tab <- prs[, .N, by=timing_category][order(timing_category)]
  # Never publish a blank-named category (an NA bucket used to leak into this file).
  timing_tab <- timing_tab[!is.na(timing_category)]
  timing_console <- copy(timing_tab)   # unsuppressed, for the local console only
  # Disclosure: blank a category count below the threshold (the internal Historical
  # exclusion below keys off the prs$timing_category column, not this written table).
  suppress_partition_all_or_none(timing_tab, count_cols = "N", threshold = opt$min_cell_count)
  fwrite(timing_tab, file.path(opt$outdir, paste0(prefix, "timing_categories.csv")))
  cat("  Timing categories:\n"); for (i in seq_len(nrow(timing_console)))
    cat(sprintf("    %-18s %d\n", timing_console$timing_category[i], timing_console$N[i]))
  # Exclude Historical (not free of T2 at index)
  n_hist <- sum(prs$timing_category == "Historical", na.rm=TRUE)
  if (n_hist > 0) cat("  Excluding", n_hist, "Historical trajectories (T2 before recruitment)\n")
  prs <- prs[timing_category != "Historical"]
  record_stage(prs, "prospective_index_T2free")
  # Optionally restrict T1 to incident as well. Default 'both' keeps prevalent T1
  # (indexed at recruitment, separated by strata(T1_STATUS) + T1_DURATION), which is
  # the long-standing behaviour and leaves every existing run unchanged.
  if (identical(opt$prospective_t1, "incident")) {
    n_prev <- sum(prs$timing_category == "Prevalent_future", na.rm=TRUE)
    cat("  --prospective_t1=incident: excluding ", n_prev,
        " Prevalent_future trajectory(ies); T1 must occur after recruitment\n", sep="")
    prs <- prs[timing_category == "Incident"]
    record_stage(prs, "prospective_incident_T1_only")
    if (!nrow(prs))
      stop("--prospective_t1=incident leaves no individuals (no incident T1 in this cohort).")
  } else if (!identical(opt$prospective_t1, "both")) {
    stop("--prospective_t1 must be 'both' or 'incident'; got '", opt$prospective_t1, "'.")
  }
  # Re-index: index = recruitment (prevalent) or T1[+incident lag] (incident).
  # AGE_INDEX = current age at prediction start; T1_DURATION = time already spent
  # in the T1 state at the index (0(+lag) for incident, recruit-T1 for prevalent);
  # T1_STATUS distinguishes the two (used as strata in the prospective models).
  prs[, AGE_INDEX  := fifelse(is_incident, AGE_T1 + lag_yr, AGE_RECRUIT)]
  prs[, T1_DURATION := AGE_INDEX - AGE_T1]
  prs[, T1_STATUS  := factor(fifelse(is_incident, "Incident", "Prevalent"))]
  prs[, START := 0]
  prs[, STOP  := AGE_EXIT - AGE_INDEX]
  n_nonpos <- sum(prs$STOP <= 0, na.rm=TRUE)
  if (n_nonpos > 0) { cat("  Dropping", n_nonpos, "rows with STOP<=0 after re-indexing\n"); prs <- prs[STOP > 0] }
}

# ------------------------------------------------------------------------------
# LAG SENSITIVITY (Layer 3, gwas_aligned only) — LEFT-TRUNCATION at the lag.
# Keep only transitions longer than the lag, then start the at-risk clock at the
# lag (time origin stays T1). This conditions on being event-free through the
# washout without changing the time origin.
# ------------------------------------------------------------------------------
if (is.finite(opt$lag_days) && opt$lag_days > 0) {
  lag_yr <- opt$lag_days / 365.25
  n_lag <- sum(prs$TIME_SINCE_T1 <= lag_yr, na.rm=TRUE)
  cat(sprintf("  Lag sensitivity (left-truncation at %g d): excluding %d rows with time_since_T1 <= lag\n",
              opt$lag_days, n_lag))
  prs <- prs[TIME_SINCE_T1 > lag_yr]
  prs[, START := lag_yr]
  prs[, STOP  := TIME_SINCE_T1]
  record_stage(prs, "lag_applied")
}

# Same-age T1/T2 exclusion, at survival's OWN time resolution. Two disjoint sets:
#   n_invalid    - STOP <= START, i.e. not a strictly ordered transition;
#   n_degenerate - ordered, but rounds to zero length under aeqSurv (float
#                  cancellation between two ages derived from the same date).
# Testing STOP - START rather than STOP == 0 also fixes a reporting bug: with a lag
# sweep START = lag_yr > 0, so STOP == 0 is unreachable and the same-age count was
# previously printed as 0 regardless of the truth.
resolvable   <- .aeqsurv_keep(prs$START, prs$STOP)
n_invalid    <- sum(prs$STOP <= prs$START, na.rm=TRUE)
n_degenerate <- sum(prs$STOP > prs$START & !resolvable, na.rm=TRUE)
n_sameage    <- n_invalid + n_degenerate
n_sameage_events <- sum((prs$STOP <= prs$START | !resolvable) & prs$PHENO == 1, na.rm=TRUE)
min_kept_followup <- suppressWarnings(min((prs$STOP - prs$START)[prs$STOP > prs$START & resolvable], na.rm=TRUE))
cat("  Same-age T1/T2 (strict-order excluded): ", n_sameage,
    " (of which ", n_degenerate, " collapse only at survival's time resolution)\n", sep="")
# The progression GWAS counted same-day T1->T2 as cases; this pipeline excludes
# them. Report how many events that costs so the discrepancy is never invisible.
cat("    ...of which events: ", n_sameage_events, "\n", sep="")
cat("    smallest RETAINED follow-up (STOP - START): ",
    if (is.finite(min_kept_followup)) format(min_kept_followup) else "n/a", "\n", sep="")
cat("  Removing rows with STOP <= START: ", n_invalid, "\n", sep="")
# Data-quality flag: integer ages collapse sub-year transitions to STOP=0 and drop
# them, which can be substantial for acute transitions (MI->HF, AF->pacemaker).
ages_are_integer <- FALSE
if (exists("AGE_T1_AVAILABLE") && isTRUE(AGE_T1_AVAILABLE) && "AGE_T1" %in% names(prs)) {
  a_ok <- prs$AGE_T1[is.finite(prs$AGE_T1)]
  # scalar test (the previous version was a length-n vector, so isTRUE was always FALSE)
  ages_are_integer <- length(a_ok) > 0 && all(abs(a_ok - round(a_ok)) < 1e-6)
  if (ages_are_integer && n_sameage > 0)
    cat("  ⚠ WARNING: ages appear INTEGER-VALUED and ", n_sameage,
        " same-age transitions were dropped. Integer ages collapse sub-year T1->T2\n",
        "    intervals (e.g. MI->HF); use fractional ages or dates matching the GWAS temporal rule.\n", sep="")
}
prs <- prs[STOP > START & resolvable]
prs[, FOLLOWUP := STOP - START]

# Post-recode integrity checks (harmonised analysis must not proceed on malformed
# survival columns).
if (!all(prs$PHENO %in% c(0L, 1L)))
  stop("Event status (PHENO) is not strictly 0/1 after recoding; check --status_col.")
if (!all(is.finite(prs$START)) || !all(is.finite(prs$STOP)))
  stop("Non-finite START/STOP after time-scale derivation; check the age/time columns.")
record_stage(prs, "valid_timing")

# ==============================================================================
# LOAD covariates (optional)
# ==============================================================================
mycovs <- character(0); clinical_covs <- character(0); dropped_covs <- character(0)
# Strict-by-default: a requested covariate absent from the file STOPS the run
# (harmonisation across biobanks). --allow_missing_covariates downgrades to a warn.
check_requested <- function(requested, available, what) {
  miss <- setdiff(requested, available)
  if (length(miss)) {
    msg <- paste0(what, " requested but not found in --cov_file: ", paste(miss, collapse=", "))
    if (isTRUE(opt$allow_missing_covariates)) {
      warning(msg, " — proceeding without them (--allow_missing_covariates).")
      cat("   ⚠ ", msg, " — dropped (--allow_missing_covariates)\n", sep="")
      dropped_covs <<- c(dropped_covs, miss)
    } else {
      stop(msg, "\n  For a harmonised multi-biobank analysis, base models must match across sites.\n",
           "  Add the column, or pass --allow_missing_covariates to proceed without it.")
    }
  }
  intersect(requested, available)
}
if (!is.null(opt$cov_file) && file.exists(opt$cov_file)) {
  cat("\n3. Loading covariates...\n")
  covf <- fread(opt$cov_file)
  if (!"IID" %in% names(covf)) setnames(covf, names(covf)[1], "IID")
  assert_unique_iid(covf, "IID", "covariate file")
  if (opt$covariates != "") {
    requested <- trimws(strsplit(opt$covariates, ",")[[1]]); requested <- requested[nzchar(requested)]
    mycovs <- check_requested(requested, names(covf), "Basic covariate(s)")
    cat("   Basic covariates:", ifelse(length(mycovs)>0, paste(mycovs, collapse=", "), "None"), "\n")
  }
  if (!is.null(opt$clinical_covariates) && opt$clinical_covariates != "") {
    requested_clin <- trimws(strsplit(opt$clinical_covariates, ",")[[1]]); requested_clin <- requested_clin[nzchar(requested_clin)]
    clinical_covs <- check_requested(requested_clin, names(covf), "Clinical covariate(s)")
    if (length(clinical_covs) > 0) cat("   Clinical risk factors:", paste(clinical_covs, collapse=", "), "\n")
  }
  all_covs <- unique(c(mycovs, clinical_covs))
  if (length(all_covs) > 0) prs <- merge(prs, covf[, c("IID", all_covs), with=FALSE], by="IID", all.x=TRUE)
} else if (opt$covariates != "") {
  if (isTRUE(opt$allow_missing_covariates)) cat("\n3. ⚠ WARNING: Covariates specified but covariate file not found\n")
  else stop("Covariates specified but --cov_file not found: ", opt$cov_file,
            "\n  Pass --allow_missing_covariates to proceed without covariates.")
}

# Complete cases on REQUIRED MODEL COLUMNS ONLY (never descriptive/time_summary
# columns): a missing descriptive value must not drop a participant from the model.
# CORE vs CLINICAL COMPLETE-CASE.
# clinical_covs are OPTIONAL exemplar variables. Including them here would let a
# missing clinical value drop a participant from the systematic Core M0-M7
# analysis, so simply supplying --clinical_covariates could silently shrink the
# manuscript's primary sample. The core filter therefore uses core columns only;
# the clinical ladder applies its own additional filter downstream.
req_cols <- c("PHENO", "START", "STOP", mycovs)
if (exists("AGE_T1_AVAILABLE") && isTRUE(AGE_T1_AVAILABLE) && "AGE_T1" %in% names(prs)) req_cols <- c(req_cols, "AGE_T1")
if (identical(opt$analysis_mode, "prospective"))
  req_cols <- c(req_cols, intersect(c("AGE_INDEX","T1_DURATION","T1_STATUS"), names(prs)))
req_cols <- intersect(unique(req_cols), names(prs))
n_before_cc <- nrow(prs)
prs <- prs[complete.cases(prs[, ..req_cols])]
if (n_before_cc - nrow(prs) > 0)
  cat("   Complete-case on CORE model columns only: dropped", n_before_cc - nrow(prs), "rows\n")
if (length(clinical_covs) > 0) {
  clin_present <- intersect(clinical_covs, names(prs))
  n_clin_complete <- if (length(clin_present)) sum(complete.cases(prs[, ..clin_present])) else nrow(prs)
  cat("   Clinical covariates are NOT part of the core filter; the Clinical ladder will use ",
      n_clin_complete, " of ", nrow(prs), " core rows.\n", sep="")
}
record_stage(prs, "complete_covariates")

# ==============================================================================
# FOLLOW-UP QC with configurable threshold (FOLLOWUP)
# ==============================================================================
cat("\n4. Checking for very short follow-up (FOLLOWUP < ", opt$min_followup_threshold, "):\n", sep="")
n_very_short <- sum(prs$FOLLOWUP < opt$min_followup_threshold, na.rm=TRUE)
if (n_very_short > 0) {
  cat("   ⚠ Removing", n_very_short, "rows with follow-up <", opt$min_followup_threshold, "\n")
  prs <- prs[FOLLOWUP >= opt$min_followup_threshold]
} else {
  cat("   ✓ None\n")
}

cat("\n5. Checking for outliers (very long follow-up):\n")
fu_99th <- quantile(prs$FOLLOWUP, 0.99, na.rm=TRUE); fu_max <- max(prs$FOLLOWUP, na.rm=TRUE)
cat("   99th percentile FOLLOWUP:", round(fu_99th, 2), "\n")
cat("   Maximum FOLLOWUP:", round(fu_max, 2), "\n")
if (fu_max > 3 * fu_99th) {
  cat("   ⚠ Max FOLLOWUP > 3×99th; N with FOLLOWUP > 2×99th:", sum(prs$FOLLOWUP > 2 * fu_99th, na.rm=TRUE), "\n")
} else cat("   ✓ No extreme outliers detected\n")

# ==============================================================================
# TIME-TO-EVENT DESCRIPTIVE SUMMARY (for raw columns), NO reverse-KM here
# ==============================================================================
cat("\n==============================================================================\nTIME-TO-EVENT SUMMARY (descriptive only for columns)\n==============================================================================\n\n")

calc_time_summary <- function(time_vec, label) {
  data.table(
    Group = label,
    N = length(time_vec),
    Mean = mean(time_vec, na.rm=TRUE),
    Median = median(time_vec, na.rm=TRUE),
    SD = sd(time_vec, na.rm=TRUE),
    Min = min(time_vec, na.rm=TRUE),
    Max = max(time_vec, na.rm=TRUE),
    Q25 = quantile(time_vec, 0.25, na.rm=TRUE),
    Q75 = quantile(time_vec, 0.75, na.rm=TRUE)
  )
}

time_summary <- data.table()
# Determine which columns are actually available for descriptive summary *after* renaming
actual_desc_cols <- character(0)

# 'STOP' is always the renamed --time_exit
actual_desc_cols <- c(actual_desc_cols, "STOP")

# 'START' is the renamed --time_entry, if provided
if (!is.null(opt$time_entry) && nzchar(opt$time_entry)) {
  actual_desc_cols <- c(actual_desc_cols, "START")
}

# Add any user-specified extra columns from --time_summary_cols. These retain their original names.
# Only include those that exist in the 'prs' data table after merge.
if (length(extra_sum_cols) > 0) {
  actual_desc_cols <- c(actual_desc_cols, intersect(extra_sum_cols, names(prs)))
}

desc_cols <- unique(actual_desc_cols) # Ensure uniqueness
desc_cols <- desc_cols[nzchar(desc_cols)] # Remove any potential empty strings

label_for_col <- function(col) paste0("[", col, "]")

for (col in desc_cols) {
  if (!col %in% names(prs)) { cat("  ⚠ Column not found, skipping: ", col, "\n", sep=""); next }
  vec <- prs[[col]]
  ts <- calc_time_summary(vec, paste0(label_for_col(col), " Overall"))
  time_summary <- rbind(time_summary, ts)
  cat(sprintf("%s Overall: N=%s, Mean=%.2f (SD %.2f), Median=%.2f (IQR %.2f–%.2f), Range %.2f–%.2f\n",
              label_for_col(col), format(ts$N, big.mark=","), ts$Mean, ts$SD, ts$Median, ts$Q25, ts$Q75, ts$Min, ts$Max))
  for (i in 0:1) {
    ph <- prs[PHENO == i]
    if (nrow(ph) > 0) {
      tmp <- calc_time_summary(ph[[col]], paste0(label_for_col(col), " ", if (i==0) "Censored (PHENO=0)" else "Events (PHENO=1)"))
      time_summary <- rbind(time_summary, tmp)
    }
  }
  if (!is.null(opt$sex_col) && nzchar(opt$sex_col) && opt$sex_col %in% names(prs)) {
    sex_map <- parse_sex_labels(opt$sex_labels)
    sex_values <- sort(unique(prs[[opt$sex_col]])); sex_values <- sex_values[!is.na(sex_values)]
    sex_labels <- sapply(sex_values, function(s) {
      s_str <- as.character(s)
      if (!is.null(sex_map) && s_str %in% names(sex_map)) sex_map[[s_str]]
      else if (!is.na(suppressWarnings(as.numeric(s_str)))) {
        num_str <- as.character(as.numeric(s_str))
        if (!is.null(sex_map) && num_str %in% names(sex_map)) sex_map[[num_str]] else paste0("Sex=", s)
      } else paste0("Sex=", s)
    })
    for (k in seq_along(sex_values)) {
      sd <- prs[get(opt$sex_col) == sex_values[k]]
      if (nrow(sd) > 0) {
        tmp <- calc_time_summary(sd[[col]], paste0(label_for_col(col), " ", sex_labels[k]))
        time_summary <- rbind(time_summary, tmp)
      }
    }
    for (k in seq_along(sex_values)) {
      for (j in 0:1) {
        ss <- prs[get(opt$sex_col) == sex_values[k] & PHENO == j]
        if (nrow(ss) > 0) {
          tmp <- calc_time_summary(ss[[col]], paste0(label_for_col(col), " ", sex_labels[k], " - ", if (j==0) "Censored" else "Events"))
          time_summary <- rbind(time_summary, tmp)
        }
      }
    }
  }
}
numeric_cols <- c("Mean","Median","SD","Min","Max","Q25","Q75")
time_summary[, (numeric_cols) := lapply(.SD, function(x) round(x, 2)), .SDcols=numeric_cols]
time_summary_file <- file.path(opt$outdir, paste0(prefix, "time_summary.csv"))
fwrite(time_summary, time_summary_file)
cat("\n✓ Saved time summary table:", basename(time_summary_file), "\n\n")

# ==============================================================================
# REVERSE KM FOLLOW-UP SUMMARY (FOLLOWUP = STOP - START)
# ==============================================================================
cat("==============================================================================\nREVERSE KM FOLLOW-UP SUMMARY (FOLLOWUP)\n==============================================================================\n\n")

reverse_km_quantiles <- function(time_vec, status_vec, probs = c(0.25, 0.5, 0.75)) {
  sv <- survival::survfit(survival::Surv(time_vec, 1L - status_vec) ~ 1)
  get_q <- function(p) {
    if (length(sv$time) == 0L) return(NA_real_)
    target <- 1 - p
    idx <- which(sv$surv <= target)
    if (length(idx) == 0L) return(NA_real_)
    min(sv$time[idx], na.rm = TRUE)
  }
  qs <- vapply(probs, get_q, numeric(1))
  names(qs) <- paste0("Q", sprintf("%02d", as.integer(probs*100)))
  list(q=qs, max_followup = if (length(sv$time)) max(sv$time, na.rm=TRUE) else NA_real_)
}

# One-row cohort summary on the FOLLOW-UP scale (col `TIME`/`FOLLOWUP`), reusing
# reverse_km_quantiles. Expects PHENO and a follow-up column; AGE_T1/TIME_SINCE_T1
# are optional (age-at-T1/T2 and the T1->T2 interval are added when present).
summarise_cohort <- function(dt, label, fu_col = "FOLLOWUP") {
  q <- function(x, p) if (length(x[is.finite(x)])) as.numeric(stats::quantile(x, p, na.rm=TRUE)) else NA_real_
  fu <- if (fu_col %in% names(dt)) dt[[fu_col]] else if ("TIME" %in% names(dt)) dt$TIME else rep(NA_real_, nrow(dt))
  rk <- tryCatch(reverse_km_quantiles(fu, dt$PHENO), error=function(e) list(q=c(Q25=NA,Q50=NA,Q75=NA), max_followup=NA_real_))
  ev <- dt$PHENO == 1
  out <- data.table::data.table(
    summary = label,
    N = nrow(dt), N_events = sum(ev, na.rm=TRUE), N_censored = sum(dt$PHENO == 0, na.rm=TRUE),
    event_percent = round(100 * mean(dt$PHENO == 1, na.rm=TRUE), 2),
    reverseKM_median = rk$q["Q50"], reverseKM_q25 = rk$q["Q25"], reverseKM_q75 = rk$q["Q75"],
    reverseKM_max = rk$max_followup
  )
  if ("AGE_T1" %in% names(dt)) out[, `:=`(age_T1_median=q(dt$AGE_T1,.5), age_T1_q25=q(dt$AGE_T1,.25), age_T1_q75=q(dt$AGE_T1,.75))]
  if (all(c("AGE_T1","TIME_SINCE_T1") %in% names(dt))) {
    age_t2_ev <- (dt$AGE_T1 + dt$TIME_SINCE_T1)[ev]
    int_ev    <- dt$TIME_SINCE_T1[ev]
    out[, `:=`(age_T2_event_median=q(age_t2_ev,.5),
               t1_t2_event_median=q(int_ev,.5), t1_t2_event_q25=q(int_ev,.25), t1_t2_event_q75=q(int_ev,.75))]
  }
  if ("T1_STATUS" %in% names(dt)) {
    out[, N_prevalent_T1 := sum(dt$T1_STATUS=="Prevalent", na.rm=TRUE)]
    out[, N_incident_T1  := sum(dt$T1_STATUS=="Incident",  na.rm=TRUE)]
    out[, events_prevalent_T1 := sum(ev & dt$T1_STATUS=="Prevalent", na.rm=TRUE)]
    out[, events_incident_T1  := sum(ev & dt$T1_STATUS=="Incident",  na.rm=TRUE)]
    # Secondary-disclosure control: suppress BOTH sides of the prevalent/incident
    # partition when EITHER side is small — hiding only the small side leaves it
    # reconstructable from the retained overall N/events. (opt is a script global.)
    thr <- if (!is.null(opt$min_cell_count)) opt$min_cell_count else 10L
    split_cols <- c("N_prevalent_T1","N_incident_T1","events_prevalent_T1","events_incident_T1")
    small_split <- (min(out$N_prevalent_T1, out$N_incident_T1) < thr) ||
                   (min(out$events_prevalent_T1, out$events_incident_T1) < thr)
    out[, minimum_cell_count_pass_t1_split := !small_split]
    if (small_split) out[, (split_cols) := lapply(.SD, function(x) NA_integer_), .SDcols = split_cols]
  }
  out[]
}

revkm_rows <- list()
fu_vec <- prs$FOLLOWUP
rk_overall <- reverse_km_quantiles(fu_vec, prs$PHENO)
revkm_rows[[length(revkm_rows)+1]] <- data.table(
  Group = "Overall (Reverse KM on FOLLOWUP)",
  N = nrow(prs),
  Censors = sum(prs$PHENO == 0),
  Events = sum(prs$PHENO == 1),
  Median_FU = rk_overall$q["Q50"],
  Q25_FU = rk_overall$q["Q25"],
  Q75_FU = rk_overall$q["Q75"],
  Max_FU = rk_overall$max_followup
)
cat(sprintf("Reverse KM (overall): Median follow-up = %.2f (Q25 %.2f, Q75 %.2f), Max %.2f | N=%s, Censors=%s\n",
            revkm_rows[[1]]$Median_FU, revkm_rows[[1]]$Q25_FU, revkm_rows[[1]]$Q75_FU, revkm_rows[[1]]$Max_FU,
            format(nrow(prs), big.mark=","), format(sum(prs$PHENO==0), big.mark=",")))

# By sex
if (!is.null(opt$sex_col) && nzchar(opt$sex_col) && opt$sex_col %in% names(prs)) {
  sex_map <- parse_sex_labels(opt$sex_labels)
  sex_values <- sort(unique(prs[[opt$sex_col]])); sex_values <- sex_values[!is.na(sex_values)]
  sex_labels <- sapply(sex_values, function(s) {
    s_str <- as.character(s)
    if (!is.null(sex_map) && s_str %in% names(sex_map)) sex_map[[s_str]]
    else if (!is.na(suppressWarnings(as.numeric(s_str)))) {
      num_str <- as.character(as.numeric(s_str))
      if (!is.null(sex_map) && num_str %in% names(sex_map)) sex_map[[num_str]] else paste0("Sex=", s)
    } else paste0("Sex=", s)
  })
  for (i in seq_along(sex_values)) {
    sd <- prs[get(opt$sex_col) == sex_values[i]]
    if (nrow(sd) > 0) {
      rk <- reverse_km_quantiles(sd$FOLLOWUP, sd$PHENO)
      revkm_rows[[length(revkm_rows)+1]] <- data.table(
        Group = paste0(sex_labels[i], " (Reverse KM on FOLLOWUP)"),
        N = nrow(sd),
        Censors = sum(sd$PHENO == 0),
        Events = sum(sd$PHENO == 1),
        Median_FU = rk$q["Q50"],
        Q25_FU = rk$q["Q25"],
        Q75_FU = rk$q["Q75"],
        Max_FU = rk$max_followup
      )
    }
  }
  cat("\n")
}

reverse_km_table <- rbindlist(revkm_rows, use.names=TRUE, fill=TRUE)
num_cols_rkm <- c("Median_FU","Q25_FU","Q75_FU","Max_FU")
reverse_km_table[, (num_cols_rkm) := lapply(.SD, function(x) round(as.numeric(x), 2)), .SDcols=num_cols_rkm]
# Disclosure: suppress a small SUBGROUP row (per group; the large overall row
# passes). Blanks N/Events/Censors and the follow-up estimates for that row.
suppress_small_cells(reverse_km_table,
                     count_cols = c("N","Events","Censors"),
                     threshold = if (!is.null(opt$min_cell_count)) opt$min_cell_count else 10L,
                     estimate_cols = num_cols_rkm)
# Secondary disclosure across the sex rows: with the overall row (row 1) plus one
# visible sex row, a suppressed sex can be reconstructed by subtraction. If ANY
# per-sex (non-overall) row was suppressed, blank the counts + FU of ALL per-sex
# rows; the overall row is retained.
if (nrow(reverse_km_table) > 1) {
  sex_rows <- 2:nrow(reverse_km_table)
  if (any(!reverse_km_table$minimum_cell_count_pass[sex_rows])) {
    blank_cols <- intersect(c("N","Events","Censors", num_cols_rkm), names(reverse_km_table))
    reverse_km_table[sex_rows, (blank_cols) := lapply(.SD, function(x) NA), .SDcols = blank_cols]
    reverse_km_table[sex_rows, minimum_cell_count_pass := FALSE]
  }
}
reverse_km_file <- file.path(opt$outdir, paste0(prefix, "followup_reverseKM.csv"))
fwrite(reverse_km_table, reverse_km_file)
cat("✓ Saved reverse KM follow-up summary:", basename(reverse_km_file), "\n\n")

# ==============================================================================
# Sample size guard
# ==============================================================================
n_events_total <- sum(prs$PHENO == 1)
# Initialise n_events_final here so single-PRS runs (which skip the multi-PRS
# alignment block that normally sets it) do not hit an undefined variable later.
# Overwritten with the aligned-sample count after three-PRS alignment.
n_events_final <- n_events_total
if (length(unique(prs$PHENO)) < 2)
  stop("Event status has no variation after QC (all ", unique(prs$PHENO)[1],
       "); nothing is estimable.")
# Hard estimability floors, deliberately INDEPENDENT of --min_events_total (which
# is a scientific reliability threshold checked just below, and can be lowered).
# These fire first, so they must say what happened — a bare "Insufficient data
# after QC" gives no way to tell a small cohort from a QC step that ate the sample.
if (nrow(prs) < 50 || n_events_total < 10) {
  # trajectory_attrition.csv is only written at the very end of a successful run,
  # so it does not exist yet — render the cohort flow inline instead. Without it
  # there is no way to tell a genuinely small cohort from a QC step that ate the
  # sample (e.g. partial recruitment-age coverage in prospective mode).
  flow <- tryCatch(paste(utils::capture.output(print(data.table::rbindlist(attrition_log))),
                         collapse = "\n"),
                   error = function(e) "  (attrition log unavailable)")
  stop(sprintf(paste0(
    "Insufficient data after QC: N = %d, T2 events = %d.\n",
    "  Hard floors are N >= 50 and events >= 10. They are NOT controlled by\n",
    "  --min_events_total (currently %d), so lowering that flag will not bypass them.\n",
    "  Cohort flow up to this point:\n%s"),
    nrow(prs), n_events_total, opt$min_events_total, flow))
}
# Cheap PRE-check only. The authoritative min-events guard runs AFTER three-PRS
# alignment (the aligned event count is what governs estimability), and the
# events-per-parameter warning uses the ACTUAL fitted parameter count (see below,
# after model fitting) rather than a fixed "95 = EPV 5" rule, because the true
# parameter count depends on spline df, factor levels, batch, strata and PCs.
if (n_events_total < opt$min_events_total) {
  stop(sprintf(paste0("Only %d T2 events after QC (< --min_events_total = %d).\n",
                      "  Full model (M7) not reliably estimable; this is re-checked on the aligned sample.\n",
                      "  Lower --min_events_total only if you accept the instability."),
              n_events_total, opt$min_events_total))
}

# ==============================================================================
# PRS direction & standardisation
# ==============================================================================
# Direction must come from discovery-stage allele harmonisation, NOT from whether
# the score predicts the outcome positively in this cohort. Default: no flipping.
#   --prs_signs   : explicit prespecified per-PRS multipliers +1/-1 (handled in the
#                   multi-PRS block below; -1 flips that score, order onset,progression,outcome)
#   --auto_flip_prs (default FALSE, DISCOURAGED): outcome-dependent auto-flip
prs_pheno_cor_raw <- cor(prs$SCORESUM, prs$PHENO, use="complete.obs")
prs_pheno_cor <- prs_pheno_cor_raw
prs_direction <- "as_supplied"

if (isTRUE(opt$auto_flip_prs) && prs_pheno_cor < 0) {
  warning("--auto_flip_prs is target-outcome dependent (leakage). ",
          "Direction should come from discovery-stage harmonisation. Flipping anyway as requested.")
  cat("⚠ WARNING: primary PRS auto-flipped on OUTCOME correlation (--auto_flip_prs) — target-outcome dependent\n")
  prs[, SCORESUM := -SCORESUM]
  prs_pheno_cor <- cor(prs$SCORESUM, prs$PHENO)
  prs_direction <- "auto(outcome)"
}

# Standardisation: whole-sample by default (outcome-independent, comparable HR/SD).
# --zscore_by_controls uses control mean/SD; an affine rescale that changes only the
# HR-per-SD scale, not the LRT/AUC, but is weakly outcome-dependent — discouraged.
if (isTRUE(opt$zscore_by_controls)) {
  warning("--zscore_by_controls uses outcome (control) status to set the PRS scale; ",
          "this changes only the HR-per-SD scale, not the tests, but harms cross-cohort comparability.")
  cat("⚠ WARNING: standardising PRS by CONTROL mean/SD (--zscore_by_controls) — whole-sample scaling is recommended\n")
  ctrl_mean <- mean(prs[PHENO==0, SCORESUM]); ctrl_sd <- sd(prs[PHENO==0, SCORESUM])
  prs[, ZSCORE := (SCORESUM - ctrl_mean) / ctrl_sd]
} else {
  prs[, ZSCORE := scale(SCORESUM)[,1]]   # whole-sample; re-standardised after 3-PRS alignment
}

# ==============================================================================
# MODEL DATA (DT) — TIME is FOLLOWUP for evaluation metrics
# ==============================================================================
DT <- prs[, .(IID, PRS_original=SCORESUM, PRS_std=ZSCORE, PHENO, START, STOP, TIME=FOLLOWUP)]
if (length(mycovs) > 0) DT <- cbind(DT, prs[, mycovs, with=FALSE])
if (length(clinical_covs) > 0) DT <- cbind(DT, prs[, clinical_covs, with=FALSE])
# Carry age-at-T1 (spline covariate + index age) and, if present, timing category.
if (AGE_T1_AVAILABLE) DT[, AGE_T1 := prs$AGE_T1]
# Carry TIME_SINCE_T1 (the T1->T2 interval) so the post-alignment summary can
# compute age-at-T2 and the interval even in prospective mode, where DT$TIME is
# time-from-index rather than time-since-T1.
if ("TIME_SINCE_T1" %in% names(prs)) DT[, TIME_SINCE_T1 := prs$TIME_SINCE_T1]
if (!is.null(prs$timing_category)) DT[, timing_category := factor(prs$timing_category)]
# Prospective-mode covariates (index age, disease duration at index, timing status).
is_prospective <- identical(opt$analysis_mode, "prospective")
if (is_prospective) {
  DT[, AGE_INDEX   := prs$AGE_INDEX]
  DT[, T1_DURATION := prs$T1_DURATION]
  DT[, T1_STATUS   := factor(prs$T1_STATUS)]
}

# ------------------------------------------------------------------------------
# COVARIATE TERMS (formula-only; kept out of the plain-column covariate vector so
# column-subsetting/VIF code is unaffected). Two buckets:
#   extra_terms : spline covariates, SAFE for both coxph and the VIF lm().
#   strata_term : strata(), coxph formulas ONLY (invalid in lm()).
# Mode-specific:
#   gwas_aligned : ns(AGE_T1)                              [age at T1]
#   prospective  : ns(AGE_INDEX) + ns(T1_DURATION) + strata(T1_STATUS)
#                  (drop AGE_T1: AGE_INDEX = AGE_T1 + T1_DURATION is collinear)
# ------------------------------------------------------------------------------
extra_terms <- character(0)
strata_term <- character(0)
spline_or_linear <- function(col, label) {
  x  <- DT[[col]][is.finite(DT[[col]])]
  nd <- length(unique(x))
  if (nd < 2) {
    # Constant covariate: unestimable. Entering it yields an NA coefficient that
    # later breaks riskRegression's predictCox ("parameters ... have no value 'NA'").
    # This happens for T1_DURATION when the prospective sample is all-incident (e.g. a
    # register/birth-start cohort indexed at birth -> everyone incident -> T1_DURATION=0).
    cat("  ⚠ Covariate ", col, " is constant (", nd, " distinct value); dropped from the model (",
        label, ")\n", sep="")
    return(character(0))
  }
  # Use a spline only if ns() can actually build its basis on this covariate. >=4 distinct
  # values is necessary but NOT sufficient: a heavily tied / zero-inflated covariate collapses
  # ns()'s quantile knots onto the boundary -> "all interior knots match left boundary knot".
  # This bites T1_DURATION whenever most of the sample is incident (T1_DURATION = 0 for every
  # incident person, so a majority-incident prospective cohort is >2/3 zeros). Building the
  # basis here mirrors exactly what coxph will do; fall back ONLY on a genuine ns() error, so a
  # covariate whose spline already worked keeps the identical term (finished runs unchanged).
  spline_ok <- FALSE
  if (opt$age_spline_df > 0 && nd >= 4) {
    spline_ok <- tryCatch({ splines::ns(x, df = opt$age_spline_df); TRUE },
                          error = function(e) FALSE)
  }
  if (spline_ok) {
    t <- sprintf("ns(%s, df=%d)", col, opt$age_spline_df)
    cat("  Covariate: ", t, " (", label, ")\n", sep=""); t
  } else {
    reason <- if (opt$age_spline_df <= 0) "--age_spline_df 0"
              else if (nd < 4) sprintf("only %d distinct", nd)
              else "spline knots collapse on a heavily tied / zero-inflated covariate"
    cat("  Covariate: ", col, " entered linearly (", reason, ")\n", sep=""); col
  }
}
if (is_prospective) {
  if ("AGE_INDEX"   %in% names(DT)) extra_terms <- c(extra_terms, spline_or_linear("AGE_INDEX", "age at prediction index"))
  if ("T1_DURATION" %in% names(DT)) extra_terms <- c(extra_terms, spline_or_linear("T1_DURATION", "T1 duration at index"))
  if ("T1_STATUS" %in% names(DT) && nlevels(DT$T1_STATUS) >= 2) {
    strata_term <- "strata(T1_STATUS)"
    cat("  Stratifying baseline hazard by T1_STATUS (prevalent vs incident)\n")
  } else if ("T1_STATUS" %in% names(DT)) {
    cat("  ⚠ T1_STATUS has <2 levels; strata term dropped\n")
  }
} else if (AGE_T1_AVAILABLE) {
  extra_terms <- spline_or_linear("AGE_T1", "age at T1")
}

cat("✓ Data preparation completed | N:", nrow(DT), " Events:", sum(DT$PHENO), " | Start–stop: TRUE\n\n")

# ==============================================================================
# MULTI-PRS: Load additional PRS files and add standardized columns to DT
# ==============================================================================
multi_prs_active <- !is.null(opt$prs_labels) && nzchar(opt$prs_labels)
multi_prs_std_cols <- list()
prs_role_cols      <- list()   # role ("onset"/"progression"/"outcome") -> std column
model_specs <- NULL

if (multi_prs_active) {
  cat("==============================================================================\nMULTI-PRS MODE\n==============================================================================\n\n")
  # Outcome-dependent PRS options are not permitted in the three-PRS manuscript
  # path: control-based scaling is overwritten by the aligned-sample re-scaling
  # (so metadata would misreport it), and outcome-based auto-flipping is leakage.
  if (isTRUE(opt$zscore_by_controls))
    stop("--zscore_by_controls is not supported in the three-PRS manuscript analysis ",
         "(whole-sample scaling is used; control-based scaling is outcome-dependent).")
  if (isTRUE(opt$auto_flip_prs))
    stop("--auto_flip_prs (outcome-based) is not permitted in the three-PRS manuscript analysis. ",
         "Use discovery-stage harmonisation or explicit --prs_signs.")
  all_prs_labels <- trimws(strsplit(opt$prs_labels, ",")[[1]])
  # Normalise roles to lower case ONCE, so every later comparison (role_of_label,
  # model_specs auto-build) is case-insensitive and consistent.
  all_prs_types  <- if (!is.null(opt$prs_types) && nzchar(opt$prs_types)) tolower(trimws(strsplit(opt$prs_types, ",")[[1]])) else rep("unknown", length(all_prs_labels))
  all_prs_cols   <- if (!is.null(opt$prs_col_list) && nzchar(opt$prs_col_list)) trimws(strsplit(opt$prs_col_list, ",")[[1]]) else rep(opt$prs_col, length(all_prs_labels))
  extra_files    <- if (!is.null(opt$prs_files) && nzchar(opt$prs_files)) trimws(strsplit(opt$prs_files, ",")[[1]]) else character(0)
  all_prs_files  <- c(opt$prs_file, extra_files)
  if (length(all_prs_files) != length(all_prs_labels))
    stop("Number of PRS files (", length(all_prs_files), ") must equal --prs_labels count (", length(all_prs_labels), ")")
  # Require the full three-PRS set, EXACTLY ONE score per role, for the M-ladder.
  needed_roles  <- c("onset", "outcome", "progression")
  role_counts   <- table(factor(all_prs_types, levels=needed_roles))
  missing_roles <- needed_roles[role_counts == 0]
  dup_roles     <- needed_roles[role_counts > 1]
  if (length(missing_roles) > 0)
    stop("The manuscript analysis requires all three PRS roles {onset, outcome, progression}.\n",
         "  Missing: ", paste(missing_roles, collapse=", "),
         "\n  Supplied --prs_types: ", paste(all_prs_types, collapse=", "))
  if (length(dup_roles) > 0)
    stop("Exactly one PRS is required per role; duplicated role(s): ", paste(dup_roles, collapse=", "),
         "\n  Duplicate scores would be mislabelled across M1-M7. Supplied --prs_types: ",
         paste(all_prs_types, collapse=", "))
  # Exactly three scores — a 4th with an unrecognised role would otherwise load
  # and enter complete-case alignment even though it is not part of M0-M7.
  if (length(all_prs_types) != 3L || !setequal(all_prs_types, needed_roles))
    stop("Exactly three PRSs are required (onset, outcome, progression). Supplied --prs_types: ",
         paste(all_prs_types, collapse=", "))
  # Parse explicit prespecified direction multipliers (--prs_signs), if given.
  prs_signs_vec <- NULL
  if (!is.null(opt$prs_signs) && nzchar(opt$prs_signs)) {
    prs_signs_vec <- suppressWarnings(as.numeric(trimws(strsplit(opt$prs_signs, ",")[[1]])))
    if (length(prs_signs_vec) != length(all_prs_labels) || any(!prs_signs_vec %in% c(-1, 1)))
      stop("--prs_signs must be ", length(all_prs_labels), " values of +1/-1 (one per --prs_labels).")
    names(prs_signs_vec) <- all_prs_labels
  }
  for (i in seq_along(all_prs_labels)) {
    lbl <- all_prs_labels[i]; std_col <- paste0(lbl, "_std")
    if (i == 1) {
      DT[[std_col]] <- DT$PRS_std
      cat("  ✓ Primary PRS:", lbl, "-> mapped from PRS_std\n")
    } else {
      if (!file.exists(all_prs_files[i])) stop("PRS file not found: ", all_prs_files[i])
      pf_dt <- fread(all_prs_files[i])
      if (!"IID" %in% names(pf_dt)) setnames(pf_dt, names(pf_dt)[1], "IID")
      assert_unique_iid(pf_dt, "IID", paste("PRS file", basename(all_prs_files[i])))
      sc <- if (i <= length(all_prs_cols)) all_prs_cols[i] else opt$prs_col
      if (!sc %in% names(pf_dt)) {
        cand <- grep("SCORE.*SUM|SCORE1_SUM|SCORESUM", names(pf_dt), value=TRUE, ignore.case=TRUE)
        if (!length(cand)) stop("Cannot find score column '", sc, "' in ", all_prs_files[i])
        sc <- cand[1]
      }
      raw_vals <- pf_dt[match(DT$IID, pf_dt$IID), get(sc)]
      if (opt$zscore_by_controls) {
        ctrl_mean <- mean(raw_vals[DT$PHENO==0], na.rm=TRUE)
        ctrl_sd   <- sd(raw_vals[DT$PHENO==0],   na.rm=TRUE)
        DT[[std_col]] <- (raw_vals - ctrl_mean) / ctrl_sd
      } else {
        DT[[std_col]] <- as.numeric(scale(raw_vals))
      }
      cat("  ✓ Loaded:", lbl, "from", basename(all_prs_files[i]), "->", std_col, "\n")
    }
    multi_prs_std_cols[[lbl]] <- std_col
    # Also key the standardized column by ROLE, so downstream display code can ask
    # for "the progression score" without knowing the biobank's label.
    if (i <= length(all_prs_types)) prs_role_cols[[all_prs_types[i]]] <- std_col
  }

  # ---------------------------------------------------------------------------
  # Direction of the ADDITIONAL PRS — DISCOURAGED outcome-based auto-flip.
  # Off by default; direction should come from discovery-stage harmonisation
  # or explicit --prs_signs (applied after alignment, below).
  # ---------------------------------------------------------------------------
  extra_prs_flipped <- character(0)
  if (isTRUE(opt$auto_flip_prs)) {
    warning("--auto_flip_prs also flips additional PRS on outcome correlation (leakage).")
    for (lbl in names(multi_prs_std_cols)) {
      std_col <- multi_prs_std_cols[[lbl]]
      if (identical(std_col, "PRS_std") || lbl == all_prs_labels[1]) next  # primary already handled
      cc <- suppressWarnings(cor(DT[[std_col]], DT$PHENO, use="complete.obs"))
      if (is.finite(cc) && cc < 0) {
        DT[[std_col]] <- -DT[[std_col]]
        extra_prs_flipped <- c(extra_prs_flipped, lbl)
        cat("  ⚠ Direction auto-flipped on OUTCOME correlation for:", lbl, "\n")
      }
    }
  }

  # ---------------------------------------------------------------------------
  # COMPLETE-CASE ALIGNMENT ACROSS ALL PRS  (required for nested model tests)
  # ---------------------------------------------------------------------------
  # Extra PRS are joined with match(), so IIDs absent from a .sscore file yield
  # NA. The complete-case filter earlier in this script runs BEFORE these
  # columns exist, so models containing an extra PRS would silently drop rows
  # that the Base model keeps. anova.coxph then refuses outright:
  #   "models were not all fit to the same size of dataset"
  # Every downstream comparison (LRT, delta C-index, Score() contrasts) assumes
  # a common sample, so align here.
  std_cols_all <- unlist(multi_prs_std_cols, use.names = FALSE)
  std_cols_all <- intersect(std_cols_all, names(DT))
  # Cohort summary on the sample entering three-PRS alignment (pre-alignment).
  summary_pre <- summarise_cohort(DT, "pre_alignment")
  if (length(std_cols_all)) {
    n_before <- nrow(DT)
    keep <- Reduce(`&`, lapply(std_cols_all, function(cc) is.finite(DT[[cc]])))
    DT <- DT[keep]
    record_stage(DT, "complete_three_PRS")
    n_dropped <- n_before - nrow(DT)
    if (n_dropped > 0) {
      cat("\n  ⚠ Dropped", format(n_dropped, big.mark=","),
          "rows with missing values in one or more PRS (",
          sprintf("%.2f%%", 100*n_dropped/n_before), ")\n", sep=" ")
      cat("    All models are now fitted on an identical set of",
          format(nrow(DT), big.mark=","), "individuals.\n")
    } else {
      cat("\n  ✓ All", format(nrow(DT), big.mark=","),
          "individuals have complete data for every PRS\n")
    }
    if (nrow(DT) < 50 || sum(DT$PHENO == 1) < 10)
      stop("Insufficient data after aligning PRS across models (N=", nrow(DT),
           ", events=", sum(DT$PHENO == 1), ")")
    multi_prs_n_dropped <- n_dropped

    # ------------------------------------------------------------------------
    # RE-STANDARDISE on the aligned sample, then apply prespecified signs.
    # Standardising BEFORE the three-PRS complete-case restriction leaves the
    # scores no longer exactly mean 0 / SD 1 in the final sample, weakening
    # HR-per-SD comparability. Canonical order: scale -> sign -> model.
    # This runs before quantile creation and fitting, so quantiles and the
    # Minimal model use the final scores.
    # ------------------------------------------------------------------------
    prs_signs_applied <- character(0)
    for (lbl in names(multi_prs_std_cols)) {
      cc <- multi_prs_std_cols[[lbl]]
      if (!cc %in% names(DT)) next
      DT[[cc]] <- as.numeric(scale(DT[[cc]]))
      if (!is.null(prs_signs_vec) && lbl %in% names(prs_signs_vec) && prs_signs_vec[[lbl]] == -1) {
        DT[[cc]] <- -DT[[cc]]
        prs_signs_applied <- c(prs_signs_applied, lbl)
      }
    }
    # Keep the primary alias PRS_std synced to the (scaled, signed) onset column.
    primary_std <- multi_prs_std_cols[[all_prs_labels[1]]]
    if (!is.null(primary_std) && primary_std %in% names(DT)) DT[, PRS_std := get(primary_std)]
    cat("  ✓ Re-standardised all PRS on the aligned sample (mean 0, SD 1)\n")
    if (length(prs_signs_applied))
      cat("  ✓ Applied prespecified sign flip (--prs_signs) to:", paste(prs_signs_applied, collapse=", "), "\n")

    # PRS orientation is checked AFTER model fitting, off each single-PRS model's
    # adjusted coefficient (the published estimand), not the marginal correlation here.

    # ------------------------------------------------------------------------
    # RE-CHECK the minimum-event guard on the ALIGNED sample (the count that
    # matters for model estimability). The early guard used the pre-alignment
    # count and could pass a stratum that falls below the floor after alignment.
    # ------------------------------------------------------------------------
    n_events_final <- sum(DT$PHENO == 1)
    if (n_events_final < opt$min_events_total)
      stop(sprintf(paste0("Only %d T2 events after three-PRS alignment (< --min_events_total = %d).\n",
                          "  The full model (M7) is not reliably estimable on the aligned sample.\n",
                          "  Lower --min_events_total only if you accept the instability."),
                  n_events_final, opt$min_events_total))
  } else {
    multi_prs_n_dropped <- 0L
    prs_signs_applied <- character(0)
    n_events_final <- sum(DT$PHENO == 1)
  }

  # Recheck prospective T1_STATUS on the ALIGNED sample: a status could vanish
  # after restricting to three-PRS-complete participants. Drop empty levels and,
  # if only one status remains, drop the strata term (strata() on a single level
  # would error). Runs before model_covs is built.
  if (exists("is_prospective") && isTRUE(is_prospective) && "T1_STATUS" %in% names(DT)) {
    DT[, T1_STATUS := droplevels(T1_STATUS)]
    if (data.table::uniqueN(DT$T1_STATUS) < 2 && length(strata_term)) {
      strata_term <- character(0)
      cat("  ⚠ Only one T1_STATUS remains after PRS alignment; strata(T1_STATUS) removed.\n")
    }
    t1_tab <- DT[, .(N = .N, Events = sum(PHENO)), by = T1_STATUS][order(T1_STATUS)]
    # Disclosure: blank small status counts (the strata-collapse guard above keys
    # off DT$T1_STATUS directly, not this written table).
    suppress_partition_all_or_none(t1_tab, count_cols = c("N","Events"), threshold = opt$min_cell_count)
    fwrite(t1_tab, file.path(opt$outdir, paste0(prefix, "t1_status_final.csv")))
  }
  cat("\n")

  # Build model_specs
  if (!is.null(opt$models) && nzchar(opt$models)) {
    model_specs <- list()
    for (spec in trimws(strsplit(opt$models, ";")[[1]])) {
      parts <- strsplit(spec, "=", fixed=TRUE)[[1]]; mn <- trimws(parts[1])
      lbls  <- if (length(parts)>1 && nzchar(trimws(parts[2]))) trimws(strsplit(trimws(parts[2]), "+", fixed=TRUE)[[1]]) else character(0)
      model_specs[[mn]] <- lbls
    }
  } else {
    onset_lbls   <- all_prs_labels[all_prs_types=="onset"]
    prog_lbls    <- all_prs_labels[all_prs_types=="progression"]
    outcome_lbls <- all_prs_labels[all_prs_types=="outcome"]
    model_specs  <- list()
    if (length(onset_lbls))   model_specs[["Onset"]]       <- onset_lbls
    if (length(prog_lbls))    model_specs[["Progression"]] <- prog_lbls
    if (length(outcome_lbls)) model_specs[["Outcome"]]     <- outcome_lbls
    if (length(onset_lbls)   && length(prog_lbls))    model_specs[["Onset_Progression"]]       <- c(onset_lbls, prog_lbls)
    if (length(outcome_lbls) && length(prog_lbls))    model_specs[["Outcome_Progression"]]     <- c(outcome_lbls, prog_lbls)
    if (length(onset_lbls)   && length(outcome_lbls)) model_specs[["Onset_Outcome"]]           <- c(onset_lbls, outcome_lbls)
    if (length(onset_lbls)   && length(prog_lbls) && length(outcome_lbls))
      model_specs[["Onset_Progression_Outcome"]] <- c(onset_lbls, prog_lbls, outcome_lbls)
  }
  cat("\nModel combinations to fit:\n")
  for (nm in names(model_specs)) cat("  ", nm, ":", paste(model_specs[[nm]], collapse=", "), "\n")
  cat("\n")

  # ---------------------------------------------------------------------------
  # MANUSCRIPT MODEL MAP  (M0-M7)  -- single source of truth
  # ---------------------------------------------------------------------------
  # Downstream scripts join this rather than re-deriving the mapping, so the
  # manuscript label for a model is defined in exactly one place.
  #
  #   {}                              -> M0   covariates only
  #   {onset}                         -> M1   + PRS_T0T1
  #   {outcome}                       -> M2   + PRS_T0T2
  #   {progression}                   -> M3   + PRS_T1T2
  #   {onset, progression}            -> M4
  #   {outcome, progression}          -> M5
  #   {onset, outcome}                -> M6   + PRS_T0T1 + PRS_T0T2  (two-PRS)
  #   {onset, outcome, progression}   -> M7   all three  (full model)
  #
  # M6 is a two-PRS model; M7 is the full model. This matches the manuscript's
  # final table. Leave-one-out tests are taken from M7.
  role_of_label <- setNames(all_prs_types, all_prs_labels)
  prs_id_of_role <- c(onset = "PRS_T0T1", outcome = "PRS_T0T2", progression = "PRS_T1T2")

  model_id_for_roles <- function(roles) {
    key <- paste(sort(unique(roles)), collapse = "+")
    if (!nzchar(key)) return("M0")   # covariate-only; switch() cannot take "" as a case name
    switch(key,
           "onset"                       = "M1",
           "outcome"                     = "M2",
           "progression"                 = "M3",
           "onset+progression"           = "M4",
           "outcome+progression"         = "M5",
           "onset+outcome"               = "M6",
           "onset+outcome+progression"   = "M7",
           NA_character_)
  }

  map_rows <- list(data.table(
    model = "Base", model_id = "M0", prs_terms = "", n_prs = 0L, ladder = "core"
  ))
  for (nm in names(model_specs)) {
    lbls  <- model_specs[[nm]]
    roles <- unname(role_of_label[lbls])
    roles <- roles[!is.na(roles)]
    ids   <- unname(prs_id_of_role[roles]); ids <- ids[!is.na(ids)]
    map_rows[[length(map_rows) + 1L]] <- data.table(
      model     = nm,
      model_id  = model_id_for_roles(roles),
      prs_terms = paste(sort(ids), collapse = "+"),
      n_prs     = length(ids),
      ladder    = "core"
    )
  }
  model_map <- rbindlist(map_rows, fill = TRUE)
  # Mirror clinical ladder rows (Clinical_<model>) when clinical covariates exist.
  if (length(clinical_covs) > 0) {
    clin_map <- copy(model_map)
    clin_map[, `:=`(model = paste0("Clinical_", model), ladder = "clinical")]
    model_map <- rbindlist(list(model_map, clin_map), use.names = TRUE)
  }
  setorder(model_map, ladder, n_prs, model_id, na.last = TRUE)

  cat("Manuscript model map:\n")
  for (i in seq_len(nrow(model_map))) {
    cat(sprintf("  %-28s %-4s %s\n",
                model_map$model[i],
                ifelse(is.na(model_map$model_id[i]), "-", model_map$model_id[i]),
                ifelse(nzchar(model_map$prs_terms[i]), model_map$prs_terms[i], "covariates only")))
  }
  unmapped <- model_map[is.na(model_id), model]
  if (length(unmapped))
    cat("  Note: no manuscript id for:", paste(unmapped, collapse=", "),
        "(reported with model_id = NA)\n")
  cat("\n")
}

# ==============================================================================
# TIME SUPPORT (analysis follow-up scale)
# ==============================================================================
cat("==============================================================================\nTIME SUPPORT (analysis follow-up scale)\n==============================================================================\n\n")
min_events_val     <- opt$min_events
min_at_risk_val    <- opt$min_at_risk
ghat_threshold_val <- opt$ghat_threshold

supported_times <- compute_supported_times_from_data(
  DT, time_col="TIME", status_col="PHENO",
  min_events=min_events_val, min_at_risk=min_at_risk_val,
  ghat_threshold=ghat_threshold_val, cap_quantile=0.95
)
ev_times <- sort(unique(DT$TIME[DT$PHENO == 1]))
extra_grid <- if (length(ev_times) > 0) stats::quantile(ev_times, probs = seq(0.1, 0.9, by = 0.1), na.rm = TRUE) else numeric(0)
candidate_times <- sort(unique(c(supported_times, extra_grid)))
time_support_table <- build_time_support_table(
  DT, candidate_times,
  time_col="TIME", status_col="PHENO",
  min_events=min_events_val, min_at_risk=min_at_risk_val, ghat_threshold=ghat_threshold_val
)
if (length(supported_times)) cat("Auto-supported time points (years): ", paste(round(supported_times, 2), collapse=", "), "\n", sep="") else cat("No auto-supported time points passed guards.\n")
time_support_file <- file.path(opt$outdir, paste0(prefix, "time_support.csv"))
fwrite(time_support_table, time_support_file)
cat("✓ Saved auto time-support report:", basename(time_support_file), "\n\n")

# Requested points
requested_time_points_raw <- trimws(opt$time_points)
requested_time_points <- numeric(0)
if (nchar(requested_time_points_raw) > 0) {
  req_parts <- trimws(strsplit(requested_time_points_raw, ",")[[1]])
  requested_time_points <- suppressWarnings(as.numeric(req_parts))
  requested_time_points <- requested_time_points[is.finite(requested_time_points)]
  requested_time_points <- sort(unique(requested_time_points))
}
requested_ghat_table <- data.table(); requested_valid_times <- numeric(0)
if (length(requested_time_points) > 0) {
  cat("Requested time points (analysis scale):", paste(requested_time_points, collapse=", "), " (years)\n\n")
  req_audit <- audit_requested_time_points(
    DT, requested_time_points, time_col="TIME", status_col="PHENO",
    min_events=min_events_val, min_at_risk=min_at_risk_val, ghat_threshold=ghat_threshold_val
  )
  requested_ghat_table <- req_audit$ghat_table
  requested_valid_times <- req_audit$valid_times
  req_audit_csv <- file.path(opt$outdir, paste0(prefix, "requested_timepoints_audit.csv"))
  fwrite(requested_ghat_table, req_audit_csv)
  cat("✓ Saved requested time-points AUDIT:", basename(req_audit_csv), "\n")
  req_pass_csv <- file.path(opt$outdir, paste0(prefix, "requested_timepoints_ghat.csv"))
  if (length(requested_valid_times)) {
    fwrite(data.table(time_years = requested_valid_times), req_pass_csv)
    cat("✓ Saved requested PASS-only:", basename(req_pass_csv), "\n\n")
  } else {
    fwrite(data.table(time_years = numeric(0)), req_pass_csv)
    cat("⚠ None of the requested time points passed; PASS-only file is empty:", basename(req_pass_csv), "\n\n")
  }
} else {
  cat("No --time_points provided; only auto-supported grid recorded.\n")
  req_pass_csv <- file.path(opt$outdir, paste0(prefix, "requested_timepoints_ghat.csv"))
  if (length(supported_times)) {
    fwrite(data.table(time_years = supported_times), req_pass_csv)
    cat("✓ Saved auto-supported PASS-only time points:", basename(req_pass_csv), "\n\n")
  } else {
    fwrite(data.table(time_years = numeric(0)), req_pass_csv)
    cat("⚠ No supported times; PASS-only file is empty:", basename(req_pass_csv), "\n\n")
  }
}
eval_time_points <- if (length(requested_valid_times) > 0) requested_valid_times else supported_times
if (length(eval_time_points) == 0) cat("⚠ eval_time_points is empty (no time passed guards). Downstream scripts must handle this.\n\n")
cat("Eval time points to use downstream (analysis scale):", ifelse(length(eval_time_points)>0, paste(eval_time_points, collapse=", "), "NONE"), "\n\n")

# ==============================================================================
# STRATIFICATION (create/use)
# ==============================================================================
if (!is.null(opt$strata_var) && !is.null(opt$strata_breaks)) {
  cat("==============================================================================\nCREATING STRATIFICATION VARIABLE\n==============================================================================\n\n")
  strata_var <- opt$strata_var
  if (!strata_var %in% names(DT)) stop("Stratification variable '", strata_var, "' not found in data.")
  breaks_str <- trimws(strsplit(opt$strata_breaks, ",")[[1]])
  breaks <- sapply(breaks_str, function(x) if (tolower(x)=="inf") Inf else if (tolower(x)=="-inf") -Inf else as.numeric(x))
  if (any(is.na(breaks[!is.infinite(breaks)]))) stop("Invalid breakpoints: ", opt$strata_breaks)
  if (!is.null(opt$strata_labels) && opt$strata_labels != "") {
    labels <- trimws(strsplit(opt$strata_labels, ",")[[1]])
    if (length(labels) != length(breaks)-1) stop("Number of labels must equal number of intervals")
  } else labels <- paste0("Stratum_", 1:(length(breaks)-1))
  strata_col <- paste0(strata_var, "_strata")
  DT[[strata_col]] <- cut(DT[[strata_var]], breaks=breaks, labels=labels, include_lowest=TRUE, right=FALSE)
  n_na <- sum(is.na(DT[[strata_col]])); if (n_na>0) { cat("⚠", n_na, "observations outside breakpoint range -> excluded\n"); DT <- DT[!is.na(get(strata_col))] }
  if (strata_var %in% mycovs) mycovs <- mycovs[mycovs != strata_var]
  if (strata_var %in% clinical_covs) clinical_covs <- clinical_covs[clinical_covs != strata_var]
  opt$strata <- strata_col
  cat("✓ Stratification variable created:", strata_col, "\n\n")
} else if (!is.null(opt$strata)) {
  cat("==============================================================================\nUSING PRE-EXISTING STRATIFICATION VARIABLE\n==============================================================================\n\n")
  if (!opt$strata %in% names(DT)) stop("Stratification variable '", opt$strata, "' not found in data.")
  if (opt$strata %in% mycovs)       mycovs       <- mycovs[mycovs       != opt$strata]
  if (opt$strata %in% clinical_covs) clinical_covs <- clinical_covs[clinical_covs != opt$strata]
  cat("Stratification variable:", opt$strata, "\n\n")
}

# ==============================================================================
# PRS QUANTILES
# ==============================================================================
cat("==============================================================================\nCREATING PRS QUANTILES\n==============================================================================\n\n")
create_standard_quantiles <- function(prs_std, n_quantiles, ref_group="lowest") {
  breaks <- quantile(prs_std, probs=seq(0, 1, length.out=n_quantiles+1), na.rm=TRUE)
  labels <- paste0("Q", 1:n_quantiles)
  quant <- cut(prs_std, breaks=breaks, labels=labels, include.lowest=TRUE)
  if (ref_group == "middle") {
    if (n_quantiles %% 2 == 0) {
      middle_low <- n_quantiles/2; middle_high <- middle_low + 1
      ref_groups <- paste0("Q", c(middle_low, middle_high))
      quant <- as.character(quant); quant[quant %in% ref_groups] <- "Q_ref"
      all_labels <- c(paste0("Q", 1:(middle_low-1)), "Q_ref", paste0("Q", (middle_high+1):n_quantiles))
      all_labels <- all_labels[nchar(all_labels) > 0]
      quant <- factor(quant, levels=all_labels); ref_q <- "Q_ref"
    } else ref_q <- paste0("Q", ceiling(n_quantiles/2))
  } else if (ref_group == "lowest") ref_q <- "Q1" else ref_q <- ref_group
  quant <- relevel(quant, ref=ref_q)
  list(quantile=quant, numeric=as.numeric(quant), n_quantiles=n_quantiles, reference=ref_q, breaks=breaks)
}
create_extreme_quantile_single <- function(prs_std, percentile) {
  groups <- rep("Middle", length(prs_std))
  top_thresh <- quantile(prs_std, 1 - percentile/100, na.rm=TRUE)
  groups[prs_std >= top_thresh] <- paste0("Top", percentile, "%")
  bottom_thresh <- quantile(prs_std, percentile/100, na.rm=TRUE)
  groups[prs_std <= bottom_thresh] <- paste0("Bottom", percentile, "%")
  all_labels <- c(paste0("Bottom", percentile, "%"), "Middle", paste0("Top", percentile, "%"))
  groups <- factor(groups, levels=all_labels); groups <- relevel(groups, ref="Middle")
  list(quantile=groups, percentile=percentile, reference="Middle", n_groups=3)
}
quantile_specs <- trimws(strsplit(opt$prs_quantiles, ",")[[1]])
quantile_list <- list()
if (length(quantile_specs) > 0) {
  cat("1. Standard quantiles:\n")
  for (spec in quantile_specs) {
    n_quant <- suppressWarnings(as.integer(spec)); if (is.na(n_quant) || n_quant < 2) { warning("Invalid quantile spec: ", spec); next }
    quant_name <- paste0("Q", n_quant); col_name <- paste0("PRS_", quant_name)
    cat("   - ", quant_name, "\n", sep="")
    qq <- create_standard_quantiles(DT$PRS_std, n_quant, opt$prs_quantile_reference)
    DT[[col_name]] <- qq$quantile; DT[[paste0(col_name, "_numeric")]] <- qq$numeric
    quantile_list[[quant_name]] <- qq
    print(table(qq$quantile))
  }
}
if (!is.null(opt$prs_extremes) && opt$prs_extremes != "") {
  cat("\n2. Extreme quantiles:\n")
  extreme_specs <- as.numeric(trimws(strsplit(opt$prs_extremes, ",")[[1]]))
  extreme_specs <- extreme_specs[!is.na(extreme_specs) & extreme_specs > 0 & extreme_specs < 50]
  for (pct in extreme_specs) {
    extreme_name <- paste0("extreme_", pct); col_name <- paste0("PRS_", extreme_name)
    cat("   - Top/Bottom ", pct, "% vs Middle\n", sep="")
    ee <- create_extreme_quantile_single(DT$PRS_std, pct)
    DT[[col_name]] <- ee$quantile; quantile_list[[extreme_name]] <- ee
    print(table(ee$quantile))
  }
}
# ------------------------------------------------------------------------------
# PER-ROLE QUANTILES (DISPLAY ONLY)
# ------------------------------------------------------------------------------
# The quantiles above are built from PRS_std — the PRIMARY score, which is the
# onset PRS as the driver supplies it. So every KM curve, extreme-group panel and
# absolute-risk stratum downstream described the ONSET score, not the outcome or
# progression score the manuscript is actually about.
#
# Build the same quantile schemes for all three roles and expose them as
# metadata$quantiles_by_role. These are DELIBERATELY NOT added to quantile_list:
# that object drives a coxph fit per entry (see the quantile-model loop below),
# each of which also lands in fitted_models, model_comparison and the C-index
# bootstrap. Tripling it would triple that cost and add meaningless rows.
quantiles_by_role <- list()
if (exists("prs_role_cols") && length(prs_role_cols) > 0) {
  role_q_specs <- quantile_specs
  cat("\n3. Per-role quantiles (display only; no extra models fitted):\n")
  for (role in names(prs_role_cols)) {
    scol <- prs_role_cols[[role]]
    if (is.null(scol) || !scol %in% names(DT)) next
    role_list <- list()
    for (spec in role_q_specs) {
      n_quant <- suppressWarnings(as.integer(spec))
      if (is.na(n_quant) || n_quant < 2) next
      qn  <- paste0("Q", n_quant)
      col <- paste0("PRS_", role, "_", qn)
      qq  <- create_standard_quantiles(DT[[scol]], n_quant, opt$prs_quantile_reference)
      DT[[col]] <- qq$quantile; DT[[paste0(col, "_numeric")]] <- qq$numeric
      role_list[[qn]] <- qq
    }
    if (exists("extreme_specs") && length(extreme_specs) > 0) {
      for (pct in extreme_specs) {
        en  <- paste0("extreme_", pct)
        col <- paste0("PRS_", role, "_", en)
        ee  <- create_extreme_quantile_single(DT[[scol]], pct)
        DT[[col]] <- ee$quantile; role_list[[en]] <- ee
      }
    }
    quantiles_by_role[[role]] <- role_list
    cat("   - ", role, " (", scol, "): ", paste(names(role_list), collapse=", "), "\n", sep="")
  }
}

if (length(quantile_specs) == 1) {
  DT$PRS_quantile <- DT[[paste0("PRS_Q", quantile_specs[1])]]
  DT$PRS_quantile_numeric <- DT[[paste0("PRS_Q", quantile_specs[1], "_numeric")]]
}
cat("✓ Quantile creation completed\n\n")

# ==============================================================================
# FIT MODELS (always start–stop Surv(START, STOP, PHENO))
# ==============================================================================
cat("==============================================================================\nFITTING MODELS\n==============================================================================\n\n")
surv_response <- "Surv(START, STOP, PHENO)"

# Final guarantee on the ACTUAL analysis sample. The float-safe exclusion applied
# earlier is not sufficient on its own: the core complete-case filter and the
# three-PRS alignment both remove rows AFTER it, and aeqSurv's second criterion
# divides by mean(abs(times)) — which shifts as rows leave. So re-ask aeqSurv here,
# where DT is finally fixed and every core model below fits on it.
{
  .n_before_aeq <- nrow(DT)
  DT <- .drop_degenerate_intervals(DT, "core M0-M7 sample")
  if (nrow(DT) < .n_before_aeq) {
    # n_events_final was computed during PRS alignment and is consumed downstream
    # (EPV report, core-sample QC row); recompute or those figures go stale.
    n_events_final <- sum(DT$PHENO == 1)
    record_stage(DT, "aeqsurv_resolvable")
  }
}

fitted_models <- list(); model_stats <- list()
has_clinical <- length(clinical_covs) > 0
# The standalone legacy Clinical/Full models (Model 3/4) are superseded by the
# systematic Clinical_M0-M7 ladder in the manuscript (multi-PRS) path, and they
# also lack strata_term. Fit them ONLY in legacy single-PRS mode.
fit_legacy_clinical <- has_clinical && !multi_prs_active

# Predictor set for every COXPH model formula = plain covariates + spline terms +
# (prospective only) the strata term. `mycovs` stays a plain-column vector (for
# cbind / VIF column selection); `extra_terms` (splines) is lm-safe; `strata_term`
# is coxph-only and must never reach the VIF lm(). model_covs is used only in
# coxph formula construction below, so folding strata in here is safe.
model_covs <- c(mycovs, extra_terms, strata_term)

cat("Model structure:\n")
cat("  Basic covariates:", ifelse(length(model_covs)>0, paste(model_covs, collapse=", "), "None"), "\n")
cat("  Clinical covariates:", ifelse(has_clinical, paste(clinical_covs, collapse=", "), "None"), "\n")
cat("  Clinical models:", ifelse(has_clinical, "YES (4 models)", "NO (2 models)"), "\n")
if (!is.null(opt$strata)) cat("  Stratification:", opt$strata, "\n")
cat("  Start–stop model: TRUE\n\n")

# Base
cat("MODEL 1: Base (basic covariates)\n---------------------------------\n")
base_formula <- if (length(model_covs) > 0) safe_formula(surv_response, model_covs, opt$strata) else as.formula(sprintf("%s ~ 1", surv_response))
cat("Formula:", deparse(base_formula), "\n")
fit_base <- coxph(base_formula, data=DT, x=TRUE, y=TRUE, model=TRUE)
cat("✓ Fitted | C-index:", round(summary(fit_base)$concordance[1], 4), "\n\n")
fitted_models$base <- fit_base; model_stats$base <- extract_comprehensive_stats(fit_base, "Base")

# Minimal
cat("MODEL 2: Minimal (PRS + basic covariates)\n-----------------------------------------\n")
minimal_formula <- if (length(model_covs) > 0) safe_formula(surv_response, c("PRS_std", model_covs), opt$strata) else safe_formula(surv_response, "PRS_std", opt$strata)
cat("Formula:", deparse(minimal_formula), "\n")
fit_minimal <- coxph(minimal_formula, data=DT, x=TRUE, y=TRUE, model=TRUE)
cat("✓ Fitted | C-index:", round(summary(fit_minimal)$concordance[1], 4), "\n\n")
fitted_models$minimal <- fit_minimal; model_stats$minimal <- extract_comprehensive_stats(fit_minimal, "Minimal")
prs_coef_minimal <- model_stats$minimal$coefficients[variable == "PRS_std"]

# Clinical & Full (legacy single-PRS only; superseded by the Clinical_M0-M7 ladder)
if (fit_legacy_clinical) {
  cat("MODEL 3: Clinical (basic + clinical)\n-------------------------------------\n")
  all_covs <- c(model_covs, clinical_covs)
  clinical_formula <- safe_formula(surv_response, all_covs, opt$strata)
  cat("Formula:", deparse(clinical_formula), "\n")
  fit_clinical <- coxph(clinical_formula, data=DT, x=TRUE, y=TRUE, model=TRUE)
  cat("✓ Fitted | C-index:", round(summary(fit_clinical)$concordance[1], 4), "\n\n")
  fitted_models$clinical <- fit_clinical; model_stats$clinical <- extract_comprehensive_stats(fit_clinical, "Clinical")
  
  cat("MODEL 4: Full (PRS + basic + clinical)\n---------------------------------------\n")
  full_formula <- safe_formula(surv_response, c("PRS_std", all_covs), opt$strata)
  cat("Formula:", deparse(full_formula), "\n")
  fit_full <- coxph(full_formula, data=DT, x=TRUE, y=TRUE, model=TRUE)
  cat("✓ Fitted | C-index:", round(summary(fit_full)$concordance[1], 4), "\n\n")
  fitted_models$full <- fit_full; model_stats$full <- extract_comprehensive_stats(fit_full, "Full")
  prs_coef_full <- model_stats$full$coefficients[variable == "PRS_std"]
} else {
  if (has_clinical && multi_prs_active)
    cat("MODEL 3 & 4 (legacy Clinical/Full): skipped — superseded by the Clinical_M0-M7 ladder\n\n")
  else
    cat("MODEL 3 & 4: Skipped (no clinical covariates specified)\n\n")
}

# Quantile models
cat("QUANTILE MODELS:\n-----------------\n")
quantile_models <- list(); quantile_covs <- model_covs
for (q_name in names(quantile_list)) {
  quant_var <- paste0("PRS_", q_name)
  if (!quant_var %in% names(DT)) next
  cat("Fitting quantile model:", q_name, "\n")
  quant_formula <- if (length(quantile_covs) > 0) safe_formula(surv_response, c(quant_var, quantile_covs), opt$strata)
  else safe_formula(surv_response, quant_var, opt$strata)
  fit_quant <- coxph(quant_formula, data=DT, x=TRUE, y=TRUE, model=TRUE)
  cat("  ✓ C-index:", round(summary(fit_quant)$concordance[1], 4), "\n")
  fitted_models[[paste0("quantile_", q_name)]] <- fit_quant
  model_stats[[paste0("quantile_", q_name)]] <- extract_comprehensive_stats(fit_quant, paste0("Quantile_", q_name))
  quantile_models[[q_name]] <- fit_quant
}
cat("\n")

# ==============================================================================
# MULTI-PRS MODEL COMBINATIONS
# ==============================================================================
if (multi_prs_active && !is.null(model_specs) && length(model_specs) > 0) {
  cat("==============================================================================\nMULTI-PRS MODEL COMBINATIONS\n==============================================================================\n\n")
  for (mn in names(model_specs)) {
    prs_cols_mn <- unname(Filter(Negate(is.null), lapply(model_specs[[mn]], function(lbl) {
      col <- multi_prs_std_cols[[lbl]]
      if (!is.null(col) && col %in% names(DT)) col else NULL
    })))
    if (!length(prs_cols_mn)) { cat("  \u26a0 Skipping", mn, ": no valid PRS columns\n"); next }
    all_preds  <- if (length(model_covs) > 0) c(prs_cols_mn, model_covs) else prs_cols_mn
    formula_mn <- safe_formula(surv_response, all_preds, opt$strata)
    cat("Fitting:", mn, "\n  Formula:", deparse(formula_mn), "\n")
    fit_mn <- tryCatch(
      coxph(formula_mn, data=DT, x=TRUE, y=TRUE, model=TRUE),
      error=function(e) { cat("  Error:", e$message, "\n"); NULL }
    )
    if (!is.null(fit_mn)) {
      cat("  \u2713 C-index:", round(summary(fit_mn)$concordance[1], 4), "\n\n")
      fitted_models[[mn]] <- fit_mn
      model_stats[[mn]]   <- extract_comprehensive_stats(fit_mn, mn)
    }
  }

  # ENFORCE the mandatory manuscript models. The loop tolerates a single failure
  # (tryCatch -> NULL) for robustness, but a distributed production run must not
  # return silently with M5 or M7 absent. Hard-stop listing any missing core model.
  required_core_models <- c("Onset","Outcome","Progression","Onset_Progression",
                            "Outcome_Progression","Onset_Outcome","Onset_Progression_Outcome")
  missing_core <- setdiff(required_core_models, names(fitted_models))
  if (length(missing_core))
    stop("Required M0-M7 core models failed to fit: ", paste(missing_core, collapse=", "),
         ". The run is incomplete; do not return partial results.")

  # ---------------------------------------------------------------------------
  # PRS ORIENTATION DETECTION (advisory; never auto-flips). A correctly aligned PRS
  # has a POSITIVE adjusted association with its outcome. Key the check off each
  # single-PRS model's fitted coefficient (the published, covariate/age-adjusted
  # estimand — the marginal correlation can disagree with it under confounding), and
  # warn only when it is SIGNIFICANTLY negative, so a truly-null score (whose sign is
  # noise) stays quiet. The fix is a PRESPECIFIED --prs_signs re-run, never an
  # outcome-driven flip. Signs are also visible centrally in all_coefficients.csv.
  for (role in names(prs_role_cols)) {
    std_col <- prs_role_cols[[role]]; label <- sub("_std$", "", std_col)
    ms <- model_stats[[label]]
    if (is.null(ms) || is.null(ms$coefficients)) next
    crow <- ms$coefficients[variable == std_col]
    if (!nrow(crow)) next
    b <- crow$coef[1]; pv <- crow$p[1]; hr <- crow$HR[1]
    if (is.finite(b) && b < 0 && is.finite(pv) && pv < 0.05) {
      cat(sprintf("  ⚠ WARNING: the %s PRS has a significantly NEGATIVE adjusted association with the outcome (HR/SD = %.3f, p = %.2g) after any --prs_signs.\n", role, hr, pv))
      cat("      A correctly aligned PRS is POSITIVELY associated with its outcome. If this reflects a known\n")
      cat("      reversed allele/sign convention for this biobank (confirm against the discovery GWAS / other\n")
      cat("      biobanks, NOT this sample's outcome), re-run with the prespecified sign, e.g. --prs_signs=-1,...\n")
      cat("      Do NOT flip on this in-sample association alone — that is target-outcome leakage.\n")
    }
  }

  # ---------------------------------------------------------------------------
  # CLINICAL LADDER (gated): parallel Clinical_M0..M7 = core M-model + clinical
  # risk factors. Dormant unless --clinical_covariates supplied. Base is included
  # as Clinical_Base so the clinical ladder has its own M0.
  # ---------------------------------------------------------------------------
  if (has_clinical) {
    cat("CLINICAL MULTI-PRS LADDER (core + clinical risk factors)\n-------------------------------------------------------\n")
    clin_covs_all <- c(model_covs, clinical_covs)  # model_covs already carries strata/splines
    # The clinical ladder runs on its OWN complete-case sample. Clinical covariates
    # are deliberately excluded from the core filter (see the complete-case block
    # above) so optional clinical missingness cannot shrink Core M0-M7; that means
    # the ladder must apply the restriction here instead.
    clin_present <- intersect(clinical_covs, names(DT))
    DT_clinical <- if (length(clin_present)) DT[complete.cases(DT[, ..clin_present])] else DT
    # DT_clinical is a strict subset of DT, so it has its own time set and its own
    # mean(abs(times)) — the guarantee applied to DT does not carry over. Without
    # this the clinical fits below, which are tryCatch-wrapped, would turn an
    # aeqSurv collapse into a silently skipped model rather than a loud failure.
    DT_clinical <- .drop_degenerate_intervals(DT_clinical, "clinical ladder sample")
    n_clin_dropped <- nrow(DT) - nrow(DT_clinical)
    cat("  Clinical sample: N=", nrow(DT_clinical), " events=", sum(DT_clinical$PHENO),
        " (", n_clin_dropped, " of ", nrow(DT), " core rows dropped for missing clinical covariates)\n", sep="")
    if (n_clin_dropped > 0)
      cat("  ⚠ Clinical_* rows are NOT directly comparable with core M0-M7 rows: different N.\n")
    if (nrow(DT_clinical) < 2 || sum(DT_clinical$PHENO) < 1) {
      cat("  ⚠ Clinical sample too small after complete-case; skipping the clinical ladder.\n\n")
      has_clinical <- FALSE
    }
  }
  if (has_clinical) {
    # Clinical Base (M0 + clinical)
    fitc_base <- tryCatch(coxph(safe_formula(surv_response, clin_covs_all, opt$strata), data=DT_clinical, x=TRUE, y=TRUE, model=TRUE),
                          error=function(e){cat("  Clinical_Base error:", e$message, "\n"); NULL})
    if (!is.null(fitc_base)) { fitted_models[["Clinical_Base"]] <- fitc_base; model_stats[["Clinical_Base"]] <- extract_comprehensive_stats(fitc_base, "Clinical_Base") }
    for (mn in names(model_specs)) {
      prs_cols_mn <- unname(Filter(Negate(is.null), lapply(model_specs[[mn]], function(lbl) {
        col <- multi_prs_std_cols[[lbl]]; if (!is.null(col) && col %in% names(DT_clinical)) col else NULL })))
      if (!length(prs_cols_mn)) next
      cmn <- paste0("Clinical_", mn)
      fit_cmn <- tryCatch(coxph(safe_formula(surv_response, c(prs_cols_mn, clin_covs_all), opt$strata), data=DT_clinical, x=TRUE, y=TRUE, model=TRUE),
                          error=function(e){cat("  ", cmn, "error:", e$message, "\n"); NULL})
      if (!is.null(fit_cmn)) {
        cat("Fitting:", cmn, "| C-index:", round(summary(fit_cmn)$concordance[1], 4), "\n")
        fitted_models[[cmn]] <- fit_cmn
        model_stats[[cmn]]   <- extract_comprehensive_stats(fit_cmn, cmn)
      }
    }
    cat("\n")
  }
}

# ==============================================================================
# MODEL COMPARISONS
# ==============================================================================
cat("MODEL COMPARISONS (vs Base)\n---------------------------\n")
comparisons_list <- list()
lrt_minimal_vs_base <- anova(fit_base, fit_minimal, test="Chisq")
comparisons_list$minimal_vs_base <- lrt_minimal_vs_base
cat("Minimal vs Base: χ²=", round(lrt_minimal_vs_base[2,"Chisq"],3),
    ", p=", format.pval(lrt_minimal_vs_base[2,"Pr(>|Chi|)"], digits=3), "\n", sep="")
if (fit_legacy_clinical) {
  lrt_clinical_vs_base <- anova(fit_base, fit_clinical, test="Chisq"); comparisons_list$clinical_vs_base <- lrt_clinical_vs_base
  cat("Clinical vs Base: χ²=", round(lrt_clinical_vs_base[2,"Chisq"],3),
      ", p=", format.pval(lrt_clinical_vs_base[2,"Pr(>|Chi|)"], digits=3), "\n", sep="")
  lrt_full_vs_base <- anova(fit_base, fit_full, test="Chisq"); comparisons_list$full_vs_base <- lrt_full_vs_base
  cat("Full vs Base: χ²=", round(lrt_full_vs_base[2,"Chisq"],3),
      ", p=", format.pval(lrt_full_vs_base[2,"Pr(>|Chi|)"], digits=3), "\n", sep="")
  lrt_full_vs_clinical <- anova(fit_clinical, fit_full, test="Chisq"); comparisons_list$full_vs_clinical <- lrt_full_vs_clinical
  cat("Full vs Clinical: χ²=", round(lrt_full_vs_clinical[2,"Chisq"],3),
      ", p=", format.pval(lrt_full_vs_clinical[2,"Pr(>|Chi|)"], digits=3), " (PRS incremental value)\n", sep="")
}

# ==============================================================================
# NESTED-MODEL LIKELIHOOD-RATIO LADDER  (M0-M7)
# ==============================================================================
# Emits every nested pair among the manuscript models, flagging the five key
# incremental tests. All models are fitted on an identical row set (enforced
# above), which anova.coxph requires.
lrt_ladder <- data.table()
if (multi_prs_active && exists("model_map") && nrow(model_map) > 1) {
  cat("==============================================================================\n")
  cat("NESTED-MODEL LIKELIHOOD-RATIO LADDER\n")
  cat("==============================================================================\n\n")

  get_fit <- function(nm) if (identical(nm, "Base")) fitted_models$base else fitted_models[[nm]]
  terms_set <- function(s) if (is.na(s) || !nzchar(s)) character(0) else strsplit(s, "+", fixed=TRUE)[[1]]

  # (reduced_id, full_id, label, is_primary). M7 is the full (three-PRS) model;
  # the leave-one-out triad drops each score in turn from M7.
  # Two CO-PRIMARY tests: M4 vs M1 (progression beyond onset) and M5 vs M2
  # (progression beyond the overall-T2 score). M7 vs M6 is the fully-adjusted
  # confirmation (progression beyond both; key main-text, not formally primary).
  # Note M4 vs M1 is NOT adjusted for the outcome PRS, so its signal can still reflect
  # general T2 liability; M7 vs M6 is the stronger conditional test.
  key_tests <- list(
    list("M1","M4","M4 vs M1: progression beyond onset",                    TRUE ),
    list("M2","M5","M5 vs M2: progression beyond outcome",                  TRUE ),
    list("M6","M7","M7 vs M6: progression beyond onset + outcome",          FALSE),
    list("M5","M7","M7 vs M5: onset beyond outcome + progression",          FALSE),
    list("M4","M7","M7 vs M4: outcome beyond onset + progression",          FALSE)
  )
  key_lookup <- setNames(
    lapply(key_tests, function(k) list(label=k[[3]], primary=k[[4]])),
    vapply(key_tests, function(k) paste0(k[[2]], "_vs_", k[[1]]), character(1))
  )

  avail <- model_map[model == "Base" | model %in% names(fitted_models)]
  rows <- list()
  for (i in seq_len(nrow(avail))) {
    for (j in seq_len(nrow(avail))) {
      if (i == j) next
      a <- avail[i]; b <- avail[j]                     # a = reduced, b = full
      if (!identical(a$ladder, b$ladder)) next         # pair only WITHIN a ladder
      sa <- terms_set(a$prs_terms); sb <- terms_set(b$prs_terms)
      if (!(length(sa) < length(sb) && all(sa %in% sb))) next   # a must nest in b

      fa <- get_fit(a$model); fb <- get_fit(b$model)
      if (is.null(fa) || is.null(fb)) next

      lr <- tryCatch(anova(fa, fb, test="Chisq"), error=function(e) {
        cat("  ⚠ LRT failed for", b$model, "vs", a$model, ":", conditionMessage(e), "\n"); NULL
      })
      if (is.null(lr)) next

      is_clin <- identical(a$ladder, "clinical")
      tag <- paste0(if (is_clin) "Clinical_" else "", b$model_id, "_vs_", a$model_id)
      kt  <- if (is_clin) NULL else key_lookup[[tag]]   # key tests defined on the core ladder

      rows[[length(rows)+1L]] <- data.table(
        test              = tag,
        ladder            = a$ladder,
        model_full        = b$model,
        model_reduced     = a$model,
        model_id_full     = b$model_id,
        model_id_reduced  = a$model_id,
        prs_added         = paste(setdiff(sb, sa), collapse="+"),
        df                = lr[2, "Df"],
        chisq             = lr[2, "Chisq"],
        p                 = lr[2, "Pr(>|Chi|)"],
        cindex_full       = summary(fb)$concordance[1],
        cindex_reduced    = summary(fa)$concordance[1],
        delta_cindex      = summary(fb)$concordance[1] - summary(fa)$concordance[1],
        rsq_RS_full       = summary(fb)$rsq[1],
        rsq_RS_reduced    = summary(fa)$rsq[1],
        delta_rsq_RS      = summary(fb)$rsq[1] - summary(fa)$rsq[1],
        n                 = fb$n,
        n_events          = fb$nevent,
        key_test          = if (!is.null(kt)) kt$label else NA_character_,
        is_primary        = if (!is.null(kt)) isTRUE(kt$primary) else FALSE
      )
    }
  }

  if (length(rows)) {
    lrt_ladder <- rbindlist(rows, fill=TRUE)
    setorder(lrt_ladder, -is_primary, key_test, model_id_full, model_id_reduced, na.last=TRUE)

    cat("Key incremental tests:\n")
    kk <- lrt_ladder[!is.na(key_test)]
    if (nrow(kk)) {
      for (i in seq_len(nrow(kk))) {
        cat(sprintf("  %-52s chisq=%8.3f  df=%d  p=%-10s %s\n",
                    kk$key_test[i], kk$chisq[i], as.integer(kk$df[i]),
                    format.pval(kk$p[i], digits=3),
                    ifelse(kk$is_primary[i], "[PRIMARY]", "")))
      }
    } else cat("  (none available with the models fitted)\n")

    lrt_file <- file.path(opt$outdir, paste0(prefix, "model_lrt_ladder.csv"))
    fwrite(lrt_ladder, lrt_file)
    cat("\n✓ Saved nested-model LRT ladder:", basename(lrt_file),
        sprintf("(%d nested pairs)\n\n", nrow(lrt_ladder)))
  } else {
    cat("No nested model pairs available.\n\n")
  }

  # run_status.csv — machine-readable completeness manifest for central QC. Reached
  # only after the mandatory-model stop() passed, so core_M0_M7_complete is TRUE.
  .avail <- function(pat) if (nrow(lrt_ladder) &&
                             "key_test" %in% names(lrt_ladder) && "p" %in% names(lrt_ladder))
                            nrow(lrt_ladder[grepl(pat, key_test) & is.finite(p)]) > 0 else FALSE
  m4vm1 <- isTRUE(.avail("M4 vs M1")); m5vm2 <- isTRUE(.avail("M5 vs M2")); m7vm6 <- isTRUE(.avail("M7 vs M6"))
  run_status <- data.table(
    core_M0_M7_complete = TRUE,
    M4_vs_M1_available  = m4vm1,   # co-primary: progression beyond onset
    M5_vs_M2_available  = m5vm2,   # co-primary: progression beyond outcome
    M7_vs_M6_available  = m7vm6,   # fully-adjusted confirmation (beyond both)
    n_core = nrow(DT), events_core = if (exists("n_events_final")) n_events_final else sum(DT$PHENO==1),
    analysis_mode = opt$analysis_mode,
    # run_complete reflects the two CO-PRIMARY tests (M4-M1, M5-M2) plus the
    # fully-adjusted confirmation (M7-M6): the mandatory-model stop() already guards
    # fitting, so this flags the case where a key LRT is degenerate/unavailable.
    run_complete = m4vm1 && m5vm2 && m7vm6
  )
  fwrite(run_status, file.path(opt$outdir, paste0(prefix, "run_status.csv")))
  cat("✓ Saved run_status.csv (M4vsM1:", m4vm1, "| M5vsM2:", m5vm2, "| M7vsM6:", m7vm6, ")\n")
  if (!run_status$run_complete)
    cat("  ⚠ run_complete=FALSE: a required key LRT (M4 vs M1 / M5 vs M2 / M7 vs M6) is unavailable.\n")
  cat("\n")
}

# ==============================================================================
# PRS CORRELATION, COLLINEARITY (VIF)
# ==============================================================================
if (multi_prs_active && length(multi_prs_std_cols) > 1) {
  cat("==============================================================================\n")
  cat("PRS CORRELATION & COLLINEARITY\n")
  cat("==============================================================================\n\n")

  lbls <- names(multi_prs_std_cols)
  cols <- unlist(multi_prs_std_cols, use.names=FALSE)
  ok   <- cols %in% names(DT)
  lbls <- lbls[ok]; cols <- cols[ok]

  prs_id_of_label <- function(l) {
    r <- role_of_label[[l]]
    if (is.null(r) || is.na(r)) return(l)
    id <- prs_id_of_role[[r]]
    if (is.null(id) || is.na(id)) l else id
  }

  # ---- pairwise correlations -------------------------------------------------
  cor_rows <- list()
  if (length(cols) > 1) {
    for (i in 1:(length(cols)-1)) for (j in (i+1):length(cols)) {
      x <- DT[[cols[i]]]; y <- DT[[cols[j]]]
      cor_rows[[length(cor_rows)+1L]] <- data.table(
        prs_a      = lbls[i],            prs_b      = lbls[j],
        prs_id_a   = prs_id_of_label(lbls[i]), prs_id_b = prs_id_of_label(lbls[j]),
        pearson    = suppressWarnings(cor(x, y, use="complete.obs", method="pearson")),
        spearman   = suppressWarnings(cor(x, y, use="complete.obs", method="spearman")),
        n          = sum(is.finite(x) & is.finite(y))
      )
    }
  }
  if (length(cor_rows)) {
    prs_cor <- rbindlist(cor_rows)
    for (i in seq_len(nrow(prs_cor)))
      cat(sprintf("  %-12s vs %-12s  Pearson %+.4f   Spearman %+.4f\n",
                  prs_cor$prs_id_a[i], prs_cor$prs_id_b[i],
                  prs_cor$pearson[i], prs_cor$spearman[i]))
    f <- file.path(opt$outdir, paste0(prefix, "prs_correlation.csv"))
    fwrite(prs_cor, f); cat("\n✓ Saved:", basename(f), "\n\n")
  }

  # ---- variance inflation factors -------------------------------------------
  # VIF_j = 1/(1 - R2_j) from regressing PRS_j on the other PRS + covariates.
  # Computed directly with lm() to avoid adding a dependency on `car`.
  if (length(cols) > 1) {
    vif_rows <- list()
    for (i in seq_along(cols)) {
      others <- setdiff(cols, cols[i])
      # CORE VIF: adjust for the other PRSs and CORE covariates only. Including
      # optional clinical_covs would make lm() complete-case on them, shrinking the
      # VIF sample and changing the core PRS collinearity estimate. (A separate
      # clinical VIF on DT_clinical is emitted below when the clinical ladder runs.)
      rhs <- c(others, mycovs)
      rhs <- intersect(rhs, names(DT))            # plain columns only
      rhs <- c(rhs, extra_terms)                  # + age spline expression (VIF adjusts for it too)
      if (!length(rhs)) next
      fml <- as.formula(paste(cols[i], "~", paste(rhs, collapse=" + ")))
      r2  <- tryCatch(summary(lm(fml, data=DT))$r.squared, error=function(e) NA_real_)
      vif_rows[[length(vif_rows)+1L]] <- data.table(
        prs = lbls[i], prs_id = prs_id_of_label(lbls[i]),
        r_squared_on_others = r2,
        VIF = if (is.finite(r2) && r2 < 1) 1/(1-r2) else NA_real_
      )
    }
    if (length(vif_rows)) {
      prs_vif <- rbindlist(vif_rows)
      for (i in seq_len(nrow(prs_vif)))
        cat(sprintf("  VIF %-12s %.3f%s\n", prs_vif$prs_id[i], prs_vif$VIF[i],
                    ifelse(is.finite(prs_vif$VIF[i]) && prs_vif$VIF[i] > 5, "   ⚠ >5", "")))
      f <- file.path(opt$outdir, paste0(prefix, "prs_vif.csv"))
      fwrite(prs_vif, f); cat("\n✓ Saved:", basename(f), "\n\n")
    }
    # Optional CLINICAL VIF: the three PRSs adjusted for core + clinical covariates,
    # on the clinical complete-case sample. Separate file so the core VIF stays on
    # the full core sample.
    if (exists("has_clinical") && isTRUE(has_clinical) && exists("DT_clinical") &&
        length(clinical_covs) && nrow(DT_clinical) > length(cols) + 2) {
      cvif_rows <- list()
      for (i in seq_along(cols)) {
        rhs_c <- c(setdiff(cols, cols[i]), mycovs, clinical_covs)
        rhs_c <- intersect(rhs_c, names(DT_clinical)); rhs_c <- c(rhs_c, extra_terms)
        if (!length(rhs_c)) next
        fml_c <- as.formula(paste(cols[i], "~", paste(rhs_c, collapse=" + ")))
        r2c <- tryCatch(summary(lm(fml_c, data=DT_clinical))$r.squared, error=function(e) NA_real_)
        cvif_rows[[length(cvif_rows)+1L]] <- data.table(
          prs = lbls[i], prs_id = prs_id_of_label(lbls[i]),
          r_squared_on_others = r2c,
          VIF = if (is.finite(r2c) && r2c < 1) 1/(1-r2c) else NA_real_,
          n = nrow(DT_clinical))
      }
      if (length(cvif_rows)) {
        fc <- file.path(opt$outdir, paste0(prefix, "prs_vif_clinical.csv"))
        fwrite(rbindlist(cvif_rows), fc); cat("✓ Saved:", basename(fc), "(clinical sample VIF)\n\n")
      }
    }
  }
}

# ==============================================================================
# HARRELL'S CIs & (optional) UNO's C
# ==============================================================================
cat("==============================================================================\nCALCULATING CONFIDENCE INTERVALS FOR HARRELL'S C-INDEX\n==============================================================================\n\n")
comparison_table <- rbindlist(lapply(model_stats, function(x) x$fit_stats))
# Attach manuscript labels (M0-M7) from the single source of truth built above.
if (exists("model_map"))
  comparison_table <- merge(comparison_table, model_map[, .(model, model_id, prs_terms, ladder)],
                            by="model", all.x=TRUE)
# Deltas are taken against the ladder-specific M0: core rows vs Base, clinical
# rows vs Clinical_Base (so Clinical_M* deltas reflect PRS increments over the
# clinical baseline, not over the core Base).
core_c0 <- comparison_table[model=="Base", concordance]
core_a0 <- comparison_table[model=="Base", AIC]
clin_c0 <- if ("Clinical_Base" %in% comparison_table$model) comparison_table[model=="Clinical_Base", concordance] else NA_real_
clin_a0 <- if ("Clinical_Base" %in% comparison_table$model) comparison_table[model=="Clinical_Base", AIC] else NA_real_
is_clin_row <- if ("ladder" %in% names(comparison_table)) comparison_table$ladder == "clinical" else rep(FALSE, nrow(comparison_table))
is_clin_row[is.na(is_clin_row)] <- FALSE
comparison_table[, `:=`(
  delta_concordance = concordance - fifelse(is_clin_row, clin_c0, core_c0),
  delta_AIC         = AIC         - fifelse(is_clin_row, clin_a0, core_a0)
)]

ci_results <- list()
main_models <- list(Base=fit_base, Minimal=fit_minimal)
if (fit_legacy_clinical) { main_models$Clinical <- fit_clinical; main_models$Full <- fit_full }
# Include the manuscript M0-M7 models (and, when present, the Clinical_M* ladder),
# so M2/M5/M6/M7 — the central comparisons — get C-index CIs and Uno's C, not just
# raw concordance. calculate_harrells_ci uses the fast SE path for large samples
# (n_case_control_pairs > 20000) and only bootstraps small ones, so this stays
# cheap on well-powered trajectories.
if (exists("model_map")) {
  for (mn in model_map$model) {
    fit_mn <- fitted_models[[mn]]
    if (!is.null(fit_mn) && inherits(fit_mn, "coxph") && is.null(main_models[[mn]]))
      main_models[[mn]] <- fit_mn
  }
}
for (model_name in names(main_models)) {
  cat(sprintf("\n%s model:\n", model_name))
  # Each model must be evaluated on the data it was FITTED on. Clinical_* models
  # use DT_clinical (a subset of DT), so scoring them against DT would mix samples.
  eval_dt <- if (grepl("^Clinical", model_name) && exists("DT_clinical")) DT_clinical else DT
  ci_results[[model_name]] <- calculate_harrells_ci(main_models[[model_name]], eval_dt, n_boot=1000, seed=42, max_pairs_for_bootstrap=20000)
}
if (length(quantile_models) > 0) {
  cat("\nQuantile models:\n")
  for (q_name in names(quantile_models)) {
    cat(sprintf("\n%s:\n", q_name))
    ci_results[[paste0("Quantile_", q_name)]] <- calculate_harrells_ci(quantile_models[[q_name]], DT, n_boot=1000, seed=42, max_pairs_for_bootstrap=20000)
  }
}
ci_table <- rbindlist(lapply(names(ci_results), function(m) {
  data.table(
    model = m,
    harrells_c = ci_results[[m]]$c_index,
    harrells_CI_lower = ci_results[[m]]$ci_lower,
    harrells_CI_upper = ci_results[[m]]$ci_upper,
    harrells_se = ci_results[[m]]$se,
    n_case_control_pairs = ci_results[[m]]$n_case_control_pairs,
    harrells_CI_method = ci_results[[m]]$method
  )
}))
comparison_table <- merge(comparison_table, ci_table, by="model", all.x=TRUE)

# Pooled Uno's C is invalid for stratified prospective models: predict(type="lp")
# ignores the per-stratum baselines, so a single pooled UnoC misrepresents the
# stratified risk. Skip it (LRTs/coefficients stay valid); use downstream
# stratum-specific metrics instead.
uno_skipped_stratified <- opt$calculate_uno && exists("is_prospective") && isTRUE(is_prospective) &&
                          exists("strata_term") && length(strata_term) > 0
if (uno_skipped_stratified) {
  cat("\n⚠ Skipping pooled Uno's C: prospective model is stratified by T1_STATUS\n",
      "  (pooled Uno ignores per-stratum baselines). Use downstream stratum-specific metrics.\n\n", sep="")
}
if (opt$calculate_uno && !uno_skipped_stratified) {
  cat("\n==============================================================================\nCALCULATING UNO'S C-INDEX FOR ALL MODELS\n==============================================================================\n\n")
  if (!requireNamespace("survAUC", quietly = TRUE)) stop("Package 'survAUC' is required for Uno's C.")
  uno_results <- list()
  for (model_name in names(main_models)) {
    cat(sprintf("\n  %s model:\n", model_name))
    eval_dt <- if (grepl("^Clinical", model_name) && exists("DT_clinical")) DT_clinical else DT
    uno_results[[model_name]] <- calculate_uno_cindex(
      model=main_models[[model_name]], data=eval_dt, time_col="TIME", status_col="PHENO",
      min_events=min_events_val, min_at_risk=min_at_risk_val, ghat_threshold=ghat_threshold_val,
      skip_ci=!opt$uno_bootstrap, n_boot=opt$uno_bootstrap_n, seed=100+match(model_name, names(main_models)),
      tau=opt$tau
    )
  }
  if (length(quantile_models) > 0) {
    for (i in seq_along(quantile_models)) {
      q_name <- names(quantile_models)[i]; model_obj <- quantile_models[[q_name]]
      uno_results[[paste0("Quantile_", q_name)]] <- calculate_uno_cindex(
        model=model_obj, data=DT, time_col="TIME", status_col="PHENO",
        min_events=min_events_val, min_at_risk=min_at_risk_val, ghat_threshold=ghat_threshold_val,
        skip_ci=!opt$uno_bootstrap, n_boot=opt$uno_bootstrap_n, seed=142+i, tau=opt$tau
      )
    }
  }
  uno_table <- rbindlist(lapply(names(uno_results), function(m) {
    data.table(
      model=m, uno_cindex=uno_results[[m]]$cindex, uno_se=uno_results[[m]]$se,
      uno_CI_lower=uno_results[[m]]$ci_lower, uno_CI_upper=uno_results[[m]]$ci_upper,
      uno_method=uno_results[[m]]$method, uno_tau=uno_results[[m]]$tau
    )
  }))
  comparison_table <- merge(comparison_table, uno_table, by="model", all.x=TRUE)
  cat("\n✓ Uno's C-index calculated and merged\n")
  cat("  NOTE: Uno's C (and Harrell's C) here is APPARENT performance — models are\n",
      "  fit and evaluated in the same sample. Optimism correction / held-out-biobank\n",
      "  external validation is a downstream step, not done in this fitting script.\n\n", sep="")
}

# Column order
core_cols <- c("model","model_id","ladder","prs_terms","n","n_events","concordance","harrells_CI_lower","harrells_CI_upper","concordance_se","harrells_CI_method",
               "n_case_control_pairs","delta_concordance","AIC","delta_AIC","BIC","loglik_null","loglik_model",
               "rsquare","max_rsquare","wald_test","wald_p","logrank_test","logrank_p")
if (opt$calculate_uno) core_cols <- append(core_cols, c("uno_cindex","uno_CI_lower","uno_CI_upper","uno_se","uno_method","uno_tau"), after=7)
setcolorder(comparison_table, core_cols[core_cols %in% names(comparison_table)])

# ------------------------------------------------------------------------------
# Events-per-parameter, from the ACTUAL fitted full model (M7). This replaces the
# fixed "95 events = EPV 5" heuristic: the real parameter count includes spline df,
# factor levels, PCs, batch and (prospective) strata.
# ------------------------------------------------------------------------------
full_fit_for_epv <- fitted_models[["Onset_Progression_Outcome"]]
if (!is.null(full_fit_for_epv)) {
  n_params_full <- length(stats::coef(full_fit_for_epv))
  epv <- if (n_params_full > 0) n_events_final / n_params_full else NA_real_
  cat(sprintf("\nEvents-per-parameter (full model M7): %d events / %d params = %.1f\n",
              n_events_final, n_params_full, epv))
  if (is.finite(epv) && epv < 5)
    cat("  ⚠ EPV < 5: M7-based estimates are unstable; report M5 vs M2 with caution.\n")
  else if (is.finite(epv) && epv < 10)
    cat("  ⚠ EPV < 10: borderline power for the fully adjusted M7 tests.\n")
  metadata_epv <- list(n_events_final=n_events_final, n_params_full=n_params_full, events_per_parameter=epv)
} else metadata_epv <- list(n_events_final=n_events_final, n_params_full=NA_integer_, events_per_parameter=NA_real_)

# Clinical-ladder EPV, on the clinical complete-case sample (which may be smaller).
clin_full <- fitted_models[["Clinical_Onset_Progression_Outcome"]]
if (!is.null(clin_full) && exists("DT_clinical")) {
  ev_clin <- sum(DT_clinical$PHENO == 1)
  np_clin <- length(stats::coef(clin_full))
  epv_clin <- if (np_clin > 0) ev_clin / np_clin else NA_real_
  cat(sprintf("Events-per-parameter (Clinical M7): %d events / %d params = %.1f\n",
              ev_clin, np_clin, epv_clin))
  if (is.finite(epv_clin) && epv_clin < 5)
    cat("  ⚠ Clinical EPV < 5: the Clinical ladder is unstable after clinical complete-case restriction.\n")
  else if (is.finite(epv_clin) && epv_clin < 10)
    cat("  ⚠ Clinical EPV < 10: borderline power for the fully adjusted Clinical M7.\n")
  metadata_epv$clinical <- list(events=ev_clin, n_params=np_clin, events_per_parameter=epv_clin)
}

# ==============================================================================
# PROPORTIONAL HAZARDS TESTS
# ==============================================================================
cat("PROPORTIONAL HAZARDS TESTS\n---------------------------\n")
check_ph <- function(fit, name) {
  ph_test <- tryCatch(cox.zph(fit), error=function(e) NULL)
  if (is.null(ph_test)) return(NULL)
  ph_dt <- data.table(
    model=name, variable=rownames(ph_test$table),
    chisq=ph_test$table[, "chisq"], df=ph_test$table[, "df"], p=ph_test$table[, "p"]
  )
  ph_dt[, violates := p < 0.05]
  if (any(ph_dt$violates)) cat("⚠ ", name, ": PH violated for ", paste(ph_dt[violates==TRUE, variable], collapse=", "), "\n", sep="")
  else cat("✓ ", name, ": PH assumption satisfied\n", sep="")
  list(test=ph_test, table=ph_dt)
}
ph_tests <- lapply(names(fitted_models), function(mn) check_ph(fitted_models[[mn]], mn))
names(ph_tests) <- names(fitted_models); cat("\n")

# ==============================================================================
# (OPTIONAL) PRS × STRATA INTERACTION
# ==============================================================================
if (opt$test_strata_interaction && !is.null(opt$strata)) {
  cat("PRS × STRATA INTERACTION TESTING\n---------------------------------\n")
  strata_var_for_interaction <- opt$strata
  minimal_no_strata_formula <- if (length(model_covs) > 0)
    safe_formula(surv_response, c("PRS_std", strata_var_for_interaction, model_covs), strata=NULL)
  else safe_formula(surv_response, c("PRS_std", strata_var_for_interaction), strata=NULL)
  fit_minimal_no_strata <- coxph(minimal_no_strata_formula, data=DT, x=TRUE, y=TRUE)
  interaction_formula <- if (length(model_covs) > 0)
    safe_formula(surv_response, c("PRS_std", strata_var_for_interaction, model_covs, paste0("PRS_std:", strata_var_for_interaction)), strata=NULL)
  else safe_formula(surv_response, c("PRS_std", strata_var_for_interaction, paste0("PRS_std:", strata_var_for_interaction)), strata=NULL)
  cat("   Formula:", deparse(interaction_formula), "\n")
  fit_minimal_interaction <- coxph(interaction_formula, data=DT, x=TRUE, y=TRUE)
  lrt_interaction <- anova(fit_minimal_no_strata, fit_minimal_interaction, test="Chisq")
  cat("   Interaction test: χ²=", round(lrt_interaction[2,"Chisq"],3),
      ", p=", format.pval(lrt_interaction[2,"Pr(>|Chi|)"], digits=3), "\n", sep="")
  if (has_clinical) {
    all_covs <- c(mycovs, clinical_covs, extra_terms, strata_term)
    full_no_strata_formula <- safe_formula(surv_response, c("PRS_std", strata_var_for_interaction, all_covs), strata=NULL)
    fit_full_no_strata <- coxph(full_no_strata_formula, data=DT, x=TRUE, y=TRUE)
    interaction_formula_full <- safe_formula(surv_response, c("PRS_std", strata_var_for_interaction, all_covs, paste0("PRS_std:", strata_var_for_interaction)), strata=NULL)
    cat("   Full formula:", deparse(interaction_formula_full), "\n")
    fit_full_interaction <- coxph(interaction_formula_full, data=DT, x=TRUE, y=TRUE)
    lrt_interaction_full <- anova(fit_full_no_strata, fit_full_interaction, test="Chisq")
    cat("   Full interaction test: χ²=", round(lrt_interaction_full[2,"Chisq"],3),
        ", p=", format.pval(lrt_interaction_full[2,"Pr(>|Chi|)"], digits=3), "\n", sep="")
  }
  cat("\n")
}

# ==============================================================================
# ASSEMBLE RESULTS
# ==============================================================================
cat("ASSEMBLING RESULTS\n------------------\n")
metadata <- list(
  timestamp = Sys.time(),
  prs_file = opt$prs_file,
  pheno_file = opt$pheno_file,
  n_obs = nrow(DT),
  n_events = sum(DT$PHENO),
  columns = list(
    start = "START",
    stop  = "STOP",
    time  = "TIME",  # analysis follow-up scale (STOP − START)
    status = "PHENO",
    prs_original = "PRS_original", prs_standardized = "PRS_std", prs_use = "PRS_std",
    basic_covariates = mycovs, clinical_covariates = clinical_covs,
    age_col = opt$age_col, sex_col = opt$sex_col, strata = opt$strata,
    time_exit_arg  = opt$time_exit,
    time_entry_arg = opt$time_entry
  ),
  # quantiles: legacy field, built from the PRIMARY score (kept for backward
  # compatibility). quantiles_by_role carries the same schemes for all three PRS
  # roles, which is what 02/06 should display.
  quantiles = quantile_list,
  quantiles_by_role = if (exists("quantiles_by_role") && length(quantiles_by_role)) quantiles_by_role else NULL,
  prs_role_columns  = if (exists("prs_role_cols") && length(prs_role_cols)) prs_role_cols else NULL,
  prs_info = list(
    correlation_with_pheno = prs_pheno_cor,
    correlation_with_pheno_raw = prs_pheno_cor_raw,
    # Direction provenance: as_supplied | prs_signs | auto(outcome).
    direction = if (exists("prs_signs_applied") && length(prs_signs_applied)) "prs_signs" else prs_direction,
    prs_signs = if (exists("prs_signs_applied")) prs_signs_applied else NULL,
    zscore_method = if (isTRUE(opt$zscore_by_controls)) "controls" else "whole_sample",
    restandardised_after_alignment = TRUE
  ),
  model_structure = list(has_clinical = has_clinical, n_models = length(fitted_models), start_stop = TRUE,
                         n_core = nrow(DT), events_core = sum(DT$PHENO),
                         n_clinical = if (exists("DT_clinical")) nrow(DT_clinical) else NA_integer_,
                         events_clinical = if (exists("DT_clinical")) sum(DT_clinical$PHENO) else NA_integer_),
  # Manuscript model labels M0-M7. Downstream scripts join this rather than
  # re-deriving the mapping, so it is defined in exactly one place.
  model_map = if (exists("model_map")) model_map else NULL,
  lrt_ladder = if (exists("lrt_ladder") && nrow(lrt_ladder)) lrt_ladder else NULL,
  n_dropped_incomplete_prs = if (exists("multi_prs_n_dropped")) multi_prs_n_dropped else 0L,
  analysis_mode = opt$analysis_mode,
  lag_days = opt$lag_days,
  n_events_total = if (exists("n_events_total")) n_events_total else NA_integer_,
  n_events_final = if (exists("n_events_final")) n_events_final else NA_integer_,
  events_per_parameter = if (exists("metadata_epv")) metadata_epv else NULL,
  dropped_covariates = if (exists("dropped_covs")) dropped_covs else character(0),
  strata_term = if (exists("strata_term") && length(strata_term)) strata_term else NA_character_,
  age_term = if (length(extra_terms)) paste(extra_terms, collapse=" + ") else NA_character_,
  formulas = lapply(fitted_models, formula),
  time_support = list(
    t_pass_auto = supported_times,
    auto_table = time_support_table,
    time_points_requested_raw = requested_time_points_raw,
    time_points_requested = requested_time_points,
    requested_audit_table = requested_ghat_table,
    requested_valid_time_points = requested_valid_times,
    eval_time_points = eval_time_points
  ),
  time_points_supported = supported_times,
  uno_supported_times   = supported_times,
  options = opt
)

results <- list(
  models = fitted_models,
  stats = model_stats,
  prs = list(
    coefficient_minimal = prs_coef_minimal,
    coefficient_full = if (exists("prs_coef_full")) prs_coef_full else NULL,
    quantile_coefficients = lapply(quantile_models, function(m) coef(summary(m)))
  ),
  comparisons = comparisons_list,
  ph_tests = ph_tests,
  metadata = metadata
)
cat("✓ Results assembled\n\n")

# ==============================================================================
# STUDY-DESIGN SUMMARIES (final analysis sample = the M0-M7 sample)
# ==============================================================================
record_stage(DT, "final_analysis")
summary_final <- summarise_cohort(DT, "final_analysis")
if (exists("summary_pre"))
  fwrite(summary_pre,   file.path(opt$outdir, paste0(prefix, "pre_alignment_summary.csv")))
fwrite(summary_final,   file.path(opt$outdir, paste0(prefix, "final_analysis_summary.csv")))
cat("✓ Saved final_analysis_summary.csv (authoritative Table 1 source)\n")

# One-row per-trajectory summary (identity + final-sample target-analysis metrics).
traj <- data.table(
  cohort = if (!is.null(opt$cohort)) opt$cohort else NA_character_,
  ancestry = if (!is.null(opt$ancestry)) opt$ancestry else NA_character_,
  trajectory = if (!is.null(opt$pheno_name)) opt$pheno_name else NA_character_,
  analysis_mode = opt$analysis_mode,
  lag_days = opt$lag_days
)
traj <- cbind(traj, summary_final[, setdiff(names(summary_final), "summary"), with=FALSE])
setnames(traj, "N", "N_final")
fwrite(traj, file.path(opt$outdir, paste0(prefix, "trajectory_summary.csv")))
cat("✓ Saved trajectory_summary.csv\n")

# Temporal-order QC among T2 events. Two files, on two samples:
#   source_phenotype_timing_qc.csv - the PRE-exclusion snapshot (keeps same-day
#     transitions, for phenotype QC);
#   final_analysis_timing_qc.csv   - the FINAL M0-M7 sample (no same-day, since
#     strict temporal order excludes them by design).
# Day-level bins are approximate when ages are integer-valued.
.timing_qc <- function(s, label, N_col) {
  ev <- s$PHENO == 1
  d <- s$TIME_SINCE_T1[ev] * 365.25   # interval in days, among events
  nE <- length(d)
  qc <- data.table(
    sample = label,
    trajectory = if (!is.null(opt$pheno_name)) opt$pheno_name else NA_character_,
    analysis_mode = opt$analysis_mode,
    age_precision = if (isTRUE(ages_are_integer)) "approximate (integer ages)" else "fractional"
  )
  qc[[N_col]] <- nE
  qc[, `:=`(
    N_same_day   = sum(d == 0, na.rm=TRUE),
    N_1_30d      = sum(d > 0    & d <= 30,  na.rm=TRUE),
    N_31_90d     = sum(d > 30   & d <= 90,  na.rm=TRUE),
    N_91_365d    = sum(d > 90   & d <= 365, na.rm=TRUE),
    N_over_365d  = sum(d > 365, na.rm=TRUE),
    pct_within_30d = if (nE) round(100*sum(d <= 30, na.rm=TRUE)/nE, 2) else NA_real_,
    pct_within_90d = if (nE) round(100*sum(d <= 90, na.rm=TRUE)/nE, 2) else NA_real_
  )]
  if (all(c("AGE_RECRUIT","AGE_T1","AGE_EXIT") %in% names(s))) {
    ar <- s$AGE_RECRUIT; a1 <- s$AGE_T1; a2 <- s$AGE_EXIT
    qc[, `:=`(
      N_T1_before_recruit = sum(a1 <= ar, na.rm=TRUE),
      N_T1_after_recruit  = sum(a1 >  ar, na.rm=TRUE),
      N_T2_before_recruit = sum(ev & a2 <= ar, na.rm=TRUE),
      N_historical        = sum(a1 <= ar & ev & a2 <= ar, na.rm=TRUE),
      N_prevalent_future  = sum(a1 <= ar & !(ev & a2 <= ar), na.rm=TRUE),
      N_incident          = sum(a1 >  ar, na.rm=TRUE)
    )]
    # Recruitment-age coverage of the ORIGINAL merge. Every count above is
    # conditional on having a recruitment age, so without this the table looks
    # complete even when the recruitment file covered a small minority. A percentage
    # is not a cell count, so it carries no disclosure risk.
    if (exists("recruit_coverage_pct"))
      qc[, `:=`(recruit_coverage_pct = round(recruit_coverage_pct, 2))]
  }
  # Disclosure: the interval bins are an EXHAUSTIVE partition that sums to the
  # retained event total, and the recruit bins likewise — so a single hidden bin is
  # reconstructable from the total and the other bins. Suppress ALL bins in a
  # partition when ANY bin in it is small (all-or-none), keeping the overall total.
  thr <- if (!is.null(opt$min_cell_count)) opt$min_cell_count else 10L
  interval_bins <- intersect(c("N_same_day","N_1_30d","N_31_90d","N_91_365d","N_over_365d"), names(qc))
  recruit_bins  <- intersect(c("N_T1_before_recruit","N_T1_after_recruit","N_T2_before_recruit",
                               "N_historical","N_prevalent_future","N_incident"), names(qc))
  any_small <- FALSE
  for (grp in list(interval_bins, recruit_bins)) {
    if (!length(grp)) next
    vals <- suppressWarnings(as.numeric(unlist(qc[1, ..grp])))
    if (any(is.finite(vals) & vals > 0 & vals < thr)) {
      for (cc in grp) qc[[cc]] <- NA_integer_
      any_small <- TRUE
    }
  }
  qc[, minimum_cell_count_pass := !any_small]
  qc
}
if (exists("timing_qc_snapshot") && nrow(timing_qc_snapshot) && "TIME_SINCE_T1" %in% names(timing_qc_snapshot)) {
  qc_src <- .timing_qc(timing_qc_snapshot, "source_phenotype", "N_events_snapshot")
  fwrite(qc_src, file.path(opt$outdir, paste0(prefix, "source_phenotype_timing_qc.csv")))
  cat("✓ Saved source_phenotype_timing_qc.csv (pre-exclusion; retains same-day)\n")

  # Final-sample QC: intersect the snapshot with the M0-M7 IIDs (DT), so the
  # intervals describe exactly the participants who entered the models. AGE_EXIT
  # etc. come from the snapshot; TIME_SINCE_T1 is carried in DT for prospective.
  if ("IID" %in% names(timing_qc_snapshot) && "IID" %in% names(DT)) {
    fin <- timing_qc_snapshot[IID %in% DT$IID]
    if (nrow(fin)) {
      qc_fin <- .timing_qc(fin, "final_analysis", "N_events_final")
      fwrite(qc_fin, file.path(opt$outdir, paste0(prefix, "final_analysis_timing_qc.csv")))
      cat("✓ Saved final_analysis_timing_qc.csv (M0-M7 sample; same-day excluded by design)\n")
    }
  }
}

# Attrition flow table. N_remaining/events_remaining are overall running totals
# (permitted); the per-filter removals can be small, so suppress an explicit small
# N_removed/events_removed (1..thr-1).
if (length(attrition_log)) {
  att <- rbindlist(attrition_log, fill=TRUE)
  att[, N_removed := c(NA_integer_, -diff(N_remaining))]
  att[, events_removed := c(NA_integer_, -diff(events_remaining))]
  thr_att <- if (!is.null(opt$min_cell_count)) opt$min_cell_count else 10L
  att[is.finite(N_removed) & N_removed > 0 & N_removed < thr_att, N_removed := NA_integer_]
  att[is.finite(events_removed) & events_removed > 0 & events_removed < thr_att, events_removed := NA_integer_]
  fwrite(att, file.path(opt$outdir, paste0(prefix, "trajectory_attrition.csv")))
  cat("✓ Saved trajectory_attrition.csv\n")
}
cat("\n")

# ==============================================================================
# SAVE OUTPUTS
# ==============================================================================
cat("SAVING OUTPUTS\n--------------\n")
data_file <- file.path(opt$outdir, paste0(prefix, "data_processed.rds")); saveRDS(DT, data_file); cat("✓", basename(data_file), "\n")
models_file <- file.path(opt$outdir, paste0(prefix, "fitted_models.rds")); saveRDS(results, models_file); cat("✓", basename(models_file), "\n")
metadata_file <- file.path(opt$outdir, paste0(prefix, "model_metadata.rds")); saveRDS(metadata, metadata_file); cat("✓", basename(metadata_file), "\n")

comparison_file <- file.path(opt$outdir, paste0(prefix, "model_comparison.csv"))
fwrite(comparison_table, comparison_file); cat("✓", basename(comparison_file), "\n")

all_coefs <- rbindlist(lapply(model_stats, function(x) x$coefficients))
if (exists("model_map"))
  all_coefs <- merge(all_coefs, model_map[, .(model, model_id, prs_terms, ladder)], by="model", all.x=TRUE)
fwrite(all_coefs, file.path(opt$outdir, paste0(prefix, "all_coefficients.csv"))); cat("✓", paste0(prefix, "all_coefficients.csv"), "\n")
ph_combined <- rbindlist(lapply(ph_tests, function(x) if(!is.null(x)) x$table else NULL))
if (nrow(ph_combined) > 0) { fwrite(ph_combined, file.path(opt$outdir, paste0(prefix, "ph_tests.csv"))); cat("✓", paste0(prefix, "ph_tests.csv"), "\n") }

ci_details <- data.table(
  model = names(ci_results),
  harrells_c = sapply(ci_results, function(x) x$c_index),
  harrells_ci_lower = sapply(ci_results, function(x) x$ci_lower),
  harrells_ci_upper = sapply(ci_results, function(x) x$ci_upper),
  harrells_se = sapply(ci_results, function(x) x$se),
  n_case_control_pairs = sapply(ci_results, function(x) x$n_case_control_pairs),
  harrells_method = sapply(ci_results, function(x) x$method)
)
# harrells C here is APPARENT (same-sample fit+eval); optimism/external validation downstream.
fwrite(ci_details, file.path(opt$outdir, paste0(prefix, "harrells_ci_details.csv"))); cat("✓", paste0(prefix, "harrells_ci_details.csv"), "\n")

if (opt$calculate_uno && exists("uno_results")) {
  uno_details <- data.table(
    model       = names(uno_results),
    uno_c       = sapply(uno_results, `[[`, "cindex"),
    uno_ci_lower= sapply(uno_results, `[[`, "ci_lower"),
    uno_ci_upper= sapply(uno_results, `[[`, "ci_upper"),
    uno_se      = sapply(uno_results, `[[`, "se"),
    uno_method  = sapply(uno_results, `[[`, "method"),
    uno_tau     = sapply(uno_results, `[[`, "tau")
  )
  fwrite(uno_details, file.path(opt$outdir, paste0(prefix, "uno_cindex_details.csv"))); cat("✓", paste0(prefix, "uno_cindex_details.csv"), "\n")
}

# Helper loader
helper_script <- sprintf('
# Load fitted models
load_fitted_models <- function(outdir="%s") { readRDS(file.path(outdir, "%sfitted_models.rds")) }
load_processed_data <- function(outdir="%s") { readRDS(file.path(outdir, "%sdata_processed.rds")) }
load_metadata <- function(outdir="%s") { readRDS(file.path(outdir, "%smodel_metadata.rds")) }
', opt$outdir, prefix, opt$outdir, prefix, opt$outdir, prefix)
helper_file <- file.path(opt$outdir, paste0(prefix, "load_models.R"))
writeLines(helper_script, helper_file); cat("✓", basename(helper_file), "\n\n")

# ==============================================================================
# SUMMARY
# ==============================================================================
cat("==============================================================================\nSUMMARY\n==============================================================================\n\n")
cat("✓ Analysis completed!\n\n")
cat("Models fitted:\n"); for (i in seq_along(fitted_models)) cat("  ", i, ". ", names(fitted_models)[i], "\n", sep=""); cat("\n")
if (!is.null(opt$strata)) cat("Stratification variable:", opt$strata, "\n\n")
cat("Start–stop used: TRUE (Surv(START, STOP, PHENO)); FOLLOWUP = STOP − START\n")
if (length(supported_times))        cat("Auto-supported time points (analysis scale): ", paste(round(supported_times,2), collapse=", "), "\n", sep="")
if (length(requested_time_points))  cat("Requested time points (analysis scale): ", paste(requested_time_points, collapse=", "), "\n", sep="")
if (length(requested_valid_times))  cat("Requested VALID time points: ", paste(round(requested_valid_times,2), collapse=", "), "\n", sep="")
if (length(eval_time_points))       cat("Eval time points (downstream): ", paste(round(eval_time_points,2), collapse=", "), "\n", sep="") else cat("Eval time points: NONE\n")
if (!is.null(opt$time_summary_cols) && nzchar(opt$time_summary_cols)) cat("Extra descriptive time columns: ", opt$time_summary_cols, "\n", sep="")
cat("\n")
cat("==============================================================================\nEnd time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("==============================================================================\n")
sink(type="message"); sink(type="output"); close(log_con)
cat("\nOutputs saved to:", opt$outdir, "\n")
cat("Log file:", log_file, "\n")
