# ==============================================================================
# prs_risk_utils.R - Shared risk-prediction helpers for the PRS evaluation pipeline
# ==============================================================================
# Sourced by 04_reclassification.R, 05_calibration.R and 06_cumulative_incidence.R.
#
# WHY THIS FILE EXISTS
#   Predicted absolute risk at a horizon was previously computed three different
#   ways, which could disagree on identical input:
#     04 : basehaz(fit)          + which.min(abs(t - .))   [nearest neighbour]
#     05 : survfit(fit)          + findInterval            [correct step lookup]
#     06 : survfit(fit, newdata) + which.min, indexing the wrong matrix axis
#
#   Two defects were present:
#
#   (1) TIME-SCALE MISMATCH. 01_fit_models.R fits start-stop models,
#       Surv(START, STOP, PHENO), so the baseline hazard is indexed on the
#       START/STOP scale (e.g. age). Evaluation horizons, however, are on the
#       follow-up scale TIME = STOP - START. Indexing the first with the second
#       silently returns the wrong risk whenever START != 0 - i.e. whenever
#       --time_entry is used for a landmark analysis.
#
#   (2) WRONG MATRIX AXIS. survfit(fit, newdata=)$surv is times x individuals.
#       Selecting S[, idx] with a *time* index returns one individual's whole
#       time-course instead of every individual's risk at t.
#
#   predict_risk_at_horizon() below fixes both. It returns risk conditional on
#   being event-free at entry:
#
#       r_i(h) = 1 - [ S0(START_i + h) / S0(START_i) ] ^ exp(lp_i)
#
#   When START = 0 this reduces exactly to 1 - S0(h)^exp(lp), so results from
#   runs without --time_entry are unchanged.
#
# Author: Ying Wang (yiwang@broadinstitute.org)
# ==============================================================================

suppressPackageStartupMessages(library(survival))
# Models fitted by 01 may contain ns(age) terms; splines must be attached wherever
# those models are predicted/re-evaluated (predict, basehaz, riskRegression::Score).
suppressPackageStartupMessages(library(splines))

# ------------------------------------------------------------------------------
# save_plot_cairo() — save a ggplot/base plot, preferring Cairo (nicer, needs
# X11/XQuartz) but falling back to base grDevices so headless machines still
# produce figures. Never aborts the caller. Shared by 04/05/06 (03 has its own
# copy since it does not source this helper).
# ------------------------------------------------------------------------------
save_plot_cairo <- function(plot_obj, filename, width, height, dpi = 300, formats = c("png")) {
  have_cairo <- requireNamespace("Cairo", quietly = TRUE)
  for (fmt in formats) {
    if (!fmt %in% c("png", "pdf")) { warning("Unsupported format: ", fmt); next }
    # Strip any existing extension, then append. Substituting only when the name
    # already ends in .png/.pdf meant an extensionless filename produced the same
    # path for every requested format (each overwriting the last).
    stem    <- sub("\\.(png|pdf)$", "", filename, ignore.case = TRUE)
    outfile <- paste0(stem, ".", fmt)
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

# ------------------------------------------------------------------------------
# Baseline survival S0(t) for a fitted coxph model, as step function(s).
# Non-stratified: returns list(strata = FALSE, base = list(time, surv)).
# Stratified:     returns list(strata = TRUE, var = <strata variable name>,
#                              bases = named list of per-stratum list(time, surv)),
#                 keyed by the basehaz strata labels (e.g. "T1_STATUS=Incident").
# Returns NULL if unavailable.
# ------------------------------------------------------------------------------
.baseline_survival <- function(fit) {
  bh <- tryCatch(survival::basehaz(fit, centered = TRUE), error = function(e) NULL)
  if (is.null(bh) || !nrow(bh)) return(NULL)
  if (!"strata" %in% names(bh)) {
    ord <- order(bh$time)
    return(list(strata = FALSE, base = list(time = bh$time[ord], surv = exp(-bh$hazard[ord]))))
  }
  # Stratified: one baseline per stratum label.
  labs <- as.character(bh$strata)
  bases <- lapply(split(bh, labs), function(d) {
    ord <- order(d$time); list(time = d$time[ord], surv = exp(-d$hazard[ord]))
  })
  svar <- tryCatch(survival::untangle.specials(stats::terms(fit), "strata")$vars,
                   error = function(e) NULL)
  list(strata = TRUE, var = svar, bases = bases)
}

# ------------------------------------------------------------------------------
# Step-function lookup: value of S0 at each requested time.
# Uses findInterval (last observed time <= t), never a nearest-neighbour match,
# so a horizon can never be served by a LATER event time.
# S0(t) = 1 for t before the first event time.
#
# extrapolate = FALSE (default): times beyond the LAST baseline event time return
# NA rather than the carried-forward final value. Carrying forward makes
# S0(exit)/S0(entry) collapse to 1 once both fall past the observed range, which
# reports an UNSUPPORTED horizon as exactly zero risk — indistinguishable from a
# genuinely zero-risk person. Use horizon_support() to see which strata/horizons
# are affected.
# ------------------------------------------------------------------------------
.S0_at <- function(base, t, extrapolate = FALSE) {
  if (is.null(base)) return(rep(NA_real_, length(t)))
  idx <- findInterval(t, base$time)
  out <- ifelse(idx == 0L, 1, base$surv[pmax(idx, 1L)])
  if (!extrapolate) {
    tmax <- if (length(base$time)) max(base$time, na.rm = TRUE) else NA_real_
    out[is.finite(t) & t > tmax] <- NA_real_
  }
  as.numeric(out)
}

# ------------------------------------------------------------------------------
# horizon_support() — per-stratum diagnostics for a requested horizon, so an
# unsupported stratum is reported rather than silently returning zero risk.
# Returns a data.frame: stratum, n, n_at_risk_at_end, n_events_by_end,
# censoring_survival, max_baseline_time, supported.
# ------------------------------------------------------------------------------
horizon_support <- function(fit, newdata, horizon, start_col = "START",
                            time_col = "TIME", status_col = "PHENO") {
  base <- .baseline_survival(fit)
  n <- nrow(newdata)
  strat <- rep("(all)", n)
  if (!is.null(base) && isTRUE(base$strata)) {
    inner <- sub("^strata\\((.*)\\)$", "\\1",
                 if (is.null(base$var)) "" else base$var[1])
    if (nzchar(inner) && inner %in% names(newdata))
      strat <- as.character(newdata[[inner]])
  }
  start <- if (start_col %in% names(newdata)) as.numeric(newdata[[start_col]]) else rep(0, n)
  tt    <- if (time_col   %in% names(newdata)) as.numeric(newdata[[time_col]])  else rep(NA_real_, n)
  ev    <- if (status_col %in% names(newdata)) as.numeric(newdata[[status_col]]) else rep(NA_real_, n)

  do.call(rbind, lapply(sort(unique(strat)), function(s) {
    ix <- which(strat == s)
    b  <- if (!is.null(base) && isTRUE(base$strata)) {
            base$bases[[s]] %||% base$bases[[paste0(sub("^strata\\((.*)\\)$", "\\1", base$var[1]), "=", s)]]
          } else if (!is.null(base)) base$base else NULL
    tmax <- if (!is.null(b) && length(b$time)) max(b$time, na.rm = TRUE) else NA_real_
    # follow-up scale: horizon is measured from each person's own entry
    endpoint <- start[ix] + horizon
    cs <- tryCatch({
      kmc <- survival::survfit(survival::Surv(tt[ix], 1 - ev[ix]) ~ 1)
      .S0_at(list(time = kmc$time, surv = kmc$surv), horizon, extrapolate = TRUE)
    }, error = function(e) NA_real_)
    data.frame(
      stratum = s, n = length(ix),
      n_at_risk_at_end = sum(tt[ix] >= horizon, na.rm = TRUE),
      n_events_by_end  = sum(tt[ix] <= horizon & ev[ix] == 1, na.rm = TRUE),
      censoring_survival = as.numeric(cs)[1],
      max_baseline_time  = tmax,
      supported = isTRUE(is.finite(tmax) && all(endpoint <= tmax, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  }))
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# ------------------------------------------------------------------------------
# predict_risk_at_horizon()
#
#   fit      : coxph model (start-stop or single-time)
#   newdata  : data.frame / data.table of individuals
#   horizon  : scalar time, on the FOLLOW-UP scale (years since entry)
#   start_col: column in newdata holding entry time on the model's time scale.
#              Defaults to "START"; if absent or all zero, entry is treated as 0
#              and the result is the ordinary 1 - S0(h)^exp(lp).
#
#   Returns a numeric vector of length nrow(newdata): probability of the event
#   by `horizon` years after entry, conditional on being event-free at entry.
# ------------------------------------------------------------------------------
predict_risk_at_horizon <- function(fit, newdata, horizon,
                                    start_col = "START") {
  n <- nrow(newdata)
  if (is.null(n) || n == 0L) return(numeric(0))

  # -- horizon must be a single finite non-negative number. A vector would be
  #    silently recycled against the per-person START and produce nonsense.
  if (length(horizon) != 1L || !is.finite(horizon) || horizon < 0)
    stop("horizon must be one finite, non-negative number; got length ",
         length(horizon), ".")

  # A zero-length interval carries no risk, whatever the baseline support is.
  # Returned BEFORE the support check below by design.
  if (horizon == 0) return(rep(0, n))

  # -- CENTERING (critical). basehaz(centered=TRUE) is centered at the OVERALL
  #    fit$means, but predict(type="lp") defaults to reference="strata", which
  #    centers within each stratum. Combining them multiplies risk by a
  #    stratum-dependent constant (measured: 2.16x in one stratum, 0.52x in the
  #    other). reference="sample" puts the linear predictor on the same overall
  #    reference as basehaz. For an unstratified model the two are identical, so
  #    this cannot move non-stratified results.
  # Covariate-free models (e.g. ~ strata(X)) have no coefficients; lp is 0 for all.
  # Some survival builds error inside predict() for such a fit, which the tryCatch
  # below would turn into an all-NA risk vector; others return zeros. Branch first
  # so the result is the same either way.
  if (length(stats::coef(fit)) == 0L) {
    lp <- rep(0, n)
  } else {
    lp <- tryCatch(as.numeric(stats::predict(fit, newdata = newdata, type = "lp",
                                             reference = "sample")),
                   error = function(e) NULL)
    if (is.null(lp) || !length(lp)) return(rep(NA_real_, n))
  }

  base <- .baseline_survival(fit)
  if (is.null(base)) return(rep(NA_real_, n))

  # Entry time on the model's own time scale. A missing start_col means the model
  # has no delayed entry (entry at 0); a PRESENT but non-finite value is a data or
  # programming error and must not be silently reinterpreted as entry at time 0.
  start <- if (!is.null(start_col) && start_col %in% names(newdata)) {
    s <- as.numeric(newdata[[start_col]])
    if (any(!is.finite(s)))
      stop("Column '", start_col, "' contains ", sum(!is.finite(s)),
           " missing or non-finite entry time(s); cannot compute conditional risk.")
    s
  } else {
    rep(0, n)
  }

  # S0 at entry and at entry + horizon; the ratio is the conditional baseline.
  if (!isTRUE(base$strata)) {
    S0_entry <- .S0_at(base$base, start)
    S0_exit  <- .S0_at(base$base, start + horizon)
  } else {
    # Stratified model: use each individual's stratum baseline. Match basehaz
    # labels like "T1_STATUS=Incident" against the newdata factor column.
    # Only a SINGLE strata() term is supported: with two the basehaz labels are
    # interaction cells and base$var[1] would silently select the wrong baseline.
    if (length(base$var) != 1L)
      stop("Exactly one strata() term is supported; this model has ",
           length(base$var), " (", paste(base$var, collapse = ", "), ").")
    inner <- sub("^strata\\((.*)\\)$", "\\1", if (is.null(base$var)) "" else base$var[1])
    if (!nzchar(inner) || !inner %in% names(newdata))
      stop("Could not identify the model strata variable ('", inner,
           "') in newdata; cannot select the correct baseline hazard.")
    labels  <- names(base$bases)
    row_lab <- paste0(inner, "=", as.character(newdata[[inner]]))
    # Fall back: if labels don't carry "var=", match on the bare level.
    if (!any(row_lab %in% labels)) row_lab <- as.character(newdata[[inner]])
    unmatched <- !row_lab %in% labels
    if (any(unmatched))
      stop("Unmatched strata level(s) in newdata: ",
           paste(unique(row_lab[unmatched]), collapse = ", "),
           ". Baseline hazards are available for: ", paste(labels, collapse = ", "), ".")
    S0_entry <- rep(NA_real_, n); S0_exit <- rep(NA_real_, n)
    for (lb in unique(row_lab)) {
      b <- base$bases[[lb]]
      idx <- which(row_lab == lb)
      if (is.null(b)) next               # stays NA
      S0_entry[idx] <- .S0_at(b, start[idx])
      S0_exit[idx]  <- .S0_at(b, start[idx] + horizon)
    }
  }

  # NA denominator => entry outside baseline support; NA exit => horizon beyond
  # support (see .S0_at). Both propagate as NA rather than a fabricated 0 risk.
  ratio <- ifelse(is.finite(S0_entry) & S0_entry > 0 & is.finite(S0_exit),
                  S0_exit / S0_entry, NA_real_)
  ratio <- pmin(pmax(ratio, 0), 1)

  # risk = 1 - ratio^exp(lp), computed as -expm1(exp(lp)*log(ratio)) to avoid
  # cancellation when ratio is very close to 1 and overflow for large lp.
  risk <- rep(NA_real_, length(ratio))
  ok <- is.finite(ratio) & ratio > 0 & ratio < 1 & is.finite(lp)
  risk[ok] <- -expm1(exp(lp[ok]) * log(ratio[ok]))
  risk[is.finite(ratio) & ratio == 1] <- 0
  risk[is.finite(ratio) & ratio == 0 & is.finite(lp)] <- 1
  as.numeric(pmin(pmax(risk, 0), 1))
}

# ------------------------------------------------------------------------------
# predict_survival_at_horizon() - convenience complement, S(h) = 1 - risk(h).
# ------------------------------------------------------------------------------
predict_survival_at_horizon <- function(fit, newdata, horizon, start_col = "START") {
  1 - predict_risk_at_horizon(fit, newdata, horizon, start_col = start_col)
}

# ------------------------------------------------------------------------------
# ipcw_weights() — inverse-probability-of-censoring weights for reducing a
# censored survival outcome to a binary outcome AT A HORIZON.
#
# Evaluating "did the event happen by time h?" with the eventual status is wrong
# twice over: someone whose event occurs after h is counted as an event at h, and
# someone censored before h is counted as a non-event. IPCW fixes this by
# restricting to individuals whose status at h is KNOWN and up-weighting them by
# the probability of remaining uncensored:
#
#   event by h            -> D = 1, w = 1 / G(T_i-)
#   event-free through h  -> D = 0, w = 1 / G(h)
#   censored before h     -> status unknown, w = 0 (dropped)
#
# G is the Kaplan-Meier estimate of the CENSORING distribution (the "reverse KM"),
# the same quantity 03_discrimination_metrics.R computes for its G-hat guards.
#
# Returns list(D, w, keep, n_dropped, G_at_horizon).
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# suppress_small_cells() — consortium disclosure control for published subgroup
# tables. A row is suppressed when ANY of its count columns is a small non-zero
# integer (1..threshold-1): both the count columns AND the named estimate columns
# are set to NA, and a minimum_cell_count_pass flag records the outcome. Zero and
# NA counts are not themselves disclosive and pass. Modifies and returns a
# data.table. Defined once here (sourced by 05/06) and inlined identically in
# 01/02, which do not source this file.
# ------------------------------------------------------------------------------
suppress_small_cells <- function(x, count_cols, threshold = 10L, estimate_cols = character()) {
  count_cols    <- intersect(count_cols, names(x))
  estimate_cols <- intersect(estimate_cols, names(x))
  if (!length(count_cols)) { x[, minimum_cell_count_pass := TRUE]; return(x[]) }
  small <- Reduce(`|`, lapply(count_cols, function(cc) {
    v <- suppressWarnings(as.numeric(x[[cc]])); is.finite(v) & v > 0 & v < threshold
  }))
  small[is.na(small)] <- FALSE
  x[, minimum_cell_count_pass := !small]
  if (any(small)) for (cc in unique(c(count_cols, estimate_cols)))
    x[small, (cc) := NA]
  x[]
}

ipcw_weights <- function(time, status, horizon) {
  time <- as.numeric(time); status <- as.numeric(status)
  n <- length(time)
  if (length(horizon) != 1L || !is.finite(horizon) || horizon <= 0)
    stop("horizon must be one finite, positive number.")

  # KM of the censoring distribution: censoring is the "event" here.
  kmc <- survival::survfit(survival::Surv(time, 1 - status) ~ 1)
  Gstep <- list(time = kmc$time, surv = kmc$surv)
  # G just BEFORE t (left limit): use the last censoring time strictly < t.
  G_minus <- function(t) {
    idx <- findInterval(t - .Machine$double.eps^0.5, Gstep$time)
    as.numeric(ifelse(idx == 0L, 1, Gstep$surv[pmax(idx, 1L)]))
  }
  G_h <- G_minus(horizon)

  event_by_h <- status == 1 & time <= horizon
  free_at_h  <- time > horizon
  cens_pre_h <- status == 0 & time <= horizon

  D <- ifelse(event_by_h, 1, 0)
  w <- rep(0, n)
  gi <- G_minus(time)
  w[event_by_h] <- ifelse(gi[event_by_h] > 0, 1 / gi[event_by_h], 0)
  w[free_at_h]  <- if (is.finite(G_h) && G_h > 0) 1 / G_h else 0
  # censored before the horizon keep weight 0 — status at h is unknowable

  list(D = D, w = w, keep = !cens_pre_h & is.finite(w) & w > 0,
       n_dropped = sum(cens_pre_h), G_at_horizon = G_h)
}

# ------------------------------------------------------------------------------
# source_prs_risk_utils() is not defined here; callers locate this file with:
#
#   .this_dir <- tryCatch(dirname(normalizePath(sub("^--file=", "",
#                  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))),
#                  error = function(e) ".")
#   source(file.path(.this_dir, "prs_risk_utils.R"))
#
# so the helper is found next to the calling script regardless of working dir.
# ------------------------------------------------------------------------------
