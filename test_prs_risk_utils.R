#!/usr/bin/env Rscript
# ==============================================================================
# test_prs_risk_utils.R - unit tests for the shared risk helper
# ==============================================================================
# Run standalone:  Rscript test_prs_risk_utils.R
#
# The oracle is survival's own native conditional prediction,
#   predict(fit, newdata = <START, START+h, status=0>, type = "survival"),
# which handles centering, strata, splines and delayed entry internally. Where
# the oracle can express a case, the helper must agree to < 1e-8.
#
# The oracle CANNOT express two cases (verified): horizon = 0 produces an invalid
# zero-length Surv interval, and a covariate-free model errors inside predict().
# Those are asserted on the helper alone.
#
# WHY THIS FILE EXISTS: the helper previously combined basehaz(centered = TRUE)
# (overall fit$means reference) with predict(type = "lp") (per-stratum reference
# by default). On stratified models that inflated risk in one stratum and
# shrank it in the other. Test 4 is the regression guard for that bug.
# ==============================================================================

suppressPackageStartupMessages({
  library(survival); library(splines)
})

# Locate prs_risk_utils.R robustly. Prefer PRS_SCRIPT_DIR (exported by the driver) so
# we never parse --file, whose spaces some macOS/OneDrive setups encode as ~+~ (which
# breaks normalizePath()). Fall back to --file (with ~+~ decoded), then to getwd().
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
if (!file.exists(.utils))
  stop("Required helper not found: ", .utils, "\nResolved script directory: ", .this_dir)
source(.utils)

TOL <- 1e-8
pass <- 0L; fail <- 0L
ok <- function(label, cond, extra = "") {
  if (isTRUE(cond)) { pass <<- pass + 1L; cat(sprintf("  PASS  %s %s\n", label, extra)) }
  else              { fail <<- fail + 1L; cat(sprintf("  FAIL  %s %s\n", label, extra)) }
}

# Oracle: native conditional survival over [START, START + h].
native_risk <- function(fit, dd, h, start_col = "START",
                        stop_col = "STOP", status_col = "PHENO") {
  nd <- as.data.frame(dd)
  nd[[stop_col]]   <- nd[[start_col]] + h
  nd[[status_col]] <- 0L
  1 - as.numeric(predict(fit, newdata = nd, type = "survival"))
}
agree <- function(fit, dd, h, label) {
  a <- predict_risk_at_horizon(fit, dd, h)
  b <- native_risk(fit, dd, h)
  keep <- is.finite(a) & is.finite(b)
  d <- if (any(keep)) max(abs(a[keep] - b[keep])) else NA_real_
  ok(label, is.finite(d) && d < TOL, sprintf("(max|diff| = %.3g, n = %d)", d, sum(keep)))
}

# ------------------------------------------------------------------ fixtures --
mk <- function(seed = 11, n = 3000, delayed = c("none","constant","varying"),
               strat = FALSE) {
  delayed <- match.arg(delayed)
  set.seed(seed)
  g   <- factor(sample(c("Prevalent","Incident"), n, TRUE))
  # deliberately divergent per-stratum covariate means: this is what makes the
  # strata-vs-sample centering bug visible
  x   <- rnorm(n) + if (strat) ifelse(g == "Prevalent", 1.5, -1.5) else 0
  age <- rnorm(n, 60, 8)
  START <- switch(delayed,
                  none     = rep(0, n),
                  constant = rep(0.25, n),
                  varying  = runif(n, 0, 3))
  rate <- 0.05 * exp(0.5 * x) * if (strat) ifelse(g == "Prevalent", 2, 1) else 1
  d <- data.frame(START = START, STOP = START + rexp(n, rate),
                  PHENO = rbinom(n, 1, 0.6), x = x, AGE_INDEX = age, T1_STATUS = g)
  d$TIME <- d$STOP - d$START
  d[d$STOP > d$START, ]
}
fitm <- function(d, form) coxph(form, data = d, x = TRUE, y = TRUE, model = TRUE)

cat("\n=== prs_risk_utils.R unit tests ===\n\n")

# 1-3. non-stratified, three delayed-entry regimes
cat("Agreement with native oracle:\n")
for (dl in c("none","constant","varying")) {
  d <- mk(delayed = dl)
  f <- fitm(d, Surv(START, STOP, PHENO) ~ x)
  agree(f, d, 5, sprintf("1-3. non-stratified, delayed entry = %-8s", dl))
}

# 4. THE REGRESSION GUARD: stratified, both statuses present
d4 <- mk(delayed = "varying", strat = TRUE)
f4 <- fitm(d4, Surv(START, STOP, PHENO) ~ x + strata(T1_STATUS))
agree(f4, d4, 5, "4.   strata(T1_STATUS), both statuses  ")
# and per-stratum ratio must be 1 (the bug showed as 2.16x / 0.52x)
r_h <- predict_risk_at_horizon(f4, d4, 5); r_n <- native_risk(f4, d4, 5)
kp <- is.finite(r_h) & is.finite(r_n) & r_n > 0
rat <- tapply((r_h / r_n)[kp], d4$T1_STATUS[kp], median)
ok("4b.  per-stratum helper/native ratio == 1",
   all(abs(rat - 1) < 1e-6), sprintf("(%s)", paste(names(rat), round(rat, 6), collapse = ", ")))

# 5. spline term
d5 <- mk(delayed = "varying", strat = TRUE)
f5 <- fitm(d5, Surv(START, STOP, PHENO) ~ x + ns(AGE_INDEX, df = 3) + strata(T1_STATUS))
agree(f5, d5, 5, "5.   ns(AGE_INDEX, df=3) + strata     ")

# 12. equivalence also holds at a second horizon
agree(f5, d5, 10, "12.  same model, horizon = 10        ")

cat("\nCases the native oracle cannot express (helper-only assertions):\n")

# 7. horizon = 0 must be exactly zero risk
r0 <- predict_risk_at_horizon(f5, d5, 0)
ok("7.   horizon = 0 -> exactly 0 risk", all(r0 == 0), sprintf("(max = %g)", max(abs(r0))))
nat0 <- suppressWarnings(tryCatch(native_risk(f5, d5, 0), error = function(e) NA))
ok("7b.  (oracle indeed cannot: NA/error)", all(is.na(nat0)) || inherits(nat0, "try-error"))

# 6. covariate-free model (strata only)
d6 <- mk(delayed = "varying", strat = TRUE)
f6 <- fitm(d6, Surv(START, STOP, PHENO) ~ strata(T1_STATUS))
r6 <- tryCatch(predict_risk_at_horizon(f6, d6, 5), error = function(e) e)
ok("6.   covariate-free model returns risk",
   is.numeric(r6) && any(is.finite(r6)) && all(r6 >= 0 & r6 <= 1, na.rm = TRUE),
   if (is.numeric(r6)) sprintf("(median = %.4f)", median(r6, na.rm = TRUE)) else "(errored)")

cat("\nValidation and support semantics:\n")

# 8. missing START must stop(), not silently become 0
d8 <- d5; d8$START[1:5] <- NA
ok("8.   non-finite START -> stop()",
   inherits(tryCatch(predict_risk_at_horizon(f5, d8, 5), error = function(e) e), "error"))

# 9. unmatched stratum level must stop()
d9 <- d5; levels(d9$T1_STATUS) <- c(levels(d5$T1_STATUS), "Unseen")
d9$T1_STATUS[1:10] <- "Unseen"
ok("9.   unmatched strata level -> stop()",
   inherits(tryCatch(predict_risk_at_horizon(f5, d9, 5), error = function(e) e), "error"))

# 10. horizon beyond baseline support -> NA, never 0
tmax <- max(basehaz(f5, centered = TRUE)$time)
r10 <- predict_risk_at_horizon(f5, d5, ceiling(tmax) + 50)
ok("10.  unsupported horizon -> NA (not 0)",
   all(is.na(r10)), sprintf("(NA: %d/%d, zeros: %d)", sum(is.na(r10)), length(r10),
                            sum(r10 == 0, na.rm = TRUE)))
hs <- horizon_support(f5, d5, ceiling(tmax) + 50)
ok("10b. horizon_support flags unsupported", all(!hs$supported))
hs_ok <- horizon_support(f5, d5, 5)
ok("10c. horizon_support flags supported", all(hs_ok$supported),
   sprintf("(strata: %s)", paste(hs_ok$stratum, collapse = ", ")))

# 11. vector horizon must stop()
ok("11.  vector horizon -> stop()",
   inherits(tryCatch(predict_risk_at_horizon(f5, d5, c(5, 10)), error = function(e) e), "error"))

# 13. predictions survive a saveRDS/readRDS round trip
tmp <- tempfile(fileext = ".rds"); saveRDS(f5, tmp); f5b <- readRDS(tmp); unlink(tmp)
ok("13.  stable across saveRDS/readRDS",
   isTRUE(all.equal(predict_risk_at_horizon(f5, d5, 5),
                    predict_risk_at_horizon(f5b, d5, 5))))

cat("\nIPCW weights:\n")
# with no censoring before the horizon, IPCW must collapse to the plain counts
d14 <- mk(delayed = "none"); h <- 5
d14$PHENO[d14$TIME <= h & d14$PHENO == 0] <- 1        # remove pre-horizon censoring
w <- ipcw_weights(d14$TIME, d14$PHENO, h)
ok("14.  no pre-horizon censoring -> none dropped", w$n_dropped == 0,
   sprintf("(dropped = %d)", w$n_dropped))
ok("14b. weighted event count ~ naive count",
   abs(sum(w$D * w$w) - sum(d14$TIME <= h & d14$PHENO == 1)) /
     max(1, sum(d14$TIME <= h & d14$PHENO == 1)) < 0.05)
d15 <- mk(delayed = "none")
w15 <- ipcw_weights(d15$TIME, d15$PHENO, h)
ok("15.  censored-before-horizon are dropped", w15$n_dropped > 0 &&
     all(w15$w[d15$PHENO == 0 & d15$TIME <= h] == 0),
   sprintf("(dropped = %d, G(h) = %.4f)", w15$n_dropped, w15$G_at_horizon))

cat(sprintf("\n=== %d passed, %d failed ===\n", pass, fail))
if (fail > 0) quit(status = 1)
