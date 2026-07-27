#!/usr/bin/env bash
# =============================================================================
# smoke_test.sh — synthetic end-to-end test for the PRS-metrics pipeline.
#
# Generates a small synthetic dataset (no real individuals), runs one
# GWAS-aligned (Layer 1) and one prospective (Layer 2) job through the driver,
# and asserts:
#   * both jobs print "RUN COMPLETE: yes";
#   * canonical output files exist with the expected prefix;
#   * the three key contrasts (M4:M1, M5:M2, M7:M6) are labelled in the
#     discrimination contrast table;
#   * run_status has M4/M5/M7 *_available = TRUE;
#   * disclosure sanitisation stripped the count columns from published
#     horizon_support.csv;
#   * no files were withheld as "unexpected".
#
# Exits 0 on success, non-zero on the first failed assertion. This is a
# functional smoke test, NOT a scientific validation.
#
# Usage:  bash smoke_test.sh [work_dir]
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK="${1:-$(mktemp -d 2>/dev/null || echo /tmp/prs_smoke_$$)}"
DATA="$WORK/data"; OUT1="$WORK/gwas_aligned"; OUT2="$WORK/prospective"
mkdir -p "$DATA"
echo "Smoke test working directory: $WORK"

fail() { echo "  ✗ SMOKE TEST FAILED: $*" >&2; exit 1; }
pass() { echo "  ✓ $*"; }

# ---- 1. Generate synthetic data (deterministic; ~3000 people, ~600 events) ---
Rscript --vanilla - "$DATA" <<'GEN'
set.seed(42)
args <- commandArgs(trailingOnly = TRUE); d <- args[1]
N <- 3000L
iid <- sprintf("S%05d", seq_len(N))
onset <- rnorm(N); prog <- rnorm(N); outc <- rnorm(N)
diagAge <- runif(N, 40, 70)
# weak signal so coefficients are estimable (not for inference)
lp  <- 0.25*prog + 0.15*outc + 0.10*onset
event <- rbinom(N, 1, plogis(-0.4 + lp))
gap_event  <- pmax(0.2, rexp(N, rate = 1/6))       # T1->T2 time for events
gap_censor <- pmax(0.2, runif(N, 1, 22))           # follow-up for non-events
gap <- ifelse(event == 1L, gap_event, gap_censor)
ageExit <- diagAge + gap
sex <- rbinom(N, 1, 0.5)
pcs <- matrix(rnorm(N*5), N, 5); colnames(pcs) <- paste0("PC", 1:5)
birthyear <- round(runif(N, 1940, 1975))
write.table(data.frame(IID = iid, SCORE1_SUM = onset),
            file.path(d, "onset.sscore"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(IID = iid, SCORE1_SUM = prog),
            file.path(d, "prog.sscore"),  sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(IID = iid, SCORE1_SUM = outc),
            file.path(d, "outcome.sscore"), sep = "\t", row.names = FALSE, quote = FALSE)
pheno <- data.frame(IID = iid, secondEvent = event, diagAge = round(diagAge, 3),
                    ageExit = round(ageExit, 3), secondTime = round(gap, 3),
                    sex = sex, round(pcs, 4), birthyear = birthyear)
write.table(pheno, file.path(d, "pheno.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
# recruitment age -> a mix of prevalent (T1 before recruit) and incident (T1 after)
recruitAge <- diagAge + runif(N, -8, 8)
write.table(data.frame(IID = iid, AgeOfConsent = round(recruitAge, 3)),
            file.path(d, "recruit.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("Synthetic events:", sum(event), "of", N, "\n")
GEN
pass "generated synthetic data"

COMMON=( --trait=SYN_to_TEST --script-dir="$SCRIPT_DIR"
         --onset-prs-file="$DATA/onset.sscore" --progression-prs-file="$DATA/prog.sscore"
         --outcome-prs-file="$DATA/outcome.sscore" --pheno-file="$DATA/pheno.txt"
         --cohort=SMOKE --ancestry=EUR --age-t1-col=diagAge --age-exit-col=ageExit
         --recruit-file="$DATA/recruit.txt" --recruit-age-col=AgeOfConsent
         --prs-method=PRScs --discovery-ancestry=EUR --time-points="1,2,5"
         --covariates="sex,PC1,PC2,PC3,PC4,PC5,birthyear" )

# ---- 2. Layer 1 (GWAS-aligned, no lag sweep for speed) ----------------------
echo "Running Layer 1 (gwas_aligned)…"
bash "$SCRIPT_DIR/00_run_pipeline_gbmi.sh" "${COMMON[@]}" --out-dir="$OUT1" --lag-days=0 \
  > "$WORK/l1.log" 2>&1 || { tail -30 "$WORK/l1.log"; fail "Layer 1 run exited non-zero"; }
grep -q "RUN COMPLETE: yes" "$WORK/l1.log" || { tail -20 "$WORK/l1.log"; fail "Layer 1 not RUN COMPLETE"; }
pass "Layer 1 RUN COMPLETE: yes"

# ---- 3. Layer 2 (prospective) -----------------------------------------------
echo "Running Layer 2 (prospective)…"
bash "$SCRIPT_DIR/00_run_pipeline_gbmi.sh" "${COMMON[@]}" --out-dir="$OUT2" \
  --analysis-mode=prospective --lag-days=0 --run-tests=0 \
  > "$WORK/l2.log" 2>&1 || { tail -30 "$WORK/l2.log"; fail "Layer 2 run exited non-zero"; }
grep -q "RUN COMPLETE: yes" "$WORK/l2.log" || { tail -20 "$WORK/l2.log"; fail "Layer 2 not RUN COMPLETE"; }
pass "Layer 2 RUN COMPLETE: yes"

# ---- 4. Assertions on the Layer-1 final/ ------------------------------------
FINAL="$OUT1/final"; PFX="SMOKE_EUR_PRScs_SYN_to_TEST"
for f in run_status run_configuration session_manifest model_lrt_ladder \
         model_comparison rr_auc_contrasts final_analysis_summary horizon_support; do
  ls "$FINAL/${PFX}"*"_${f}.csv" >/dev/null 2>&1 || fail "missing canonical output: *_${f}.csv"
done
pass "canonical outputs present"

# no unexpected withheld files
if grep -q "matched no allowlist entry" "$WORK/l1.log"; then
  grep -A20 "matched no allowlist entry" "$WORK/l1.log"; fail "unexpected withheld file(s) in Layer 1"
fi
pass "no unexpected withheld files"

Rscript --vanilla - "$FINAL" "$PFX" <<'CHK'
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE); final <- a[1]; pfx <- a[2]
g <- function(suffix) list.files(final, pattern = paste0(pfx, ".*", suffix, "$"), full.names = TRUE)[1]
rs <- fread(g("run_status.csv"))
for (col in c("M4_vs_M1_available","M5_vs_M2_available","M7_vs_M6_available","run_complete"))
  if (!isTRUE(as.logical(rs[[col]][1]))) stop("run_status ", col, " is not TRUE")
ct <- fread(g("rr_auc_contrasts.csv"))
have <- unique(ct$contrast_label[!is.na(ct$contrast_label)])
for (lb in c("M4:M1","M5:M2","M7:M6"))
  if (!(lb %in% have)) stop("contrast_label ", lb, " missing from rr_auc_contrasts.csv")
if (!("report_tier" %in% names(ct))) stop("rr_auc_contrasts.csv lacks report_tier column")
hs <- fread(g("horizon_support.csv"))
leaked <- intersect(c("n","n_at_risk_at_end","n_events_by_end"), names(hs))
if (length(leaked)) stop("horizon_support.csv still exposes count column(s): ", paste(leaked, collapse=", "))
# published time audit must carry status but NOT the raw count columns
au <- fread(g("requested_timepoints_audit.csv"))
aul <- intersect(c("n_events_by_t","n_at_risk"), names(au))
if (length(aul)) stop("requested_timepoints_audit.csv still exposes count column(s): ", paste(aul, collapse=", "))
# unit tests were run (Layer 1 keeps them on) and recorded PASS in the manifest
mf <- fread(g("session_manifest.csv"), header = FALSE)
utp <- mf[grepl("unit_tests_passed", V1)]
if (nrow(utp) && !any(grepl("TRUE", unlist(utp)))) stop("session_manifest unit_tests_passed is not TRUE")
cat("R assertions OK\n")
CHK
[[ $? -eq 0 ]] || fail "R assertions failed"
pass "run_status co-primary flags TRUE; M4:M1/M5:M2/M7:M6 labelled; horizon_support + time-audit sanitised; unit_tests_passed"

# ---- 5. Layer-2 canonical outputs + no unexpected withheld ------------------
FINAL2="$OUT2/final"
for f in run_status model_lrt_ladder rr_auc_contrasts final_analysis_summary timing_categories; do
  ls "$FINAL2/${PFX}"*"_${f}.csv" >/dev/null 2>&1 || fail "Layer 2 missing canonical output: *_${f}.csv"
done
grep -q "matched no allowlist entry" "$WORK/l2.log" && { grep -A20 "matched no allowlist entry" "$WORK/l2.log"; fail "unexpected withheld file(s) in Layer 2"; }
pass "Layer 2 canonical outputs present; no unexpected withheld files"

# ---- 6. Downstream-only rerun preserves the original fit provenance ----------
echo "Running downstream-only rerun (--run-01=0)…"
sha() { shasum -a 256 "$1" 2>/dev/null | awk '{print $1}'; }
man="$OUT1/intermediate/${PFX}_session_manifest.csv"; cfg="$OUT1/intermediate/${PFX}_run_configuration.csv"
man_before="$(sha "$man")"; cfg_before="$(sha "$cfg")"
bash "$SCRIPT_DIR/00_run_pipeline_gbmi.sh" --trait=SYN_to_TEST --out-dir="$OUT1" --script-dir="$SCRIPT_DIR" \
  --cohort=SMOKE --ancestry=EUR --prs-method=PRScs --run-01=0 --run-tests=0 \
  > "$WORK/rerun.log" 2>&1 || { tail -30 "$WORK/rerun.log"; fail "downstream rerun exited non-zero"; }
[[ "$(sha "$man")" == "$man_before" ]] || fail "downstream rerun OVERWROTE session_manifest.csv"
[[ "$(sha "$cfg")" == "$cfg_before" ]] || fail "downstream rerun OVERWROTE run_configuration.csv"
ls "$OUT1/intermediate/${PFX}_downstream_rerun_manifest.csv" >/dev/null 2>&1 || fail "downstream_rerun_manifest.csv not written"
pass "downstream rerun preserved session_manifest + run_configuration; wrote downstream_rerun_manifest"

echo ""
echo "✅ SMOKE TEST PASSED  (artifacts in $WORK)"
