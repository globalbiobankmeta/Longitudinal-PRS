#!/usr/bin/env bash
# ==============================================================================
# 00_run_pipeline_gbmi.sh
#
# GBMI-aware multi-PRS pipeline. Given a trait name (e.g. T1D_to_CAD or
# PARKINSONtoDEMENTIA), this script:
#   1. Parses the trait name into ONSET, PROGRESSION, and OUTCOME sub-trait names
#      - Onset      = left side of separator  (e.g. PARKINSON)
#      - Progression = the full compound trait (e.g. PARKINSONtoDEMENTIA)
#      - Outcome    = right side of separator (e.g. DEMENTIA)
#   2. Resolves sscore and phenotype files from the standard GBMI directory tree
#   3. Runs 01_fit_models.R in multi-PRS mode — auto-builds all 8 models:
#        Base, Onset, Progression, Outcome,
#        Onset_Progression, Outcome_Progression, Onset_Outcome,
#        Onset_Progression_Outcome
#   4. Chains scripts 02-06 on the output RDS files
#
# Supported trait name formats:
#   ONSET_to_OUTCOME   e.g. T1D_to_CAD, T2D_to_CAD
#   ONSETtoOUTCOME     e.g. PARKINSONtoDEMENTIA, CADtoHEARTFAIL
#
# Inputs are supplied explicitly (required): --onset/progression/outcome-prs-file,
# --pheno-file and --out-dir.
# ==============================================================================
set -euo pipefail

# ============================== Defaults ==============================
TRAIT="${TRAIT:-}"              # required; e.g. T1D_to_CAD
OUT_DIR="${OUT_DIR:-}"          # required; each run gets its own output directory
SCRIPT_DIR="${SCRIPT_DIR:-$(dirname "$(realpath "$0")")}"

PIPELINE_VERSION="1.0.1-rc3"
COHORT="${COHORT:-MGBB}"
ANCESTRY="${ANCESTRY:-EUR}"
PRS_METHOD="${PRS_METHOD:-PRScs}"                # goes into the output prefix (P+T vs PRS-CS)

# Optional provenance recorded verbatim in the session manifest (no analysis effect)
DISCOVERY_ANCESTRY="${DISCOVERY_ANCESTRY:-}"
RUN_TESTS="${RUN_TESTS:-1}"                      # run helper unit tests in preflight

STATUS_COL="${STATUS_COL:-secondEvent}"
SEX_COL="${SEX_COL:-sex}"

# ---- Time scale (age-derived, preferred) --------------------------------------
# STOP = age_at_exit - age_at_T1 is derived in 01 from these two AGE columns.
# age_at_T1 also serves as the ns() age covariate (one input, two uses).
AGE_T1_COL="${AGE_T1_COL:-diagAge}"        # age at T1 diagnosis (index age + ns() covariate)
# age at T2 / last follow-up. The time scale must be stated explicitly: pass EITHER
# --age-exit-col (recommended: STOP = age_exit - age_t1 derived in 01) OR --time-exit
# (a precomputed time-since-T1 column). There is no implicit default — the driver
# stops otherwise, so two biobanks cannot silently analyse different time definitions.
AGE_EXIT_COL="${AGE_EXIT_COL:-}"
AGE_SPLINE_DF="${AGE_SPLINE_DF:-3}"
TIME_ENTRY="${TIME_ENTRY:-}"
TIME_EXIT="${TIME_EXIT:-secondTime}"        # precomputed time-since-T1 col; used only when --time-exit is passed explicitly
TIME_EXIT_EXPLICIT=0                        # set to 1 when --time-exit is passed on the CLI
COV_FILE="${COV_FILE:-}"                    # optional separate covariate file
TIME_SUMMARY_COLS="${TIME_SUMMARY_COLS:-diagAge}"

# age at T1 is entered as ns(); do NOT also list it (or age.x) among plain covariates
COVARIATES="${COVARIATES:-sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,birthyear}"

# ---- Analysis layers ----------------------------------------------------------
ANALYSIS_MODE="${ANALYSIS_MODE:-gwas_aligned}"   # gwas_aligned (Layer 1) | prospective (Layer 2)
RECRUIT_FILE="${RECRUIT_FILE:-}"                 # Layer 2 re-indexing; Layer 1 timing QC
RECRUIT_AGE_COL="${RECRUIT_AGE_COL:-age_at_recruitment}"
POP_FILE="${POP_FILE:-}"                         # optional keep-list (FID IID); restrict the run to these
                                                 # IIDs (e.g. unrelated individuals, or one target ancestry)
CALC_RMST="${CALC_RMST:-1}"                      # restricted mean progression-free time (stage 02)
CALC_UNO="${CALC_UNO:-0}"                        # Uno's C in stage 01 (off; pooled Uno skipped in prospective)
RMST_N_BOOT="${RMST_N_BOOT:-200}"                # bootstrap reps for RMST CIs/p (0 = no CIs)
CLINICAL_COVARIATES="${CLINICAL_COVARIATES:-}"   # enables the Clinical M0-M7 ladder in stage 01
ENABLE_CLINICAL_UTILITY="${ENABLE_CLINICAL_UTILITY:-0}"  # stage 06 net benefit / PAF / NNS / lifetime risk
BOOTSTRAP_RISK_DIFFERENCES="${BOOTSTRAP_RISK_DIFFERENCES:-0}"  # stage 06 resampling CIs (slow; off by default)
MIN_CELL_COUNT="${MIN_CELL_COUNT:-10}"           # stage 06 disclosure: suppress groups below this at-risk
RECLASS_PAIR="${RECLASS_PAIR:-M5:M2}"            # stage 04 comparison:reference (manuscript primary)
RISK_MODEL="${RISK_MODEL:-M5}"                   # stage 06 model for predicted risk / decision curves
SCORE_ROLE="${SCORE_ROLE:-all}"                  # which PRS drives 02/06 displays
ROLE_QUANTILES="${ROLE_QUANTILES:-Q5,Q10}"       # quantile schemes per role
INCIDENT_LAG_DAYS="${INCIDENT_LAG_DAYS:-0}"
PROSPECTIVE_T1="${PROSPECTIVE_T1:-both}"   # prospective: include both T1 states, or incident T1 only
LAG_DAYS="${LAG_DAYS:-0,30,90,365}"              # Layer 3 sweep; lag 0 = primary
MIN_EVENTS_TOTAL="${MIN_EVENTS_TOTAL:-50}"

TIME_POINTS="${TIME_POINTS:-1,2,5,10,15,20,25,30,35}"
RISK_THRESHOLDS="${RISK_THRESHOLDS:-0.05,0.075,0.1,0.2}"
USE_SUPPORTED_TIMES="${USE_SUPPORTED_TIMES:-1}"
PRS_QUANTILES="${PRS_QUANTILES:-5,10}"
PRS_EXTREMES="${PRS_EXTREMES:-1,2,5,10,20}"

# Score column name inside .sscore files (SCORE1_SUM is PRScs default)
PRS_SCORE_COL="${PRS_SCORE_COL:-SCORE1_SUM}"

# The manuscript analysis always requires the full three-PRS set.
PROGRESSION_PRS_FILE_OVERRIDE=""

# Optional explicit model spec — if set, overrides auto-combination
# Format: "Name=Label1+Label2;Name2=Label3"
MODELS_SPEC="${MODELS_SPEC:-}"

# Stage toggles (1=run, 0=skip)
RUN_01="${RUN_01:-1}"
RUN_02="${RUN_02:-1}"
RUN_03="${RUN_03:-1}"
# Stage 04 (IDI/NRI) is OFF by default. It is now censoring-aware (IPCW at the
# horizon), but reclassification is not needed for the manuscript's primary
# claims and its interpretation is contested; enable deliberately with --run-04=1.
RUN_04="${RUN_04:-0}"
RUN_05="${RUN_05:-1}"
RUN_06="${RUN_06:-1}"

# ============================== CLI Parser ==============================
usage() {
  cat <<'USAGE'
Usage: ./00_run_pipeline_gbmi.sh --trait=TRAIT --out-dir=DIR \
         --onset-prs-file=FILE --progression-prs-file=FILE --outcome-prs-file=FILE \
         --pheno-file=FILE [options]

Required:
  --trait=STR               Trait name. Supports two formats:
                              ONSET_to_OUTCOME  (e.g. T1D_to_CAD, T2D_to_CAD)
                              ONSETtoOUTCOME   (e.g. PARKINSONtoDEMENTIA)
  --out-dir=DIR             Output directory (created if absent). Give each
                            trajectory / ancestry / PRS-method / analysis-layer its own.
  --onset-prs-file=FILE        Onset (T0->T1) PRS .sscore
  --progression-prs-file=FILE  Progression (T1->T2) PRS .sscore
  --outcome-prs-file=FILE      Outcome (T0->T2) PRS .sscore
  --pheno-file=FILE            Phenotype file
  (All three PRS are mandatory. The four input files are needed only for a fitting
   run; a downstream-only rerun with --run-01=0 works from existing RDS outputs.)

Trait name parts (optional; parsed from TRAIT for the run summary):
  --onset-name=STR          Onset sub-trait   (left of separator, e.g. PARKINSON)
  --outcome-name=STR        Outcome sub-trait (right of separator, e.g. DEMENTIA)

Identity:
  --cohort=STR              (default: MGBB)
  --ancestry=STR            (default: EUR)
  --script-dir=DIR          Directory containing R scripts (default: same dir as this script)

Columns:
  --status-col=NAME         Event indicator column (default: secondEvent)
  --sex-col=NAME            (default: sex)
  --covariates=CSV          Plain covariates; do NOT include age (default: sex,PC1..PC10,birthyear)
  --cov-file=FILE           Covariate file if separate from --pheno-file (default: the phenotype file)
  --min-cell-count=INT      Disclosure: suppress subgroup counts/risks below this at-risk (default: 10)
  --bootstrap-risk-differences  Stage 06: enable the (slow) risk-difference bootstrap CIs (off by default)
  --prs-method=STR          PRS method, embedded in the output prefix (default: PRScs)
  --run-tests=0|1           Run helper unit tests in preflight (default: 1)
  --discovery-ancestry=STR    Provenance recorded in the session manifest (verbatim)
  --prs-score-col=NAME      Column name inside .sscore files (default: SCORE1_SUM)

Time scale (supply exactly ONE, explicitly — there is no implicit default):
  --age-t1-col=NAME         Age at T1 diagnosis; index age AND ns() covariate (default: diagAge)
  --age-exit-col=NAME       Age at T2 (events) / last follow-up (non-events). RECOMMENDED.
                            STOP = age_exit - age_t1 is then derived in-script.
  --age-spline-df=INT       df for ns(age_at_T1); 0 = linear (default: 3)
  --time-exit=NAME          Alternative to --age-exit-col: a precomputed time-since-T1 column.
                            Either --age-exit-col or --time-exit must be passed explicitly;
                            the driver stops rather than assuming a built-in default.
  --time-entry=NAME         Legacy left-truncation column (optional)
  --clinical-covariates=CSV Extra clinical covariates; enables the Clinical M0-M7 ladder

Downstream model / display selection:
  --reclass-pair=A:B        Stage 04 comparison:reference (default: M5:M2, the primary test)
  --risk-model=ID           Stage 06 model for predicted risk / decision curves (default: M5)
  --score-role=STR          Which PRS drives 02/06 displays: all (default) | onset | outcome |
                            progression (comma-separated subset allowed)
  --role-quantiles=CSV      Quantile schemes per role (default: Q5,Q10; 'all' for every scheme)
  --enable-clinical-utility Emit net benefit / PAF / NNS / lifetime risk from stage 06.
                            These are APPARENT (same-sample) estimates; off by default.
  --no-rmst                 Skip restricted mean progression-free time in stage 02
  --calculate-uno           Also compute Uno's C in stage 01 (off by default; the pooled
                            Uno's C is skipped with a warning in prospective/stratified mode)

Analysis layers:
  --analysis-mode=STR       gwas_aligned (Layer 1, default) | prospective (Layer 2)
  --recruit-file=FILE       prospective: IID + age at recruitment
  --recruit-age-col=NAME    (default: age_at_recruitment)
  --pop-file=FILE           Optional keep-list (FID IID); restrict the analysis to these IIDs — e.g.
                            unrelated individuals, or one target ancestry. Default: use the whole input.
  --incident-lag-days=NUM   prospective: lag added to the T1 index for incident T1 (default: 0)
  --prospective-t1=WHICH    prospective: "both" (default) keeps prevalent T1 alongside incident T1,
                            indexing prevalent people at recruitment and separating the two with
                            strata(T1_STATUS) + a T1-duration term. "incident" restricts to T1 that
                            occurs after recruitment: one clean time-since-T1 clock, smaller sample.
  --lag-days=CSV            Layer 3 washout sweep; 0 = primary (default: 0,30,90,365)
  --min-events-total=INT    Hard floor on total T2 events (default: 50)

Models:
  --models-spec=STR         Explicit model combinations (overrides auto-build). MUST still include
                            all seven M1-M7 PRS models — 01 stops if any core model is absent.
                            Format: "Name=Label1+Label2;Name2=Label3"
                            Valid labels: Onset, Progression, Outcome

Analysis:
  --time-points=CSV         (default: 1,2,5,10,15,20,25,30,35)
  --risk-thresholds=CSV     (default: 0.05,0.075,0.1,0.2)
  --prs-quantiles=CSV       (default: 5,10)
  --prs-extremes=CSV        (default: 1,2,5,10,20)
  --use_supported_times=0|1 (default: 1)

Stage toggles (1=run, 0=skip):
  --run-01=0|1  --run-02=0|1  --run-03=0|1
  --run-04=0|1  --run-05=0|1  --run-06=0|1

  -h, --help                Show this help and exit.

Examples:
  # Primary (Layer 1), age-derived time scale, all 8 models M0-M7 + lag sweep
  ./00_run_pipeline_gbmi.sh --trait=T1D_to_CAD --out-dir=/results/T1D_to_CAD \
    --age-t1-col=diagAge --age-exit-col=ageExit

  # Prospective validation (Layer 2) — a time scale must still be given explicitly
  ./00_run_pipeline_gbmi.sh --trait=T1D_to_CAD --out-dir=/results/T1D_to_CAD/prospective \
    --age-t1-col=diagAge --age-exit-col=ageExit \
    --analysis-mode=prospective --recruit-file=/data/recruit.txt --lag-days=0

  # Override the score column (if not SCORE1_SUM)
  ./00_run_pipeline_gbmi.sh --trait=T2D_to_CAD --out-dir=/results --prs-score-col=T2D_to_CAD

  # Skip model fitting, re-run downstream only
  ./00_run_pipeline_gbmi.sh --trait=T1D_to_CAD --out-dir=/results --run-01=0

Environment variable equivalents:
  TRAIT=T1D_to_CAD OUT_DIR=/results ./00_run_pipeline_gbmi.sh
USAGE
}

ONSET_NAME_OVERRIDE=""
OUTCOME_NAME_OVERRIDE=""
ONSET_PRS_FILE_OVERRIDE=""
PROGRESSION_PRS_FILE_OVERRIDE=""
OUTCOME_PRS_FILE_OVERRIDE=""
PHENO_FILE_OVERRIDE=""

while (( $# )); do
  case "$1" in
    -h|--help) usage; exit 0 ;;

    --trait=*)           TRAIT="${1#*=}" ;;
    --trait)             shift; TRAIT="$1" ;;
    --out-dir=*)         OUT_DIR="${1#*=}" ;;
    --out-dir)           shift; OUT_DIR="$1" ;;
    --script-dir=*)      SCRIPT_DIR="${1#*=}" ;;
    --script-dir)        shift; SCRIPT_DIR="$1" ;;


    --onset-name=*)      ONSET_NAME_OVERRIDE="${1#*=}" ;;
    --onset-name)        shift; ONSET_NAME_OVERRIDE="$1" ;;
    --outcome-name=*)    OUTCOME_NAME_OVERRIDE="${1#*=}" ;;
    --outcome-name)      shift; OUTCOME_NAME_OVERRIDE="$1" ;;

    --onset-prs-file=*)        ONSET_PRS_FILE_OVERRIDE="${1#*=}" ;;
    --onset-prs-file)          shift; ONSET_PRS_FILE_OVERRIDE="$1" ;;
    --progression-prs-file=*)  PROGRESSION_PRS_FILE_OVERRIDE="${1#*=}" ;;
    --progression-prs-file)    shift; PROGRESSION_PRS_FILE_OVERRIDE="$1" ;;
    --outcome-prs-file=*)      OUTCOME_PRS_FILE_OVERRIDE="${1#*=}" ;;
    --outcome-prs-file)        shift; OUTCOME_PRS_FILE_OVERRIDE="$1" ;;
    --pheno-file=*)            PHENO_FILE_OVERRIDE="${1#*=}" ;;
    --pheno-file)              shift; PHENO_FILE_OVERRIDE="$1" ;;
    --no-progression)          echo "ERROR: --no-progression is no longer supported: the analysis requires all three PRS (onset, outcome, progression)." >&2; exit 2 ;;

    --age-t1-col=*)      AGE_T1_COL="${1#*=}" ;;
    --age-t1-col)        shift; AGE_T1_COL="$1" ;;
    --age-exit-col=*)    AGE_EXIT_COL="${1#*=}" ;;
    --cov-file=*)        COV_FILE="${1#*=}" ;;
    --cov-file)          shift; COV_FILE="$1" ;;
    --age-exit-col)      shift; AGE_EXIT_COL="$1" ;;
    --age-spline-df=*)   AGE_SPLINE_DF="${1#*=}" ;;
    --age-spline-df)     shift; AGE_SPLINE_DF="$1" ;;

    --analysis-mode=*)   ANALYSIS_MODE="${1#*=}" ;;
    --analysis-mode)     shift; ANALYSIS_MODE="$1" ;;
    --recruit-file=*)    RECRUIT_FILE="${1#*=}" ;;
    --recruit-file)      shift; RECRUIT_FILE="$1" ;;
    --recruit-age-col=*) RECRUIT_AGE_COL="${1#*=}" ;;
    --recruit-age-col)   shift; RECRUIT_AGE_COL="$1" ;;
    --pop-file=*)        POP_FILE="${1#*=}" ;;
    --pop-file)          shift; POP_FILE="$1" ;;
    --incident-lag-days=*) INCIDENT_LAG_DAYS="${1#*=}" ;;
    --incident-lag-days)   shift; INCIDENT_LAG_DAYS="$1" ;;
    --prospective-t1=*)  PROSPECTIVE_T1="${1#*=}" ;;
    --prospective-t1)    shift; PROSPECTIVE_T1="$1" ;;
    --lag-days=*)        LAG_DAYS="${1#*=}" ;;
    --lag-days)          shift; LAG_DAYS="$1" ;;
    --min-events-total=*) MIN_EVENTS_TOTAL="${1#*=}" ;;
    --calculate-rmst=*)  CALC_RMST="${1#*=}" ;;
    --no-rmst)           CALC_RMST=0 ;;
    --calculate-uno=*)   CALC_UNO="${1#*=}" ;;
    --calculate-uno)     CALC_UNO=1 ;;
    --clinical-covariates=*) CLINICAL_COVARIATES="${1#*=}" ;;
    --clinical-covariates)   shift; CLINICAL_COVARIATES="$1" ;;
    --enable-clinical-utility) ENABLE_CLINICAL_UTILITY=1 ;;
    --bootstrap-risk-differences) BOOTSTRAP_RISK_DIFFERENCES=1 ;;
    --min-cell-count=*)  MIN_CELL_COUNT="${1#*=}" ;;
    --prs-method=*)      PRS_METHOD="${1#*=}" ;;
    --run-tests=*)       RUN_TESTS="${1#*=}" ;;
    --discovery-ancestry=*)   DISCOVERY_ANCESTRY="${1#*=}" ;;
    --reclass-pair=*)    RECLASS_PAIR="${1#*=}" ;;
    --risk-model=*)      RISK_MODEL="${1#*=}" ;;
    --score-role=*)      SCORE_ROLE="${1#*=}" ;;
    --role-quantiles=*)  ROLE_QUANTILES="${1#*=}" ;;
    --min-events-total)   shift; MIN_EVENTS_TOTAL="$1" ;;

    --cohort=*)          COHORT="${1#*=}"; COHORT_EXPLICIT=1 ;;
    --cohort)            shift; COHORT="$1"; COHORT_EXPLICIT=1 ;;
    --ancestry=*)        ANCESTRY="${1#*=}"; ANCESTRY_EXPLICIT=1 ;;
    --ancestry)          shift; ANCESTRY="$1"; ANCESTRY_EXPLICIT=1 ;;

    --status-col=*)      STATUS_COL="${1#*=}" ;;
    --status-col)        shift; STATUS_COL="$1" ;;
    --sex-col=*)         SEX_COL="${1#*=}" ;;
    --sex-col)           shift; SEX_COL="$1" ;;
    --covariates=*)      COVARIATES="${1#*=}" ;;
    --covariates)        shift; COVARIATES="$1" ;;
    --prs-score-col=*)   PRS_SCORE_COL="${1#*=}" ;;
    --prs-score-col)     shift; PRS_SCORE_COL="$1" ;;

    --time-entry=*)      TIME_ENTRY="${1#*=}" ;;
    --time-entry)        shift; TIME_ENTRY="$1" ;;
    --time-exit=*)       TIME_EXIT="${1#*=}"; TIME_EXIT_EXPLICIT=1 ;;
    --time-exit)         shift; TIME_EXIT="$1"; TIME_EXIT_EXPLICIT=1 ;;
    --time-summary-cols=*)  TIME_SUMMARY_COLS="${1#*=}" ;;
    --time-summary-cols)    shift; TIME_SUMMARY_COLS="$1" ;;

    --models-spec=*)     MODELS_SPEC="${1#*=}" ;;
    --models-spec)       shift; MODELS_SPEC="$1" ;;

    --time-points=*)     TIME_POINTS="${1#*=}" ;;
    --time-points)       shift; TIME_POINTS="$1" ;;
    --risk-thresholds=*) RISK_THRESHOLDS="${1#*=}" ;;
    --risk-thresholds)   shift; RISK_THRESHOLDS="$1" ;;
    --prs-quantiles=*)   PRS_QUANTILES="${1#*=}" ;;
    --prs-quantiles)     shift; PRS_QUANTILES="$1" ;;
    --prs-extremes=*)    PRS_EXTREMES="${1#*=}" ;;
    --prs-extremes)      shift; PRS_EXTREMES="$1" ;;
    --use_supported_times=*)  USE_SUPPORTED_TIMES="${1#*=}" ;;
    --use_supported_times)    USE_SUPPORTED_TIMES=1 ;;

    --run-01=*)  RUN_01="${1#*=}" ;;
    --run-01)    shift; RUN_01="$1" ;;
    --run-02=*)  RUN_02="${1#*=}" ;;
    --run-02)    shift; RUN_02="$1" ;;
    --run-03=*)  RUN_03="${1#*=}" ;;
    --run-03)    shift; RUN_03="$1" ;;
    --run-04=*)  RUN_04="${1#*=}" ;;
    --run-04)    shift; RUN_04="$1" ;;
    --run-05=*)  RUN_05="${1#*=}" ;;
    --run-05)    shift; RUN_05="$1" ;;
    --run-06=*)  RUN_06="${1#*=}" ;;
    --run-06)    shift; RUN_06="$1" ;;

    --*) echo "Unknown option: $1" >&2; usage; exit 2 ;;
    *)   echo "Unexpected argument: $1" >&2; usage; exit 2 ;;
  esac
  shift || true
done

# ============================== Helpers ==============================
die()  { echo -e "ERROR: $*" >&2; exit 1; }
info() { echo -e "\n==> $*\n"; }
uset() { [[ "$USE_SUPPORTED_TIMES" == "1" ]] && echo "--use_supported_times" || true; }
opt()  { local name="$1" value="${2:-}"; [[ -n "$value" ]] && printf -- "%s %q " "$name" "$value" || true; }

# ---- run_stage: invoke a stage and, on failure, SHOW WHY -----------------------
# Every stage redirects its own stderr into a log with sink(type="message"), so a
# fatal R error never reaches this console. Combined with `set -e` that made any
# stage failure silent: the run just stopped mid-output with no message and no
# hint that a log existed. Print the tail of the right log before dying.
#
# Two naming conventions are in use and both must be searched: 01 writes
# <prefix>fit_models.log and 03 writes <prefix>discrimination_metrics.log, while
# 02/04/05/06 all write <prefix>log.txt.
run_stage() {
  local label="$1" script="$2"; shift 2
  # Capture the stage's REAL exit code. `if Rscript ...; then return 0; fi` would
  # leave $? = 0 on failure (a false if with no else returns 0), so the message
  # below reported "exited 0" for a stage that actually died. `&& return 0` keeps
  # the failing code in $? for the rc capture.
  Rscript "$script" "$@" && return 0
  local rc=$?
  echo "" >&2
  echo "ERROR: stage ${label} exited ${rc}. Its R error is in the stage log, not on this console." >&2
  local newest="" f
  for f in "${WORK_DIR}"/*fit_models.log "${WORK_DIR}"/*discrimination_metrics.log "${WORK_DIR}"/*log.txt; do
    [[ -f "$f" ]] || continue
    [[ -z "$newest" || "$f" -nt "$newest" ]] && newest="$f"
  done
  if [[ -n "$newest" ]]; then
    echo "--- last 40 lines of $(basename "$newest") ---" >&2
    tail -40 "$newest" >&2
    echo "--- full log: $newest ---" >&2
  else
    echo "  (no stage log found under $WORK_DIR — the stage may have failed before opening one)" >&2
  fi
  die "stage ${label} failed."
}

# ============================== Validate required args ==============================
[[ -n "$TRAIT" ]] || die "--trait is required. E.g.: --trait=T1D_to_CAD"

# ============================== Parse trait name ==============================
# Formats supported:
#   ONSET_to_OUTCOME  (underscore delimited, e.g. T1D_to_CAD)
#   ONSETtoOUTCOME    (camel-style, e.g. PARKINSONtoDEMENTIA, CADtoHEARTFAIL)
# Check _to_ first so T1D_to_CAD does not split on the lowercase 'to' inside '_to_'

if [[ -n "$ONSET_NAME_OVERRIDE" && -n "$OUTCOME_NAME_OVERRIDE" ]]; then
  ONSET_NAME="$ONSET_NAME_OVERRIDE"
  OUTCOME_NAME="$OUTCOME_NAME_OVERRIDE"
elif [[ "$TRAIT" == *_to_* ]]; then
  ONSET_NAME="${TRAIT%%_to_*}"
  OUTCOME_NAME="${TRAIT##*_to_}"
elif [[ "$TRAIT" == *to* ]]; then
  # lowercase 'to' is the separator; uppercase disease names never contain 'to'
  ONSET_NAME="${TRAIT%%to*}"
  OUTCOME_NAME="${TRAIT##*to}"
else
  die "Cannot split trait '$TRAIT' into onset/outcome parts.
  Expected format: ONSET_to_OUTCOME  (e.g. T1D_to_CAD)
               or: ONSETtoOUTCOME   (e.g. PARKINSONtoDEMENTIA)
  Use --onset-name and --outcome-name to override manually."
fi

[[ -n "$ONSET_NAME"   ]] || die "Parsed onset name is empty for trait '$TRAIT'. Use --onset-name to set manually."
[[ -n "$OUTCOME_NAME" ]] || die "Parsed outcome name is empty for trait '$TRAIT'. Use --outcome-name to set manually."

# ============================== Resolve file paths ==============================
# Inputs are supplied explicitly (required for a fitting run; validated below when RUN_01=1).
ONSET_PRS_PATH="$ONSET_PRS_FILE_OVERRIDE"
PROGRESSION_PRS_PATH="$PROGRESSION_PRS_FILE_OVERRIDE"
OUTCOME_PRS_PATH="$OUTCOME_PRS_FILE_OVERRIDE"
PHENO_PATH="$PHENO_FILE_OVERRIDE"
# Covariates default to the phenotype file, but may live in a separate file.
if [[ -n "$COV_FILE" ]]; then
  [[ -f "$COV_FILE" ]] || die "--cov-file not found: $COV_FILE"
  COV_PATH="$COV_FILE"
else
  COV_PATH="$PHENO_PATH"   # covariate columns live in the same phenotype file
fi

# ============================== Validate files ==============================
# Only stage 01 consumes the PRS/phenotype inputs. A downstream-only rerun
# (--run-01=0) works purely from the RDS outputs, so skip these checks then and
# rely on the per-stage "Missing 01 outputs" guards below.
if [[ "$RUN_01" == "1" ]]; then
[[ -f "$ONSET_PRS_PATH" ]] || die "Onset PRS file not provided or not found: '$ONSET_PRS_PATH'
  Pass it explicitly with --onset-prs-file=FILE"

# All three PRS are required (onset, progression, outcome) — no partial runs.
[[ -f "$PROGRESSION_PRS_PATH" ]] || die "Progression PRS file not provided or not found: '$PROGRESSION_PRS_PATH'
  Pass it explicitly with --progression-prs-file=FILE
  (the three-PRS set is mandatory)"

[[ -f "$OUTCOME_PRS_PATH" ]] || die "Outcome PRS file not provided or not found: '$OUTCOME_PRS_PATH'
  Pass it explicitly with --outcome-prs-file=FILE"

[[ -f "$PHENO_PATH" ]] || die "Phenotype file not provided or not found: '$PHENO_PATH'
  Pass it explicitly with --pheno-file=FILE"
fi

# R scripts
S01="${SCRIPT_DIR}/01_fit_models.R"
S02="${SCRIPT_DIR}/02_kaplan_meier.R"
S03="${SCRIPT_DIR}/03_discrimination_metrics.R"
S04="${SCRIPT_DIR}/04_reclassification.R"
S05="${SCRIPT_DIR}/05_calibration.R"
S06="${SCRIPT_DIR}/06_cumulative_incidence.R"
SUTIL="${SCRIPT_DIR}/prs_risk_utils.R"
for s in "$S01" "$S02" "$S03" "$S04" "$S05" "$S06" "$SUTIL"; do
  [[ -f "$s" ]] || die "R script not found: $s  (set --script-dir)"
done

# A supplied --recruit-file path must exist, in EITHER layer. Otherwise a typo is
# silently ignored (Layer 1 forwards it only when the file exists), quietly dropping
# the recruitment-timing QC. Fail loudly instead. (Empty = not requested, is fine.)
if [[ -n "$RECRUIT_FILE" && ! -f "$RECRUIT_FILE" ]]; then
  die "--recruit-file was given but does not exist: $RECRUIT_FILE"
fi

# A supplied --pop-file keep-list must exist (else the run would silently use the whole input).
if [[ -n "$POP_FILE" && ! -f "$POP_FILE" ]]; then
  die "--pop-file was given but does not exist: $POP_FILE"
fi

# Prospective mode needs a recruitment-age file — but only for stage 01, which
# derives the prospective cohort. A downstream-only rerun (--run-01=0) works from
# the RDS, which already holds the derived prospective variables.
if [[ "$RUN_01" == "1" && "$ANALYSIS_MODE" == "prospective" ]]; then
  [[ -n "$RECRUIT_FILE" && -f "$RECRUIT_FILE" ]] || die "--analysis-mode=prospective requires --recruit-file=FILE (IID + $RECRUIT_AGE_COL)"
  # Prospective re-indexing and a disease-duration lag define different time origins.
  if [[ "$LAG_DAYS" != "0" && -n "$LAG_DAYS" ]]; then
    die "--analysis-mode=prospective cannot be combined with a lag sweep. Set --lag-days=0 for prospective, and run the lag sensitivity as a separate gwas_aligned job."
  fi
fi
# Validate here as well as in stage 01, so a typo fails before the run starts
# rather than after the PRS files have been read.
case "$PROSPECTIVE_T1" in
  both|incident) ;;
  *) die "--prospective-t1 must be 'both' or 'incident' (got '$PROSPECTIVE_T1')." ;;
esac
if [[ "$ANALYSIS_MODE" != "prospective" && "$PROSPECTIVE_T1" != "both" ]]; then
  die "--prospective-t1=$PROSPECTIVE_T1 only applies to --analysis-mode=prospective."
fi

# --out-dir is required (give every run its own directory). Two subfolders are created:
#   intermediate/  individual-level + working files (RDS design matrices, per-person
#                  risk movement, logs) — do NOT share these outside your biobank.
#   final/         aggregate, shareable summary tables and figures.
[[ -n "$OUT_DIR" ]] || die "--out-dir is required (give each trajectory/ancestry/method/layer run its own directory)."
WORK_DIR="${OUT_DIR}/intermediate"
FINAL_DIR="${OUT_DIR}/final"
mkdir -p "$WORK_DIR" "$FINAL_DIR"
echo "Output directory:        $OUT_DIR"
echo "  intermediate (private): $WORK_DIR"
echo "  final (candidate aggregate outputs; local review required): $FINAL_DIR"

# Single base stem, with PRS_METHOD embedded so a P+T run cannot overwrite a PRS-CS
# run for the same cohort/trajectory. PREFIX, LAG_PREFIX and the lag-summary stem
# all derive from BASE_STEM so they cannot drift.
# Provenance defaults + completeness (post-parse). Target ancestry defaults to the
# analysis ancestry; warn (do not fail at rc2) on blank provenance so a site knows
# these fields are MISSING in the returned manifest.
_prov_missing=()
[[ -z "$DISCOVERY_ANCESTRY" ]]   && _prov_missing+=("--discovery-ancestry")
if (( ${#_prov_missing[@]} )); then
  echo "⚠ Provenance fields not supplied (recorded as MISSING in the manifest): ${_prov_missing[*]}" >&2
  echo "  Supply them for a production run so the manifest records the discovery ancestry." >&2
fi

# Identity safety: the MGBB/EUR defaults are convenient locally but mislabel outputs at
# another site. Warn loudly when an identity flag was NOT passed and its default is used
# (outputs are still computed correctly).
_id_defaults=()
[[ "${COHORT_EXPLICIT:-0}" != "1" ]]   && _id_defaults+=("--cohort (using default '$COHORT')")
[[ "${ANCESTRY_EXPLICIT:-0}" != "1" ]] && _id_defaults+=("--ancestry (using default '$ANCESTRY')")
if (( ${#_id_defaults[@]} )); then
  echo "⚠ SITE IDENTITY: using built-in defaults for: ${_id_defaults[*]}" >&2
  echo "  Outputs will be LABELLED '${COHORT}_${ANCESTRY}'. If that is not your site, pass" >&2
  echo "  --cohort / --ancestry explicitly so results are not mislabelled." >&2
fi

BASE_STEM="${COHORT}_${ANCESTRY}_${PRS_METHOD}_${TRAIT}"
PREFIX="$BASE_STEM"
# 01 appends the analysis mode to its auto-prefix for non-default modes; mirror
# that here so FM_MAIN/DT_MAIN (and the downstream stages) resolve in prospective mode.
[[ "$ANALYSIS_MODE" != "gwas_aligned" ]] && PREFIX="${PREFIX}_${ANALYSIS_MODE}"
# 01 also appends an "incidentLag<N>d" component when --incident_lag_days > 0 (passed
# only in prospective mode). Mirror it too, or the driver's FM_MAIN/DT_MAIN would point
# at "..._prospective_01_..." while 01 self-names "..._prospective_incidentLag<N>d_01_...",
# and every downstream stage would die "Missing 01 outputs" though 01 succeeded. awk gives
# a float-safe >0 test. (The main stage-01 call always runs lag_days=0, so no "lag<N>d".)
if [[ "$ANALYSIS_MODE" == "prospective" ]] && awk "BEGIN{exit !((${INCIDENT_LAG_DAYS:-0})>0)}"; then
  PREFIX="${PREFIX}_incidentLag${INCIDENT_LAG_DAYS}d"
fi

# ============================== Package preflight ==============================
# CHECK versions; do NOT install. Silent runtime installation fails in biobanks
# without internet access and lets sites drift to different package versions —
# unacceptable for a multi-site consortium. Set PRS_ALLOW_INSTALL=1 to opt in.
# Cairo is optional (figures only) and needs X11/XQuartz; absence is not fatal.
info "Checking R packages and versions"
# The authoritative fit manifest (input checksums, fit config) is written ONLY when
# stage 01 fits the models (RUN_01=1). A downstream-only rerun (RUN_01=0) must NOT
# overwrite it with a manifest that lacks the input checksums — it writes a separate
# downstream_rerun_manifest.csv instead (current package versions + script checksums),
# leaving the original session_manifest.csv (preserved in intermediate/) untouched.
if [[ "$RUN_01" == "1" ]]; then
  MANIFEST_FILE="${WORK_DIR}/${PREFIX}_session_manifest.csv"
else
  MANIFEST_FILE="${WORK_DIR}/${PREFIX}_downstream_rerun_manifest.csv"
  info "Downstream-only rerun (--run-01=0): preserving the original session_manifest.csv / run_configuration.csv; writing $(basename "$MANIFEST_FILE")"
fi
# Stage-aware required set: only demand packages the ENABLED stages/options need,
# so a site running just the core (01+03) is not blocked by, say, missing pec.
REQ_PKGS="optparse data.table survival splines riskRegression ggplot2 prodlim"
[[ "$RUN_02" == "1" ]] && REQ_PKGS+=" survminer scales RColorBrewer"
[[ "$RUN_04" == "1" ]] && REQ_PKGS+=" gridExtra ggpubr"
[[ "$RUN_05" == "1" ]] && REQ_PKGS+=" rms pec"
[[ "$RUN_06" == "1" ]] && REQ_PKGS+=" boot rlang scales"
[[ "${CALC_UNO:-0}" == "1" ]] && REQ_PKGS+=" survAUC"
[[ "${AUC_ENGINE:-riskRegression}" == "timeROC" || "${AUC_ENGINE:-riskRegression}" == "both" ]] && REQ_PKGS+=" timeROC"
# Provenance checksums, computed in bash (the manifest is written in R). sha256
# via shasum, falling back to sha256sum; absent/unreadable files record a marker.
_sha() {
  local f="$1"
  [[ -f "$f" ]] || { echo "absent"; return; }
  if command -v shasum >/dev/null 2>&1;      then shasum -a 256 "$f" | awk '{print $1}'
  elif command -v sha256sum >/dev/null 2>&1;  then sha256sum "$f"    | awk '{print $1}'
  else echo "unavailable"; fi
}
CKSUM_LINES=""
# The driver decides arguments, publishing rules and provenance, so checksum it too.
CKSUM_LINES+="driver:00_run_pipeline_gbmi.sh=$(_sha "${BASH_SOURCE[0]}")"$'\n'
for s in "$S01" "$S02" "$S03" "$S04" "$S05" "$S06" "$SUTIL"; do
  CKSUM_LINES+="script:$(basename "$s")=$(_sha "$s")"$'\n'
done
CKSUM_LINES+="test_script=$(_sha "${SCRIPT_DIR}/test_prs_risk_utils.R")"$'\n'
[[ -f "${SCRIPT_DIR}/smoke_test.sh" ]]      && CKSUM_LINES+="script:smoke_test.sh=$(_sha "${SCRIPT_DIR}/smoke_test.sh")"$'\n'
[[ -f "${SCRIPT_DIR}/RUN_PRS_METRICS_GBMI.md" ]] && CKSUM_LINES+="doc:RUN_PRS_METRICS_GBMI.md=$(_sha "${SCRIPT_DIR}/RUN_PRS_METRICS_GBMI.md")"$'\n'
if [[ "$RUN_01" == "1" ]]; then
  CKSUM_LINES+="prs_onset=$(_sha "$ONSET_PRS_PATH")"$'\n'
  CKSUM_LINES+="prs_progression=$(_sha "$PROGRESSION_PRS_PATH")"$'\n'
  CKSUM_LINES+="prs_outcome=$(_sha "$OUTCOME_PRS_PATH")"$'\n'
  CKSUM_LINES+="phenotype=$(_sha "$PHENO_PATH")"$'\n'
  [[ -n "$COV_FILE" ]] && CKSUM_LINES+="covariate_file=$(_sha "$COV_FILE")"$'\n'
  [[ -n "$RECRUIT_FILE" ]] && CKSUM_LINES+="recruit_file=$(_sha "$RECRUIT_FILE")"$'\n'
  [[ -n "$POP_FILE" ]] && CKSUM_LINES+="pop_file=$(_sha "$POP_FILE")"$'\n'
fi
PRS_MANIFEST="$MANIFEST_FILE" PRS_VERSION="$PIPELINE_VERSION" PRS_REQ_PKGS="$REQ_PKGS" \
  PRS_CKSUMS="$CKSUM_LINES" PRS_METHOD="$PRS_METHOD" \
  PRS_DISC_ANC="$DISCOVERY_ANCESTRY" PRS_TARG_ANC="$ANCESTRY" \
  Rscript - <<'RCHECK'
# Enforce R version first (was recorded but not checked).
if (getRversion() < "4.1.0") {
  cat("ERROR: R", as.character(getRversion()), "is too old; R >= 4.1.0 is required.\n")
  quit(status=1)
}
min <- c(optparse="1.7.0", data.table="1.14.0", survival="3.5", splines="0",
         riskRegression="2023.03.22", ggplot2="3.4.0", survminer="0.4.9",
         scales="1.2.0", RColorBrewer="1.1", rms="6.0", pec="2022.05.04",
         boot="1.3", rlang="1.0.0", gridExtra="2.3", ggpubr="0.6.0",
         timeROC="0.4", survAUC="1.1", prodlim="2019.11.13")
req <- strsplit(trimws(Sys.getenv("PRS_REQ_PKGS")), "\\s+")[[1]]
pkgs <- intersect(names(min), unique(req))
pkgs <- setdiff(pkgs, "splines")   # base R, no version to check
inst <- rownames(installed.packages())
missing <- setdiff(pkgs, inst)
allow  <- identical(Sys.getenv("PRS_ALLOW_INSTALL"), "1")
if (length(missing) && allow) {
  cat("Installing missing packages (PRS_ALLOW_INSTALL=1):", paste(missing, collapse=", "), "\n")
  install.packages(missing, repos="https://cloud.r-project.org/", quiet=TRUE)
  inst <- rownames(installed.packages())
  missing <- setdiff(pkgs, inst)
}
rows <- list(); bad <- character(0)
for (p in pkgs) {
  if (!p %in% inst) { cat(sprintf("  %-16s MISSING\n", p)); bad <- c(bad, p); next }
  have <- as.character(utils::packageVersion(p))
  old  <- utils::compareVersion(have, min[[p]]) < 0
  cat(sprintf("  %-16s %-12s %s\n", p, have, if (old) paste0("TOO OLD (need >= ", min[[p]], ")") else "OK"))
  if (old) bad <- c(bad, p)
  rows[[length(rows)+1]] <- data.frame(component=p, version=have, required=min[[p]],
                                       status=if (old) "TOO_OLD" else "OK")
}
cairo_ok <- requireNamespace("Cairo", quietly=TRUE)
cat(sprintf("  %-16s %s\n", "Cairo (optional)",
            if (cairo_ok) "OK" else "absent -> figures via grDevices fallback"))
rows[[length(rows)+1]] <- data.frame(component="R", version=as.character(getRversion()),
                                     required="4.1.0",
                                     status=if (getRversion() < "4.1.0") "TOO_OLD" else "OK")
rows[[length(rows)+1]] <- data.frame(component="pipeline_version",
                                     version=Sys.getenv("PRS_VERSION"), required="", status="OK")
# Provenance: PRS method, optional pass-through fields, and sha256 checksums.
prov <- c(prs_method=Sys.getenv("PRS_METHOD"),
          discovery_ancestry=Sys.getenv("PRS_DISC_ANC"),
          target_ancestry=Sys.getenv("PRS_TARG_ANC"))
for (nm in names(prov)) {
  val <- prov[[nm]]
  # Blank provenance is recorded as MISSING (not OK): prs_method always has a
  # default, so it stays OK; the optional pass-through fields flag when unset.
  st <- if (nm != "prs_method" && (is.null(val) || !nzchar(val))) "MISSING" else "OK"
  rows[[length(rows)+1]] <- data.frame(component=nm, version=val, required="", status=st)
}
ck <- Sys.getenv("PRS_CKSUMS")
if (nzchar(ck)) for (ln in strsplit(ck, "\n")[[1]]) {
  ln <- trimws(ln); if (!nzchar(ln)) next
  kv <- strsplit(ln, "=", fixed=TRUE)[[1]]
  rows[[length(rows)+1]] <- data.frame(component=paste0("sha256:", kv[1]),
                                       version=if (length(kv)>1) kv[2] else "",
                                       required="", status="OK")
}
mf <- Sys.getenv("PRS_MANIFEST")
if (nzchar(mf)) {
  dir.create(dirname(mf), recursive=TRUE, showWarnings=FALSE)
  write.csv(do.call(rbind, rows), mf, row.names=FALSE)
}
if (length(bad)) {
  cat("\nERROR: missing or out-of-date packages:", paste(unique(bad), collapse=", "), "\n")
  cat("Install/update them, then re-run:\n")
  cat('  install.packages(c(', paste(sprintf('"%s"', unique(bad)), collapse=", "), '))\n', sep="")
  cat("(or set PRS_ALLOW_INSTALL=1 to let the pipeline install them for you)\n")
  quit(status=1)
}
RCHECK
[[ $? -eq 0 ]] || die "R package preflight failed — see the versions listed above."

# Export the resolved script directory so the self-locating R scripts (test, 04-06)
# find prs_risk_utils.R WITHOUT parsing their own --file path. Some macOS/OneDrive
# setups encode the spaces in --file as ~+~, which breaks normalizePath(); this env
# var carries the real path (spaces intact) to every child Rscript. SCRIPT_DIR is
# absolute (realpath default, or the user's --script-dir).
export PRS_SCRIPT_DIR="$SCRIPT_DIR"

# Run the helper unit tests (default on): catches a broken environment before a
# long run, and makes the README's claim true. Skip with --run-tests=0.
# A MISSING test file must fail loudly, not be silently skipped — the absolute-risk
# helper it validates drives stages 04-06.
UNIT_TESTS_REQUESTED="$RUN_TESTS"; UNIT_TESTS_PASSED="not_run"
if [[ "$RUN_TESTS" == "1" ]]; then
  TEST_FILE="${SCRIPT_DIR}/test_prs_risk_utils.R"
  [[ -f "$TEST_FILE" ]] || die "Unit-test file not found: $TEST_FILE
  It is part of the distributed package. Use --run-tests=0 only for a documented rerun."
  info "Running helper unit tests (test_prs_risk_utils.R)"
  if Rscript "$TEST_FILE"; then UNIT_TESTS_PASSED="TRUE"
  else UNIT_TESTS_PASSED="FALSE"; die "Helper unit tests failed — aborting. (Set --run-tests=0 to skip.)"; fi
fi
# Record the test outcome in the manifest (written by the preflight above, before
# the tests could run).
if [[ -f "$MANIFEST_FILE" ]]; then
  printf 'unit_tests_requested,%s,,OK\nunit_tests_passed,%s,,%s\n' \
    "$UNIT_TESTS_REQUESTED" "$UNIT_TESTS_PASSED" \
    "$([[ "$UNIT_TESTS_PASSED" == "TRUE" || "$UNIT_TESTS_PASSED" == "not_run" ]] && echo OK || echo FAILED)" \
    >> "$MANIFEST_FILE"
fi

# The time scale must be stated explicitly. Falling back to the built-in
# secondTime default silently changes the analysis definition between sites,
# which breaks cross-biobank harmonisation.
if [[ "$RUN_01" == "1" && -z "$AGE_EXIT_COL" && "$TIME_EXIT_EXPLICIT" != "1" ]]; then
  die "Specify the time scale explicitly:\n  --age-exit-col=NAME   (recommended: STOP = age_exit - age_t1 derived in-script)\n  or --time-exit=NAME   (a precomputed time-since-T1 column, e.g. secondTime)\nRelying on the built-in default is not permitted for consortium runs."
fi

FM_MAIN="${WORK_DIR}/${PREFIX}_01_fit_models_fitted_models.rds"
DT_MAIN="${WORK_DIR}/${PREFIX}_01_fit_models_data_processed.rds"

# Stages 03/04 gate time-point inclusion (hence the smallest counts they publish)
# on --min_events, which the driver previously left at each stage's own default of
# 10 regardless of --min-cell-count. Make that floor the STRICTER of the reliability
# minimum (10) and the disclosure threshold, so raising --min-cell-count tightens
# 03/04 too without ever letting the reliability floor drop below 10.
MIN_EVENTS_0304=$(( MIN_CELL_COUNT > 10 ? MIN_CELL_COUNT : 10 ))

# ============================== Detect score column ==============================
# Peek at the header of the onset sscore file; fall back to PRS_SCORE_COL default.
# 01_fit_models.R also has grep fallback for the extra PRS files.
detect_score_col() {
  local file="$1" preferred="$2"
  # Read first line (header), check if preferred column exists
  local header
  header=$(head -1 "$file")
  if echo "$header" | tr '\t' '\n' | grep -qxF "$preferred"; then
    echo "$preferred"
  else
    # Try common PRScs column names in order
    for cand in SCORE1_SUM SCORESUM SCORE_SUM; do
      if echo "$header" | tr '\t' '\n' | grep -qxF "$cand"; then
        echo "$cand"
        return
      fi
    done
    # Fall back to user-supplied default
    echo "$preferred"
  fi
}

# Only stage 01 reads the PRS files. A downstream-only rerun (--run-01=0) works
# purely from the RDS outputs, so do NOT touch the PRS files then — detect_score_col
# does `head -1 "$file"` and would fail if the score files are absent/moved.
if [[ "$RUN_01" == "1" ]]; then
  ONSET_COL=$(detect_score_col "$ONSET_PRS_PATH" "$PRS_SCORE_COL")
  OUTCOME_COL=$(detect_score_col "$OUTCOME_PRS_PATH" "$PRS_SCORE_COL")
  PROGRESSION_COL=$(detect_score_col "$PROGRESSION_PRS_PATH" "$PRS_SCORE_COL")
else
  ONSET_COL="(unused: --run-01=0)"; OUTCOME_COL="$ONSET_COL"; PROGRESSION_COL="$ONSET_COL"
fi

# ============================== Print configuration ==============================
info "GBMI Multi-PRS Pipeline — Configuration"
printf "  %-26s %s\n"  "Trait:"             "$TRAIT"
printf "  %-26s %s\n"  "Onset name:"        "$ONSET_NAME"
printf "  %-26s %s\n"  "Progression name:"  "$TRAIT  (full compound trait)"
printf "  %-26s %s\n"  "Outcome name:"      "$OUTCOME_NAME"
printf "  %-26s %s\n"  "Cohort/Ancestry:"   "${COHORT}/${ANCESTRY}"
printf "  %-26s %s\n"  "Analysis mode:"     "$ANALYSIS_MODE"
echo   ""
printf "  %-26s %s\n"  "Onset sscore:"      "$ONSET_PRS_PATH  [$ONSET_COL]"
printf "  %-26s %s\n"  "Progression sscore:" "$PROGRESSION_PRS_PATH  [$PROGRESSION_COL]"
printf "  %-26s %s\n"  "Outcome sscore:"    "$OUTCOME_PRS_PATH  [$OUTCOME_COL]"
printf "  %-26s %s\n"  "Phenotype file:"    "$PHENO_PATH"
[[ "$ANALYSIS_MODE" == "prospective" ]] && printf "  %-26s %s\n" "Recruitment file:" "$RECRUIT_FILE"
printf "  %-26s %s\n"  "Output dir:"        "$OUT_DIR"
echo   ""
if [[ -n "$AGE_EXIT_COL" ]]; then
  printf "  %-26s %s\n" "Time scale:" "STOP = ${AGE_EXIT_COL} - ${AGE_T1_COL} (derived)"
else
  printf "  %-26s %s\n" "Time scale:" "STOP = ${TIME_EXIT} (column)"
fi
printf "  %-26s %s\n"  "Age covariate:"     "ns(${AGE_T1_COL}, df=${AGE_SPLINE_DF})"
printf "  %-26s %s\n"  "Covariates:"        "$COVARIATES"
printf "  %-26s %s\n"  "Lag sweep (days):"  "$LAG_DAYS"
printf "  %-26s %s\n"  "Min events (total):" "$MIN_EVENTS_TOTAL"
if [[ -n "$MODELS_SPEC" ]]; then
  printf "  %-26s %s\n" "Models spec:"      "$MODELS_SPEC"
else
  printf "  %-26s %s\n" "Models (M0-M7):"   "Base, Onset, Outcome, Progression, Onset_Progression, Outcome_Progression, Onset_Outcome, Onset_Progression_Outcome"
fi
echo ""

# ---- Shared argument groups reused by the main run and the lag sweep ----------
MULTI_PRS_ARGS=(
  --prs_labels   "Onset,Progression,Outcome"
  --prs_types    "onset,progression,outcome"
  --prs_col_list "${ONSET_COL},${PROGRESSION_COL},${OUTCOME_COL}"
  --prs_files    "${PROGRESSION_PRS_PATH},${OUTCOME_PRS_PATH}"
)
[[ -n "$MODELS_SPEC" ]] && MULTI_PRS_ARGS+=(--models "$MODELS_SPEC")

# Time-scale args. Always pass the age-at-T1 (spline) column and a --time_exit
# fallback; if an age-at-exit column is given too, 01 prefers STOP = age_exit - age_t1.
AGE_TIME_ARGS=(--age_t1_col "$AGE_T1_COL" --age_spline_df "$AGE_SPLINE_DF"
               $(opt --time_exit "$TIME_EXIT") $(opt --time_entry "$TIME_ENTRY"))
[[ -n "$AGE_EXIT_COL" ]] && AGE_TIME_ARGS+=(--age_exit_col "$AGE_EXIT_COL")

# Stage 04 previously ran with its own defaults (reference=base, comparison=minimal),
# i.e. approximately M1 vs M0 — not the manuscript's M5 vs M2. Derive the pair here.
# RECLASS_PAIR is "COMPARISON:REFERENCE" (the richer model first), e.g. M5:M2.
RECLASS_COMP="${RECLASS_PAIR%%:*}"
RECLASS_REF="${RECLASS_PAIR##*:}"
[[ -n "$RECLASS_COMP" && -n "$RECLASS_REF" && "$RECLASS_COMP" != "$RECLASS_REF" ]] \
  || die "--reclass-pair must be COMPARISON:REFERENCE with two different models (got '$RECLASS_PAIR')"
# Key paired contrasts reported by stage 03.
CONTRAST_PAIRS="${CONTRAST_PAIRS:-M4:M1,M5:M2,M7:M6}"

MODE_ARGS=(--analysis_mode "$ANALYSIS_MODE" --min_events_total "$MIN_EVENTS_TOTAL"
           --min_cell_count "$MIN_CELL_COUNT")
[[ "$CALC_UNO" == "1" ]] && MODE_ARGS+=(--calculate_uno)
[[ -n "$CLINICAL_COVARIATES" ]] && MODE_ARGS+=(--clinical_covariates "$CLINICAL_COVARIATES")
# Optional keep-list: stage 01 subsets the sample before every downstream stage, and MODE_ARGS is
# passed to both the main fit and the lag re-fits, so the whole run honours the same keep-list.
[[ -n "$POP_FILE" ]] && MODE_ARGS+=(--pop_file "$POP_FILE")
# Forward the recruitment-age file whenever it is supplied, not only in prospective
# mode: prospective needs it to re-index, but gwas_aligned uses it for the
# descriptive recruitment-timing QC (prevalent / incident / historical counts).
if [[ -n "$RECRUIT_FILE" && -f "$RECRUIT_FILE" ]]; then
  MODE_ARGS+=(--recruit_file "$RECRUIT_FILE" --recruit_age_col "$RECRUIT_AGE_COL")
fi
if [[ "$ANALYSIS_MODE" == "prospective" ]]; then
  MODE_ARGS+=(--incident_lag_days "$INCIDENT_LAG_DAYS" --prospective_t1 "$PROSPECTIVE_T1")
fi

# Shareable run configuration (published) — records exactly how this job was run, so a
# central analyst does not have to reconstruct settings from the command line. No counts.
# Written ONLY when stage 01 fits the models; a downstream-only rerun preserves the
# original (its settings describe the fit, which the rerun does not change).
RUNCFG="${WORK_DIR}/${PREFIX}_run_configuration.csv"
if [[ "$RUN_01" == "1" ]]; then
{
  echo "key,value"
  echo "pipeline_version,${PIPELINE_VERSION}"
  echo "cohort,${COHORT}"
  echo "ancestry,${ANCESTRY}"
  echo "prs_method,${PRS_METHOD}"
  echo "trait,${TRAIT}"
  echo "analysis_mode,${ANALYSIS_MODE}"
  echo "discovery_ancestry,${DISCOVERY_ANCESTRY}"
  echo "target_ancestry,${ANCESTRY}"
  echo "age_t1_col,${AGE_T1_COL}"
  echo "age_exit_col,${AGE_EXIT_COL}"
  echo "time_exit_col,${TIME_EXIT}"
  echo "status_col,${STATUS_COL}"
  echo "covariates,\"${COVARIATES}\""
  echo "clinical_covariates,\"${CLINICAL_COVARIATES}\""
  echo "time_points,\"${TIME_POINTS}\""
  echo "lag_days,\"${LAG_DAYS}\""
  echo "min_events_total,${MIN_EVENTS_TOTAL}"
  echo "min_cell_count,${MIN_CELL_COUNT}"
  echo "contrast_pairs,\"${CONTRAST_PAIRS}\""
  echo "score_roles,\"${SCORE_ROLE}\""
  echo "calculate_rmst,${CALC_RMST}"
  echo "calculate_uno,${CALC_UNO}"
  echo "stages,\"01=${RUN_01} 02=${RUN_02} 03=${RUN_03} 04=${RUN_04} 05=${RUN_05} 06=${RUN_06}\""
  echo "recruit_file_provided,$([[ -n "$RECRUIT_FILE" ]] && echo yes || echo no)"
  echo "pop_file_provided,$([[ -n "$POP_FILE" ]] && echo yes || echo no)"
  # Additional resolved options, so the config captures every setting that shapes outputs.
  echo "age_spline_df,${AGE_SPLINE_DF}"
  echo "sex_col,${SEX_COL}"
  echo "recruit_age_col,${RECRUIT_AGE_COL}"
  echo "incident_lag_days,${INCIDENT_LAG_DAYS}"
  echo "prospective_t1,${PROSPECTIVE_T1}"
  echo "prs_quantiles,\"${PRS_QUANTILES}\""
  echo "prs_extremes,\"${PRS_EXTREMES}\""
  echo "role_quantiles,\"${ROLE_QUANTILES}\""
  echo "use_supported_times,${USE_SUPPORTED_TIMES}"
  echo "risk_model,${RISK_MODEL}"
  echo "reclass_pair,${RECLASS_PAIR}"
  echo "risk_thresholds,\"${RISK_THRESHOLDS}\""
  echo "rmst_n_boot,${RMST_N_BOOT}"
  echo "bootstrap_risk_differences,${BOOTSTRAP_RISK_DIFFERENCES}"
  echo "enable_clinical_utility,${ENABLE_CLINICAL_UTILITY}"
} > "$RUNCFG"
info "Wrote run configuration: $(basename "$RUNCFG")"
fi   # end RUN_01==1 run_configuration guard

# ============================== Stage 01 (main, lag 0) ==============================
if [[ "$RUN_01" == "1" ]]; then
  info "01 — fitting multi-PRS Cox models (M0-M7)"
  run_stage 01 "$S01" \
    --prs_file        "$ONSET_PRS_PATH" \
    --prs_col         "$ONSET_COL" \
    --pheno_file      "$PHENO_PATH" \
    --cov_file        "$COV_PATH" \
    --status_col      "$STATUS_COL" \
    "${AGE_TIME_ARGS[@]}" \
    "${MODE_ARGS[@]}" \
    $(opt --time_summary_cols "$TIME_SUMMARY_COLS") \
    --cohort          "$COHORT" \
    --ancestry        "$ANCESTRY" \
    --prs_method      "$PRS_METHOD" \
    --pheno_name      "$TRAIT" \
    --covariates      "$COVARIATES" \
    --sex_col         "$SEX_COL" \
    --prs_quantiles   "$PRS_QUANTILES" \
    --prs_extremes    "$PRS_EXTREMES" \
    --time_points     "$TIME_POINTS" \
    --outdir          "$WORK_DIR" \
    "${MULTI_PRS_ARGS[@]}"
fi

# ============================== Stage 03 ==============================
if [[ "$RUN_03" == "1" ]]; then
  info "03 — discrimination metrics (AUC, Brier, iAUC, IBS)"
  [[ -f "$FM_MAIN" && -f "$DT_MAIN" ]] || die "Missing 01 outputs — run stage 01 first."
  run_stage 03 "$S03" \
    --fitted_models  "$FM_MAIN" \
    --data_processed "$DT_MAIN" \
    --prefix "$PREFIX" \
    --contrast_pairs "$CONTRAST_PAIRS" \
    --time_points "$TIME_POINTS" \
    --min_events "$MIN_EVENTS_0304" \
    $(uset) \
    --outdir "$WORK_DIR"
fi

# ============================== Stage 02 ==============================
if [[ "$RUN_02" == "1" ]]; then
  info "02 — Kaplan-Meier curves by PRS quantile groups"
  [[ -f "$FM_MAIN" && -f "$DT_MAIN" ]] || die "Missing 01 outputs — run stage 01 first."
  KM_ARGS=(--models_file "$FM_MAIN" --data_file "$DT_MAIN" --outdir "$WORK_DIR"
           --prefix "$PREFIX" --score_role "$SCORE_ROLE" --role_quantiles "$ROLE_QUANTILES"
           --min_cell_count "$MIN_CELL_COUNT")
  # RMST restricted to the shared evaluation horizons, so restricted means are
  # comparable across PRS strata (the script's default tau is each stratum's own
  # max follow-up, which is not). --rmst_tau is the CANDIDATE set; stage 02
  # intersects it with the stage-01 supported horizons and the per-group support
  # rule before integrating.
  if [[ "$CALC_RMST" == "1" ]]; then
    KM_ARGS+=(--calculate_rmst --rmst_tau "$TIME_POINTS")
    # Without bootstrap the SE/CI/p columns are all NA, so rmst_diff is
    # uninterpretable. ~25 s at 200 reps on ~6.6k rows.
    [[ "$RMST_N_BOOT" -gt 0 ]] && KM_ARGS+=(--rmst_bootstrap --rmst_n_boot "$RMST_N_BOOT")
  fi
  run_stage 02 "$S02" "${KM_ARGS[@]}"
fi

# ============================== Stage 04 ==============================
if [[ "$RUN_04" == "1" ]]; then
  info "04 — reclassification (IDI / NRI)"
  [[ -f "$FM_MAIN" && -f "$DT_MAIN" ]] || die "Missing 01 outputs — run stage 01 first."
  run_stage 04 "$S04" \
    --models_file    "$FM_MAIN" \
    --data_file      "$DT_MAIN" \
    --risk_thresholds "$RISK_THRESHOLDS" \
    --prefix "$PREFIX" \
    --reference_model "$RECLASS_REF" \
    --comparison_model "$RECLASS_COMP" \
    --time_points "$TIME_POINTS" \
    --min_events "$MIN_EVENTS_0304" \
    $(uset) \
    --outdir "$WORK_DIR"
fi

# ============================== Stage 05 ==============================
if [[ "$RUN_05" == "1" ]]; then
  info "05 — calibration"
  [[ -f "$FM_MAIN" && -f "$DT_MAIN" ]] || die "Missing 01 outputs — run stage 01 first."
  run_stage 05 "$S05" \
    --models_file "$FM_MAIN" \
    --data_file   "$DT_MAIN" \
    --n_groups 10 \
    --min_cell_count "$MIN_CELL_COUNT" \
    --prefix "$PREFIX" \
    $(uset) \
    --outdir "$WORK_DIR"
fi

# ============================== Stage 06 ==============================
if [[ "$RUN_06" == "1" ]]; then
  info "06 — cumulative incidence & clinical utility"
  [[ -f "$FM_MAIN" ]] || die "Missing 01 output — run stage 01 first."
  # 06 expects an actual stratum LABEL here, not the word "lowest" (which 01
  # uses for --prs_quantile_reference). Passing "lowest" fell through to
  # strata_levels[1] and was correct only by coincidence. Q1 is the lowest
  # quantile produced by 01, so state it explicitly.
  CI_ARGS=(--models_file "$FM_MAIN" --time_points "$TIME_POINTS"
           --reference_quantile "${REFERENCE_QUANTILE:-Q1}"
           --prefix "$PREFIX" --score_role "$SCORE_ROLE" --role_quantiles "$ROLE_QUANTILES"
           --risk_model "$RISK_MODEL" --min_cell_count "$MIN_CELL_COUNT")
  [[ "$ENABLE_CLINICAL_UTILITY" == "1" ]] && CI_ARGS+=(--enable_clinical_utility)
  [[ "$BOOTSTRAP_RISK_DIFFERENCES" == "1" ]] && CI_ARGS+=(--bootstrap_risk_differences)
  run_stage 06 "$S06" "${CI_ARGS[@]}" $(uset) --outdir "$WORK_DIR"
fi

# ============================== Lag sensitivity sweep (Layer 3) ==============================
# The main run above is lag 0 (primary). Here we re-run STAGE 01 ONLY for each
# additional lag (lag sensitivity concerns the key tests, not the metric curves),
# then assemble a compact summary. Skipped if RUN_01=0 or LAG_DAYS has only 0.
if [[ "$RUN_01" == "1" ]]; then
  LAG_DIR="${WORK_DIR}/lag_sensitivity"
  IFS=',' read -ra LAG_ARR <<< "$LAG_DAYS"
  RAN_LAGS=()
  for L in "${LAG_ARR[@]}"; do
    L="$(echo "$L" | tr -d '[:space:]')"
    [[ "$L" == "0" || -z "$L" ]] && continue
    mkdir -p "$LAG_DIR"
    info "Lag sensitivity: refitting M0-M7 with washout > ${L} days (stage 01 only)"
    LAG_PREFIX="${BASE_STEM}_lag${L}"
    run_stage "01 (lag ${L}d)" "$S01" \
      --prs_file        "$ONSET_PRS_PATH" \
      --prs_col         "$ONSET_COL" \
      --pheno_file      "$PHENO_PATH" \
      --cov_file        "$COV_PATH" \
      --status_col      "$STATUS_COL" \
      "${AGE_TIME_ARGS[@]}" \
      "${MODE_ARGS[@]}" \
      --lag_days        "$L" \
      --prefix          "$LAG_PREFIX" \
      --covariates      "$COVARIATES" \
      --sex_col         "$SEX_COL" \
      --prs_quantiles   "$PRS_QUANTILES" \
      --time_points     "$TIME_POINTS" \
      --outdir          "$LAG_DIR" \
      "${MULTI_PRS_ARGS[@]}" \
      || echo "  ⚠ lag ${L} run failed (likely too few events after washout); continuing"
    RAN_LAGS+=("$L")
  done

  # Assemble lag_sensitivity_summary.csv: N, events, and key-test stats per lag.
  if [[ ${#RAN_LAGS[@]} -gt 0 ]]; then
    info "Assembling lag sensitivity summary"
    LAG_LIST="0,$(IFS=,; echo "${RAN_LAGS[*]}")"
    Rscript - "$WORK_DIR" "$LAG_DIR" "$BASE_STEM" "$LAG_LIST" <<'RSUM'
args <- commandArgs(trailingOnly=TRUE)
suppressPackageStartupMessages(library(data.table))
main_dir <- args[1]; lag_dir <- args[2]; stem <- args[3]
lags <- as.numeric(strsplit(args[4], ",")[[1]])
rows <- list()
for (L in lags) {
  if (L == 0) { d <- main_dir; pfx <- paste0(stem, "_01_fit_models_") }
  else        { d <- lag_dir;  pfx <- paste0(stem, "_lag", L, "_") }
  lad_f <- file.path(d, paste0(pfx, "model_lrt_ladder.csv"))
  cmp_f <- file.path(d, paste0(pfx, "model_comparison.csv"))
  if (!file.exists(lad_f)) next
  lad <- fread(lad_f)
  n <- NA_integer_; ev <- NA_integer_
  if (file.exists(cmp_f)) { cm <- fread(cmp_f); n <- max(cm$n, na.rm=TRUE); ev <- max(cm$n_events, na.rm=TRUE) }
  grab <- function(tst) { r <- lad[lad$test==tst]; if (nrow(r)) c(chisq=r$chisq[1], p=r$p[1], dC=r$delta_cindex[1]) else c(chisq=NA,p=NA,dC=NA) }
  m4 <- grab("M4_vs_M1"); m5 <- grab("M5_vs_M2"); m7 <- grab("M7_vs_M6")
  rows[[length(rows)+1]] <- data.table(lag_days=L, N=n, events=ev,
    M4vsM1_chisq=m4["chisq"], M4vsM1_p=m4["p"], M4vsM1_dCindex=m4["dC"],
    M5vsM2_chisq=m5["chisq"], M5vsM2_p=m5["p"], M5vsM2_dCindex=m5["dC"],
    M7vsM6_chisq=m7["chisq"], M7vsM6_p=m7["p"], M7vsM6_dCindex=m7["dC"])
}
if (length(rows)) {
  out <- rbindlist(rows)
  setorder(out, lag_days)
  # Write to WORK_DIR (intermediate), NOT final/. The publish step clears
  # final/${PREFIX}* before republishing from WORK_DIR, so a file written straight
  # to final/ here would be deleted and never restored. The allowlist entry
  # *lag_sensitivity_summary.csv then publishes it.
  fwrite(out, file.path(main_dir, paste0(stem, "_lag_sensitivity_summary.csv")))
  cat("✓ Saved lag_sensitivity_summary.csv (intermediate/ -> published)\n"); print(out)
} else cat("⚠ No lag ladders found to summarize\n")
RSUM
  fi
fi

# ============================== Publish shareable outputs ==============================
# Copy aggregate/summary files to final/. Keep individual-level and bulky working
# files in intermediate/ only:
#   - *.rds            (fitted models store per-person design matrices; data_processed is per-person)
#   - *risk_movement*  (per-individual predicted risks from stage 04)
#   - *load_models.R   (loader referencing local RDS paths)
#   - *.log            (verbose run logs)
info "Publishing shareable summaries to: $FINAL_DIR"

# Remove this run's previously published files before republishing, so a rerun
# that no longer produces some output cannot silently return the stale copy.
# Scoped to ${PREFIX}* so a sibling layer sharing --out-dir is untouched. Safe with
# --run-0N=0: WORK_DIR keeps every earlier stage's output and all are recopied.
shopt -s nullglob
stale=("$FINAL_DIR/${PREFIX}"*)
if (( ${#stale[@]} )); then
  rm -f "${stale[@]}"
  echo "  Cleared ${#stale[@]} previously published file(s) for ${PREFIX}*"
fi

# ALLOWLIST. A blacklist ("everything except *risk_movement*") silently ships any
# future per-individual output whose name nobody remembered to exclude. Only the
# aggregate summaries named here are publishable; anything else stays private and
# is reported below so a genuinely new aggregate output is noticed, not lost.
is_publishable() {
  case "$1" in
    # per-individual or environment-specific → never publish.
    # NOTE both *.log and *log.txt: stages write logs under both conventions, and
    # the loop below walks every file, not just .csv/.png/.pdf.
    *risk_movement*|*reclassification_table.csv|*.rds|*.log|*log.txt|*load_models.R) return 1 ;;   # individual-level or count cross-tabs -> private
  esac
  case "$1" in
    *model_comparison.csv|*model_lrt_ladder.csv|*all_coefficients.csv|*ph_tests.csv) return 0 ;;
    *run_status.csv|*run_status_discrimination.csv|*prs_vif_clinical.csv) return 0 ;;
    *prs_correlation.csv|*prs_vif.csv|*harrells_ci_details.csv|*uno_cindex_details.csv) return 0 ;;
    *lag_sensitivity_summary.csv|*timing_categories.csv|*t1_status_final.csv) return 0 ;;
    *final_analysis_summary.csv|*pre_alignment_summary.csv|*trajectory_summary.csv) return 0 ;;
    *source_phenotype_timing_qc.csv|*final_analysis_timing_qc.csv) return 0 ;;
    *trajectory_attrition.csv|*horizon_support.csv) return 0 ;;
    *time_summary.csv|*time_support.csv|*followup_reverseKM.csv) return 0 ;;
    *requested_timepoints_audit.csv|*requested_timepoints_ghat.csv) return 0 ;;
    *discrimination_summary.csv|*rr_auc_by_time.csv|*rr_brier_curve.csv|*tdauc_results.csv) return 0 ;;
    *rr_iauc_by_time.csv|*rr_ibs_by_time.csv|*rr_auc_contrasts.csv|*rr_brier_contrasts.csv) return 0 ;;
    *rr_iauc.csv|*rr_ibs.csv) return 0 ;;   # single-tau iAUC/IBS (some riskRegression versions)
    *ghat_censoring.csv|*cal_times_ghat_audit.csv) return 0 ;;
    *calibration_metrics.csv|*calibration_table*.csv|*calibration_by_time.csv) return 0 ;;
    *calibration_groups.csv|*calibration_plot_points.csv|*calibration_bands.csv|*calibration_plot_bands.csv) return 0 ;;
    *survival_summary*.csv|*rmst.csv|*median_times.csv|*strata_comparison_test_*.csv|*logrank_test*.csv) return 0 ;;
    *percentile_survival.csv|*logrank_trend.csv|*pairwise_comparisons.csv|*incidence_rates.csv|*stratified_results.csv) return 0 ;;
    *logrank_tests.csv|*percentile_survival.csv) return 0 ;;
    *absolute_risks_*.csv|*risk_differences_*.csv) return 0 ;;
    # both the per-stratification (risk_table_Q5.csv) and combined (risk_table.csv) forms
    *risk_table_*.csv|*risk_table.csv) return 0 ;;
    *risk_per_sd.csv|*risk_advancement_period.csv|*bootstrap_ci.csv) return 0 ;;
    # opt-in clinical-utility outputs (aggregate; only produced with the flag)
    *net_benefit.csv|*clinical_utility_metrics.csv|*population_attributable_fraction.csv) return 0 ;;
    *nns_*.csv|*max_observed_cumulative_risk.csv) return 0 ;;
    # opt-in reclassification aggregates (stage 04; only exist with --run-04=1)
    *idi_results.csv|*nri_continuous.csv|*nri_categorical.csv|*nri_*.csv) return 0 ;;
    *reclassification_times_ghat_audit.csv) return 0 ;;   # sanitized: count columns stripped on publish
    *session_manifest.csv|*run_configuration.csv|*downstream_rerun_manifest.csv) return 0 ;;
    # AGGREGATE FIGURES ONLY, by name. A blanket *.png|*.pdf rule would publish any
    # future per-individual scatterplot the moment some stage started drawing one.
    *km_*.png|*km_*.pdf|*survival_*.png|*survival_*.pdf) return 0 ;;
    *calibration_*.png|*calibration_*.pdf|*cal_*.png|*cal_*.pdf) return 0 ;;
    *auc_*.png|*auc_*.pdf|*pe_curve*.png|*pe_curve*.pdf) return 0 ;;
    *cumulative_incidence_*.png|*cumulative_incidence_*.pdf) return 0 ;;
    *decision_curve_*.png|*decision_curve_*.pdf|*decision_curve.png|*decision_curve.pdf) return 0 ;;
    *reclassification_plot*.png|*reclassification_plot*.pdf) return 0 ;;
    *forest*.png|*forest*.pdf|*slope_trajectory*.png|*slope_trajectory*.pdf) return 0 ;;
    *ici_trajectory*.png|*ici_trajectory*.pdf) return 0 ;;
  esac
  return 1
}

published=0; withheld=0; withheld_names=()
# Publish the count-bearing QC/audit tables in SANITIZED form: strip the raw
# n_events_by_t / n_at_risk columns (a failed horizon could otherwise report a
# count < the disclosure threshold). The detailed table stays in intermediate/.
# The QC signal (time, ghat, pass/fail, thresholds) is retained.
needs_sanitize() {
  case "$1" in
    *requested_timepoints_audit.csv|*time_support.csv|*ghat_censoring.csv|*cal_times_ghat_audit.csv) return 0 ;;
    *reclassification_times_ghat_audit.csv) return 0 ;;   # same n_events_by_t / n_at_risk columns
    # horizon_support: per-stratum counts (n, n_at_risk_at_end, n_events_by_end) let a hidden
    # small stratum be recovered from the published total in a two-stratum analysis. Strip them;
    # keep stratum / horizon / censoring_survival / max_baseline_time / supported / model.
    *horizon_support.csv) return 0 ;;
  esac
  return 1
}
# trajectory_attrition.csv leaks small per-stage removals by subtraction: N_removed is
# blanked in 01, but the consecutive N_remaining running totals let a suppressed removal
# be recovered as N_remaining[i-1] - N_remaining[i]. Publish a coarsened copy (the detailed
# table stays in intermediate/): always drop the derivable *_removed columns, and if ANY
# stage removes fewer than the disclosure threshold, blank the interior running totals so
# only the start (row 1) and final (last row) endpoints remain.
is_attrition() { case "$1" in *trajectory_attrition.csv) return 0 ;; esac; return 1; }
# lag_sensitivity_summary.csv: exact N/events at consecutive lags let a small lag-to-lag
# removal be recovered by subtraction (same class as attrition). If any consecutive delta
# is below the threshold, blank N/events for the lag>0 rows (keep lag 0) while retaining
# every test statistic. The detailed copy stays in intermediate/.
is_lagsummary() { case "$1" in *lag_sensitivity_summary.csv) return 0 ;; esac; return 1; }

for f in "$WORK_DIR"/*; do
  [[ -f "$f" ]] || continue
  base="$(basename "$f")"
  if is_publishable "$base"; then
    if needs_sanitize "$base"; then
      PRS_IN="$f" PRS_OUT="$FINAL_DIR/$base" Rscript - <<'SAN' && published=$((published+1))
suppressPackageStartupMessages(library(data.table))
d <- fread(Sys.getenv("PRS_IN"))
# Union of count columns across the sanitized QC/audit tables; intersect() drops only
# those present (n_events_by_t/n_at_risk for the time audits; n/n_at_risk_at_end/
# n_events_by_end for horizon_support). The remaining columns are the QC signal.
drop <- intersect(c("n_events_by_t","n_at_risk","n","n_at_risk_at_end","n_events_by_end"), names(d))
if (length(drop)) d[, (drop) := NULL]
fwrite(d, Sys.getenv("PRS_OUT"))
SAN
    elif is_attrition "$base"; then
      PRS_IN="$f" PRS_OUT="$FINAL_DIR/$base" MIN_CELL_COUNT="$MIN_CELL_COUNT" Rscript - <<'ATT' && published=$((published+1))
suppressPackageStartupMessages(library(data.table))
d <- fread(Sys.getenv("PRS_IN"))
thr <- suppressWarnings(as.integer(Sys.getenv("MIN_CELL_COUNT"))); if (is.na(thr)) thr <- 10L
drop <- intersect(c("N_removed","events_removed"), names(d))       # always drop: derivable + leak-prone
if (length(drop)) d[, (drop) := NULL]
n <- nrow(d)
if (n >= 3L) {                                                     # need interior rows to blank
  small <- FALSE
  for (col in intersect(c("N_remaining","events_remaining"), names(d))) {
    v <- suppressWarnings(as.numeric(d[[col]])); rem <- -diff(v)   # rem[i] = row i -> row i+1 removal
    if (any(is.finite(rem) & rem > 0 & rem < thr)) small <- TRUE
  }
  if (small) {
    ir <- 2:(n - 1L)                                              # interior rows; keep start + final
    for (col in intersect(c("N_remaining","events_remaining"), names(d)))
      d[ir, (col) := NA_integer_]
  }
}
fwrite(d, Sys.getenv("PRS_OUT"))
ATT
    elif is_lagsummary "$base"; then
      PRS_IN="$f" PRS_OUT="$FINAL_DIR/$base" MIN_CELL_COUNT="$MIN_CELL_COUNT" Rscript - <<'LAG' && published=$((published+1))
suppressPackageStartupMessages(library(data.table))
d <- fread(Sys.getenv("PRS_IN"))
thr <- suppressWarnings(as.integer(Sys.getenv("MIN_CELL_COUNT"))); if (is.na(thr)) thr <- 10L
n <- nrow(d)
if (n >= 2L) {                                       # rows are sorted by lag (lag 0 first)
  small <- FALSE
  for (col in intersect(c("N","events"), names(d))) {
    v <- suppressWarnings(as.numeric(d[[col]])); rem <- -diff(v)   # lag i -> lag i+1 removal
    if (any(is.finite(rem) & rem > 0 & rem < thr)) small <- TRUE
  }
  if (small) for (col in intersect(c("N","events"), names(d))) d[2:n, (col) := NA_integer_]  # keep lag 0
}
fwrite(d, Sys.getenv("PRS_OUT"))
LAG
    else
      cp -f "$f" "$FINAL_DIR/" && published=$((published+1))
    fi
  else
    case "$base" in
      *risk_movement*|*reclassification_table.csv|*ci_curves_*|*km_curve_data*|*.rds|*.log|*log.txt|*load_models.R) ;;   # expected private (individual-level / count cross-tabs / dense event-time), no need to report
      *) withheld=$((withheld+1)); withheld_names+=("$base") ;;
    esac
  fi
done
shopt -u nullglob
echo "  Published $published shareable file(s). Individual-level files remain in $WORK_DIR"
if (( withheld )); then
  echo "  ⚠ $withheld file(s) matched no allowlist entry and were NOT published:"
  printf '      %s\n' "${withheld_names[@]}"
  echo "    If any is an aggregate summary, add it to is_publishable() in this script."
fi

# Consolidated completeness check from the machine-readable status files. Flag
# only (no hard fail) so a site notices an incomplete run before returning it.
RS="${FINAL_DIR}/${PREFIX}_01_fit_models_run_status.csv"
RSD="${FINAL_DIR}/${PREFIX}_03_discrimination_metrics_run_status_discrimination.csv"
if [[ -f "$RS" ]]; then
  Rscript - "$RS" "$RSD" <<'RSTAT'
args <- commandArgs(trailingOnly=TRUE)
suppressPackageStartupMessages(library(data.table))
g <- function(d,c) if (c %in% names(d)) as.logical(d[[c]][1]) else NA
rs <- tryCatch(fread(args[1]), error=function(e) NULL)
rsd <- if (length(args) > 1 && nzchar(args[2]) && file.exists(args[2])) tryCatch(fread(args[2]), error=function(e) NULL) else NULL
ck <- function(x) if (isTRUE(x)) "✓" else if (isFALSE(x)) "✗" else "?"
core <- g(rs,"core_M0_M7_complete"); m4 <- g(rs,"M4_vs_M1_available")
m5 <- g(rs,"M5_vs_M2_available"); m7 <- g(rs,"M7_vs_M6_available")
dauc <- g(rsd,"paired_delta_AUC_available"); dbri <- g(rsd,"paired_delta_Brier_available")
cat(sprintf("  Completeness: core_M0_M7 %s  M4vsM1 %s  M5vsM2 %s  M7vsM6 %s  dAUC %s  dBrier %s\n",
            ck(core), ck(m4), ck(m5), ck(m7), ck(dauc), ck(dbri)))
complete <- isTRUE(core) && isTRUE(m4) && isTRUE(m5) && isTRUE(m7) &&
            (is.null(rsd) || (isTRUE(dauc) && isTRUE(dbri)))
if (complete) cat("  RUN COMPLETE: yes\n") else {
  cat("  RUN COMPLETE: NO\n")
  cat("  ⚠ One or more primary tests / paired contrasts are missing — review before returning final/.\n")
}
RSTAT
fi

info "Pipeline finished: ${COHORT}_${ANCESTRY}_${TRAIT}  (mode: ${ANALYSIS_MODE})"
