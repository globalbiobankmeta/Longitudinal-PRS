# Calculate PRS metrics for disease-progression trajectories

This pipeline evaluates polygenic risk scores derived from time-to-event GWAS across clinically ordered disease trajectories:

`T0 → T1 → T2`

For example, in the `T2D_to_CAD` trajectory, T1 is type 2 diabetes and T2 is coronary artery disease.

The pipeline uses three PRS:

- **Onset PRS (`PRS_T0T1`)**: genetic risk of reaching T1
- **Outcome PRS (`PRS_T0T2`)**: genetic risk of reaching T2 in the broader population
- **Progression PRS (`PRS_T1T2`)**: genetic risk of progressing from T1 to T2 among people with T1


<a id="top"></a>

<details open>
<summary><strong>Contents</strong></summary>

- [1. Quick start](#section-1)
- [2. Models](#section-2)
- [3. Prerequisites](#section-3)
  - [3.1 Software](#section-3-1)
  - [3.2 Download the pipeline](#section-3-2)
  - [3.3 Test the installation](#section-3-3)
- [4. Input files](#section-4)
  - [4.1 Three PRS files](#section-4-1)
  - [4.2 Progression phenotype file](#section-4-2)
  - [4.3 Covariate file](#section-4-3)
  - [4.4 Recruitment-age file](#section-4-4)
  - [4.5 Population keep-list (optional)](#section-4-5)
- [5. Set variables](#section-5)
- [6. Run the pipeline](#section-6)
  - [6.1 Layer 1: GWAS-aligned analysis](#section-6-1)
  - [6.2 Layer 2: prospective validation](#section-6-2)
- [7. Confirm the run succeeded](#section-7)
- [8. What to share](#section-8)
- [9. Troubleshooting](#section-9)
- [10. Important interpretation note](#section-10)

</details>

---

<a id="section-1"></a>

# 1. Quick start

For each assigned **trajectory × target ancestry**:

1. Run **Layer 1** for the GWAS-aligned analysis.
2. Run **Layer 2** for prospective validation when recruitment age is available.
3. Run one trajectory first as a pilot.
4. Confirm all input files have a unique `IID`.
5. Run `bash smoke_test.sh` once after installing the pipeline.
6. Edit the variables in Section 5.
7. Put the variable block and pipeline command into your site's job script, or run them within an allocated compute-node session.
8. Confirm the job ends with:

   ```text
   RUN COMPLETE: yes
   ```

9. After local disclosure review, share the corresponding **`final/` folder** with the coordinating team.
10. Do **not** share `intermediate/`, phenotype files, individual-level data, or local log files.

**Only full trajectories using all three PRS are required.**

---

<a id="section-2"></a>

# 2. Models

The pipeline fits the following models:

| Model | PRSs |
|---|---|
| M0 | Covariates only |
| M1 | Onset PRS |
| M2 | Outcome PRS |
| M3 | Progression PRS |
| M4 | Onset + progression PRS |
| M5 | Outcome + progression PRS |
| M6 | Onset + outcome PRS |
| M7 | All three PRS |

The main comparisons are:

- **M4 vs M1**: progression PRS beyond onset PRS
- **M5 vs M2**: progression PRS beyond outcome PRS
- **M7 vs M6**: progression PRS beyond both onset and outcome PRS

---

<a id="section-3"></a>

# 3. Prerequisites

<a id="section-3-1"></a>

## 3.1 Software

- R ≥ 4.1
- Bash
- Access to your biobank's batch-compute environment or an allocated compute node
- Required R packages:

```r
install.packages(c(
  "optparse", "data.table", "survival", "riskRegression",
  "ggplot2", "survminer", "scales", "RColorBrewer",
  "rms", "pec", "boot", "rlang", "gridExtra",
  "ggpubr", "timeROC", "survAUC", "prodlim"
))
```

The pipeline checks package availability and minimum versions before running.

<a id="section-3-2"></a>

## 3.2 Download the pipeline

The GBMI longitudinal PRS github repository: https://github.com/globalbiobankmeta/Longitudinal-PRS

```bash
git clone https://github.com/globalbiobankmeta/Longitudinal-PRS.git
cd Longitudinal-PRS
```

Keep all scripts together in the same directory:

```text
00_run_pipeline_gbmi.sh
01_fit_models.R
02_kaplan_meier.R
03_discrimination_metrics.R
04_reclassification.R
05_calibration.R
06_cumulative_incidence.R
prs_risk_utils.R
test_prs_risk_utils.R
smoke_test.sh
```

Use the exact filenames shown above.

<a id="section-3-3"></a>

## 3.3 Test the installation

Run the installation test once after installing the pipeline. Run it in an allocated compute environment when required by your site:

```bash
bash smoke_test.sh
```

A successful test ends with:

```text
SMOKE TEST PASSED
```

When the smoke test fails, run the helper test for additional diagnosis:

```bash
Rscript --vanilla test_prs_risk_utils.R
```

Repeat the smoke test after updating the pipeline, R, or major package versions.

---

<a id="section-4"></a>

# 4. Input files

<a id="section-4-1"></a>

## 4.1 Three PRS files

The weights for per-biobank using LOBO multi-ancestry meta GWAS will be distributed to each biobank analyst.

All three PRS files are required:

| PRS | Description |
|---|---|
| `PRS_T0T1` | Onset PRS |
| `PRS_T1T2` | Progression PRS |
| `PRS_T0T2` | Outcome PRS |

Each file must contain:

- `IID`
- a PRS score column, usually `SCORE1_SUM`

Individuals missing any one of the three scores are excluded from all models.



<span style="color:red;">**Note: Use unrelated individuals for PRS calculation.**</span>

The **survival models likewise assume independent observations**, so the analysis should also be run
on unrelated individuals — see [4.5 Population keep-list](#section-4-5) for the optional `--pop-file`
that enforces this (and/or restricts to one target ancestry).

For PRS calculation and chromosome merging, see:

- https://github.com/globalbiobankmeta/PRS
- https://github.com/globalbiobankmeta/PRS/blob/main/run_prscs_auto_pipe.md

<span style="color:red;">
<strong>Note:</strong> To ensure compatibility with the pipeline, please use the following naming conventions for PRS outputs: <code>T1</code> for onset PRS, <code>T1toT2</code> or <code>T1_to_T2</code> for progression PRS, and <code>T2</code> for outcome PRS. Please refer to the <a href="https://github.com/globalbiobankmeta/Longitudinal-PRS/blob/main/gbmi_trajectory_biobank_gwas_names.tsv">trajectory definition file</a> for the complete trajectory definitions for each biobank.
</span>

<a id="section-4-2"></a>

## 4.2 Progression phenotype file

The phenotype file should contain one row per person with T1.

Required columns for the recommended age-derived analysis:

| Column | Description |
|---|---|
| `IID` | Participant identifier |
| `secondEvent` | 1 if T2 occurred after T1, otherwise 0 |
| `diagAge` | Age at T1 |
| `ageExit` | Age at T2 for events, or age at censoring for non-events |
| `sex` | Binary numeric variable coded 0/1 |
| `PC1`-`PC10` | Genetic principal components |
| `birthyear` | Birth year |

The recommended approach calculates:

```text
time since T1 = ageExit - diagAge
```

A validated precomputed time-since-T1 column may instead be supplied with `--time-exit`.

Use exactly one of these time definitions:

```bash
--age-t1-col=diagAge --age-exit-col=ageExit
```

or:

```bash
--age-t1-col=diagAge --time-exit=secondTime
```

The `--time-exit` column must already contain elapsed time from T1 to T2 or censoring, in years. It must not contain absolute age.

### Phenotype requirements

- `IID` must be unique.
- Age and follow-up variables must be numeric and measured in years.
- Fractional age is preferred.
- For events, use the first qualifying T2 after T1.
- For non-events, use the agreed consortium censoring definition.
- Use `NA` or an empty value for missing data.
- Recode sex to numeric 0/1.
- Covariate names must match the names supplied to `--covariates`.

Example:

```text
IID    secondEvent    diagAge    ageExit    sex    PC1      PC2      birthyear
1001   1              52.41      61.27      0      0.0031   -0.012   1962
1002   0              48.10      65.32      1     -0.0088    0.004   1958
1003   1              67.03      72.14      1      0.0015   -0.021   1949
```

<a id="section-4-3"></a>

## 4.3 Covariate file

The covariate file may be the same as the phenotype file. 

**Recommended: Use the same covariates as used in GWAS.**

Required columns:

- `IID`
- `sex`
- genetic PCs
- `birthyear`
- optional `batch` or other biobanks-specific covariates

Do not include age at T1 in `--covariates`; it is added automatically.

<a id="section-4-4"></a>

## 4.4 Recruitment-age file

Required for Layer 2 and recommended for Layer 1.

It should contain:

- `IID`
- age at recruitment or consent

Example:

```text
IID    age_at_recruitment
1001   54.2
1002   51.6
1003   63.8
```

---

<a id="section-4-5"></a>

## 4.5 Population keep-list (optional)

An **optional** input file that restricts the analysis to a chosen set of individuals. Leave it unset
(`POP_FILE=""`) to use the whole phenotype file — most sites, who already supply a per-ancestry,
unrelated phenotype, do **not** need it.

Use it when your phenotype file still contains individuals you want to drop, most commonly to keep only:

- **unrelated individuals** — the survival models assume independent observations (the same reason
  unrelated individuals are required for PRS calculation); and/or
- a single **target ancestry**, when you run from a shared multi-ancestry input.

Format: a **headerless PLINK keep-list** with two whitespace-separated columns, `FID IID`; matching is
on the second column (IID). Example:

```text
1001   1001
1002   1002
1005   1005
```

It is an **IID subset only** — it does not change your covariates or scores, so you are still
responsible for supplying ancestry-appropriate covariates/PCs and PRS. Pass it with `--pop-file`
(variable `POP_FILE` in Section 5); it is applied in stage 01 before every downstream stage, so both
Layer 1 (including the lag re-fits) and Layer 2 (prospective re-indexing) run on the restricted sample.

---

<a id="section-5"></a>

# 5. Set variables

Copy this block into the batch job script for each trajectory and target ancestry. Edit the values in this section rather than changing the pipeline commands in Section 6.

```bash
#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------------
# Trajectory, cohort, ancestry and PRS method
# ------------------------------------------------------------------
TRAIT="T1D_to_CAD"      #Update based on your biobank progression phenotype name
COHORT="MYBIOBANK"      #e.g, UKBB
ANCESTRY="EUR"          # Target ancestry

PRS_METHOD="PRScs"      #PRScs or PT
DISCOVERY_ANCESTRY="MULTI"

# ------------------------------------------------------------------
# Pipeline and output locations
# ------------------------------------------------------------------
SCRIPT_DIR="/path/to/Longitudinal-PRS"  # where you downloaded the github repo

OUT_ROOT="/path/to/results/${TRAIT}/${ANCESTRY}/${PRS_METHOD}" #Update to reflect your output directory
LAYER1_OUT_DIR="${OUT_ROOT}/gwas_aligned" #default
LAYER2_OUT_DIR="${OUT_ROOT}/prospective" #default

# ------------------------------------------------------------------
# Input files
# ------------------------------------------------------------------

# --- Derive onset (T1) and outcome (T2) names from TRAIT (split on 'to' / '_to_') ---
# Assumes your PRS .sscore files use the standard trajectory naming. If your files use
# other local names (e.g. AOU Birth_to_<T1>, FinnGen endpoint codes, phecodes), just edit
# the three *_PRS paths below directly.

read -r ONSET_NAME OUTCOME_NAME <<< "$(awk -v t="$TRAIT" \
  'BEGIN{ s=(t ~ /_to_/) ? "_to_" : "to"; n=split(t,a,s); print a[1], a[n] }')"
if [[ -z "$ONSET_NAME" || -z "$OUTCOME_NAME" || "$ONSET_NAME" == "$TRAIT" ]]; then
  echo "ERROR: TRAIT must contain 'to' or '_to_' (e.g. T2DtoCAD or T2D_to_CAD)." >&2
  exit 1
fi

ONSET_PRS="/path/to/prs/${ONSET_NAME}.sscore" ##note to update prefix or suffix for the onset PRS file (T0T1)
PROG_PRS="/path/to/prs/${TRAIT}.sscore" ##note to update prefix or suffix for the progression PRS file (T1T2)
OUTCOME_PRS="/path/to/prs/${OUTCOME_NAME}.sscore" ##note to update prefix or suffix for the outcome PRS file (T0T2)

PHENO="/path/to/pheno/pheno.txt" # full path to your phenotype file

# Leave empty when the covariates are stored in PHENO.
COV_FILE="" # full path to your covariate file

# Required for Layer 2 and recommended for Layer 1.
# Leave empty only when Layer 1 will be run without recruitment-timing QC.
RECRUIT="/path/to/recruit/recruit.txt" #note to update as empty when there is no recruitement age file in your biobank

# Optional keep-list (headerless PLINK "FID IID"); restricts the run to these individuals —
# unrelated set and/or one target ancestry. See section 4.5. Leave empty to use the whole input.
POP_FILE="" 

# ------------------------------------------------------------------------------
# Input column names [update to reflect the names used in your biobank]
# ------------------------------------------------------------------------------
STATUS_COL="secondEvent"            #status of T2 (0/1)
SEX_COL="sex"                     
PRS_SCORE_COL="SCORE1_SUM"
RECRUIT_AGE_COL="age_at_recruitment" #note this is age note Date 

COVARIATES="sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,birthyear"  #same as GWAS; don't include age at T1 here

# ------------------------------------------------------------------
# Time definition: [choose exactly one approach]
# ------------------------------------------------------------------
AGE_T1_COL="diagAge" #age at T1

# Recommended: derive time since T1 from ageExit - diagAge.
AGE_EXIT_COL="ageExit" #age at T2 or censored
TIME_EXIT_COL=""  

# Alternative: use a validated precomputed time-since-T1 column.
# To use this approach, comment out the two lines above and use:
# AGE_EXIT_COL=""
# TIME_EXIT_COL="secondTime" #time since T1 (in years)

if [[ -n "$AGE_EXIT_COL" && -z "$TIME_EXIT_COL" ]]; then
  TIME_ARGS=(
    --age-t1-col="$AGE_T1_COL"
    --age-exit-col="$AGE_EXIT_COL"
    --age-spline-df="$AGE_SPLINE_DF"
  )
elif [[ -z "$AGE_EXIT_COL" && -n "$TIME_EXIT_COL" ]]; then
  TIME_ARGS=(
    --age-t1-col="$AGE_T1_COL"
    --time-exit="$TIME_EXIT_COL"
    --age-spline-df="$AGE_SPLINE_DF"
  )
else
  echo "ERROR: set exactly one of AGE_EXIT_COL or TIME_EXIT_COL." >&2
  exit 1
fi

# ------------------------------------------------------------------
# Analysis layers and lag settings [keep as default]
# ------------------------------------------------------------------
LAYER1_ANALYSIS_MODE="gwas_aligned"
LAYER2_ANALYSIS_MODE="prospective"

# Layer 1 washout sensitivity.
LAYER1_LAG_DAYS="0,30,90,365"

# Keep this at 0 for Layer 2. Prospective mode defines its own prediction
# index and must not be combined with the Layer-1 lag sweep.
LAYER2_LAG_DAYS="0"

# Optional delay applied only to incident T1 cases in Layer 2.
# Keep 0 for the primary prospective analysis.
INCIDENT_LAG_DAYS="0"

# ------------------------------------------------------------------
# Evaluation and disclosure settings [keep as default]
# ------------------------------------------------------------------
TIME_POINTS="1,2,5,10,15,20,25,30,35"

MIN_EVENTS_TOTAL="50"
MIN_CELL_COUNT="10"

PRS_QUANTILES="5,10"
PRS_EXTREMES="1,2,5,10,20"

USE_SUPPORTED_TIMES="1"
RUN_TESTS="1"

# ------------------------------------------------------------------
# Optional argument arrays [keep as default]
# ------------------------------------------------------------------
COV_FILE_ARGS=()
if [[ -n "$COV_FILE" ]]; then
  COV_FILE_ARGS=(--cov-file="$COV_FILE")
fi

POP_FILE_ARGS=()
if [[ -n "$POP_FILE" ]]; then
  POP_FILE_ARGS=(--pop-file="$POP_FILE")
fi

LAYER1_RECRUIT_ARGS=()
if [[ -n "$RECRUIT" ]]; then
  LAYER1_RECRUIT_ARGS=(
    --recruit-file="$RECRUIT"
    --recruit-age-col="$RECRUIT_AGE_COL"
  )
fi
```

For the precomputed-time approach, set:

```bash
AGE_EXIT_COL=""
TIME_EXIT_COL="secondTime"
```

Do not set both.

---

<a id="section-6"></a>

# 6. Run the pipeline

<a id="section-6-1"></a>

## 6.1 Layer 1: GWAS-aligned analysis

```bash
bash "$SCRIPT_DIR/00_run_pipeline_gbmi.sh" \
  --trait="$TRAIT" \
  --out-dir="$LAYER1_OUT_DIR" \
  --script-dir="$SCRIPT_DIR" \
  --onset-prs-file="$ONSET_PRS" \
  --progression-prs-file="$PROG_PRS" \
  --outcome-prs-file="$OUTCOME_PRS" \
  --pheno-file="$PHENO" \
  ${COV_FILE_ARGS[@]+"${COV_FILE_ARGS[@]}"} \
  ${POP_FILE_ARGS[@]+"${POP_FILE_ARGS[@]}"} \
  --cohort="$COHORT" \
  --ancestry="$ANCESTRY" \
  --status-col="$STATUS_COL" \
  --sex-col="$SEX_COL" \
  --prs-score-col="$PRS_SCORE_COL" \
  --covariates="$COVARIATES" \
  ${TIME_ARGS[@]+"${TIME_ARGS[@]}"} \
  ${LAYER1_RECRUIT_ARGS[@]+"${LAYER1_RECRUIT_ARGS[@]}"} \
  --analysis-mode="$LAYER1_ANALYSIS_MODE" \
  --prs-method="$PRS_METHOD" \
  --discovery-ancestry="$DISCOVERY_ANCESTRY" \
  --time-points="$TIME_POINTS" \
  --lag-days="$LAYER1_LAG_DAYS" \
  --min-events-total="$MIN_EVENTS_TOTAL" \
  --min-cell-count="$MIN_CELL_COUNT" \
  --prs-quantiles="$PRS_QUANTILES" \
  --prs-extremes="$PRS_EXTREMES" \
  --use_supported_times="$USE_SUPPORTED_TIMES" \
  --run-tests="$RUN_TESTS"
```

Layer 1 is the primary GWAS-aligned analysis and includes the requested washout sensitivity analyses.


<a id="section-6-2"></a>

## 6.2 Layer 2: prospective validation

Layer 2 requires the recruitment-age file:

```bash
if [[ -z "$RECRUIT" ]]; then
  echo "ERROR: RECRUIT is required for Layer 2." >&2
  exit 1
fi
```

Then run:

```bash
bash "$SCRIPT_DIR/00_run_pipeline_gbmi.sh" \
  --trait="$TRAIT" \
  --out-dir="$LAYER2_OUT_DIR" \
  --script-dir="$SCRIPT_DIR" \
  --onset-prs-file="$ONSET_PRS" \
  --progression-prs-file="$PROG_PRS" \
  --outcome-prs-file="$OUTCOME_PRS" \
  --pheno-file="$PHENO" \
  ${COV_FILE_ARGS[@]+"${COV_FILE_ARGS[@]}"} \
  ${POP_FILE_ARGS[@]+"${POP_FILE_ARGS[@]}"} \
  --cohort="$COHORT" \
  --ancestry="$ANCESTRY" \
  --status-col="$STATUS_COL" \
  --sex-col="$SEX_COL" \
  --prs-score-col="$PRS_SCORE_COL" \
  --covariates="$COVARIATES" \
  ${TIME_ARGS[@]+"${TIME_ARGS[@]}"} \
  --analysis-mode="$LAYER2_ANALYSIS_MODE" \
  --recruit-file="$RECRUIT" \
  --recruit-age-col="$RECRUIT_AGE_COL" \
  --lag-days="$LAYER2_LAG_DAYS" \
  --incident-lag-days="$INCIDENT_LAG_DAYS" \
  --prs-method="$PRS_METHOD" \
  --discovery-ancestry="$DISCOVERY_ANCESTRY" \
  --time-points="$TIME_POINTS" \
  --min-events-total="$MIN_EVENTS_TOTAL" \
  --min-cell-count="$MIN_CELL_COUNT" \
  --prs-quantiles="$PRS_QUANTILES" \
  --prs-extremes="$PRS_EXTREMES" \
  --use_supported_times="$USE_SUPPORTED_TIMES" \
  --run-tests="$RUN_TESTS"
```

Layer 1 and Layer 2 must be run separately and returned as separate folders.

When recruitment age is unavailable, run Layer 1 only and notify the coordinating team.

---

<a id="section-7"></a>

# 7. Confirm the run succeeded

Before sharing results, confirm:

- the compute job exit status is 0;
- the job output ends with `RUN COMPLETE: yes`;
- the `final/` folder exists;
- no error is reported in the job output;
- the requested Layer 1 lag analyses completed;
- the local disclosure review approves the contents of `final/`.

When a job fails, keep its output directory and scheduler/job log for troubleshooting. Do not send phenotype or individual-level files by email.

---

<a id="section-8"></a>

# 8. What to share

Each run creates:

```text
final/
intermediate/
```

Share only:

```text
final/
```

Keep private:

```text
intermediate/
```

For Layer 1, share:

```text
$LAYER1_OUT_DIR/final/
```

For Layer 2, share:

```text
$LAYER2_OUT_DIR/final/
```

Send the two folders separately to the coordinating team after local disclosure review.

Do not share:

- `intermediate/`
- phenotype files
- PRS files
- participant-level data
- individual-level predictions
- local scheduler or analysis logs unless specifically requested for troubleshooting

The `final/` folder contains the aggregate outputs and run information needed by the coordinating team.

---

<a id="section-9"></a>

# 9. Troubleshooting

| Problem | Suggested action |
|---|---|
| Smoke test fails | Run `Rscript --vanilla test_prs_risk_utils.R` and review the first error |
| Missing R package | Install the package named in the error and rerun |
| Duplicate `IID` error | De-duplicate the corresponding input file |
| Missing PRS score column | Check `PRS_SCORE_COL` |
| Missing covariate | Confirm the column exists and is listed correctly in `COVARIATES` |
| Time-definition error | Set exactly one of `AGE_EXIT_COL` or `TIME_EXIT_COL` |
| Too few events | Notify the coordinating team that the trajectory cannot be analyzed |
| Many people are dropped | Check overlap among the three PRS files and the phenotype file |
| Layer 2 fails | Confirm `RECRUIT` and `RECRUIT_AGE_COL` |
| Path with spaces fails | Quote all paths; when needed, move the pipeline to a path without special characters |
| `RUN COMPLETE: no` | Keep the job output and contact the coordinating team |

---

<a id="section-10"></a>

# 10. Important interpretation note

Layer 1 and Layer 2 answer different questions:

- **Layer 1** evaluates association with recorded time from T1 to T2.
- **Layer 2** evaluates prediction of future T2 from a defined prediction index.

Do not combine or pool Layer 1 and Layer 2 results.

Please contact Ying Wang (yiwang@broadinstitute.org) if you have any questions, and share the results with her upon completion.
