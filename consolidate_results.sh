#!/usr/bin/env bash
# =============================================================================
# consolidate_results.sh — bundle a biobank's shareable PRS-metrics outputs into
# ONE file for a single sandbox export.
#
# The pipeline writes each run to  <out-dir>/final/  (aggregate, disclosure-safe:
# allowlisted + small-cell-suppressed) and  <out-dir>/intermediate/  (private,
# individual-level). This script collects ONLY the  final/  contents from every
# run under a results root and packages them as a single artifact, so a federated
# site (e.g. FinnGen) that manually reviews each export needs to shepherd through
# just one file instead of hundreds of CSVs.
#
# It NEVER reads  intermediate/ . It also scans every bundled CSV header for a
# genuine per-individual ID column and refuses to write if one is found.
#
# Usage:
#   bash consolidate_results.sh --results-root=DIR [options]
#
# Options:
#   --results-root=DIR   Directory searched recursively for run  final/  folders (required).
#   --out=PATH           Output artifact path. Default: <results-root>/consolidated_prs_results.<ext>
#   --format=tar|text    tar  = one .tar.gz (default; smallest, preserves each CSV).
#                        text = one plain-text .txt with every table concatenated and
#                               delimited (for sandboxes that only export/review text).
#   --cohort=NAME        Label recorded in the README / default filename (e.g. FinnGen).
#   --include-figures    Also include the .png figures (excluded by default; they are
#                        derivable at the coordinating centre from the CSVs).
#   --allow-flagged      Package even if the ID-column scan finds something (records it
#                        in the manifest). Off by default — the scan is a hard stop.
#   -h | --help
#
# Output is a SINGLE file that contains, at its root:
#   MANIFEST.tsv  — every table with sha256, size and row count
#   README.txt    — what this is, when it was made, counts, and the disclosure-scan result
#   <run>/...     — the per-run final/ tables, grouped by run
# =============================================================================
set -euo pipefail

RESULTS_ROOT=""; OUT=""; FORMAT="tar"; COHORT=""; INCLUDE_FIGURES=0; ALLOW_FLAGGED=0
die() { echo "ERROR: $*" >&2; exit 1; }

for arg in "$@"; do
  case "$arg" in
    --results-root=*) RESULTS_ROOT="${arg#*=}" ;;
    --out=*)          OUT="${arg#*=}" ;;
    --format=*)       FORMAT="${arg#*=}" ;;
    --cohort=*)       COHORT="${arg#*=}" ;;
    --include-figures) INCLUDE_FIGURES=1 ;;
    --allow-flagged)  ALLOW_FLAGGED=1 ;;
    -h|--help) sed -n '2,45p' "$0"; exit 0 ;;
    *) die "unknown option: $arg (see --help)" ;;
  esac
done

[[ -n "$RESULTS_ROOT" ]] || die "--results-root is required (see --help)."
[[ -d "$RESULTS_ROOT" ]] || die "results root not found: $RESULTS_ROOT"
case "$FORMAT" in tar|text) ;; *) die "--format must be 'tar' or 'text' (got '$FORMAT')." ;; esac
RESULTS_ROOT="$(cd "$RESULTS_ROOT" && pwd)"

# Default output path/name.
if [[ -z "$OUT" ]]; then
  base="consolidated_prs_results${COHORT:+_$COHORT}"
  OUT="$RESULTS_ROOT/$base.$([[ "$FORMAT" == "tar" ]] && echo tar.gz || echo txt)"
fi

# ---- Find every run's final/ directory (NEVER intermediate/) -----------------
finals=()
while IFS= read -r d; do finals+=("$d"); done < <(find "$RESULTS_ROOT" -type d -name final -not -path '*/intermediate/*' | sort)
[[ ${#finals[@]} -gt 0 ]] || die "no  final/  directories found under $RESULTS_ROOT. Point --results-root at the folder that holds your run output directories."

echo "Found ${#finals[@]} run(s) with a final/ directory under $RESULTS_ROOT"

STAGE="$(mktemp -d 2>/dev/null || echo "/tmp/prs_consolidate_$$")"
ROOTNAME="consolidated_prs_results"
DEST="$STAGE/$ROOTNAME"
mkdir -p "$DEST"
trap 'rm -rf "$STAGE"' EXIT

MAN="$DEST/MANIFEST.tsv"
printf 'run\tfile\tsha256\tbytes\tdata_rows\n' > "$MAN"

ID_RE='^(IID|FID|eid|person_id|patient_id)$'   # genuine per-individual identifiers
flagged=()
n_files=0; n_csv=0; n_png=0

for fdir in "${finals[@]}"; do
  rel="${fdir#"$RESULTS_ROOT"/}"; run="${rel%/final}"
  [[ "$run" == "final" || -z "$run" ]] && run="$(basename "$RESULTS_ROOT")"
  mkdir -p "$DEST/$run"
  for f in "$fdir"/*; do
    [[ -f "$f" ]] || continue
    b="$(basename "$f")"
    case "$b" in
      *.png|*.pdf) [[ "$INCLUDE_FIGURES" == "1" ]] || continue; n_png=$((n_png+1)) ;;
      *.csv|*.tsv) n_csv=$((n_csv+1))
        # Disclosure scan: reject a genuine per-individual ID column.
        if head -1 "$f" | tr ',\t' '\n\n' | sed 's/^"//; s/"$//' | grep -qiE "$ID_RE"; then
          flagged+=("$run/$b")
        fi ;;
      *) ;;   # skip anything unexpected
    esac
    cp -p "$f" "$DEST/$run/$b"
    sha="$(shasum -a 256 "$f" | awk '{print $1}')"
    bytes="$(wc -c < "$f" | tr -d ' ')"
    rows=""; case "$b" in *.csv|*.tsv) rows="$(awk 'END{print NR-1}' "$f")" ;; esac
    printf '%s\t%s\t%s\t%s\t%s\n' "$run" "$b" "$sha" "$bytes" "$rows" >> "$MAN"
    n_files=$((n_files+1))
  done
done

# ---- Disclosure scan verdict -------------------------------------------------
if [[ ${#flagged[@]} -gt 0 ]]; then
  echo "" >&2
  echo "⚠ DISCLOSURE SCAN: ${#flagged[@]} bundled file(s) contain a per-individual ID column:" >&2
  printf '    %s\n' "${flagged[@]}" >&2
  if [[ "$ALLOW_FLAGGED" != "1" ]]; then
    die "refusing to write. These belong in intermediate/, not final/ — check your run, or re-run --allow-flagged if you have verified they are aggregate."
  fi
  echo "  (--allow-flagged set: packaging anyway; recorded in README.)" >&2
fi

# ---- README ------------------------------------------------------------------
{
  echo "GBMI Longitudinal-PRS — consolidated shareable results"
  echo "Generated: $(date -u '+%Y-%m-%d %H:%M:%S UTC')"
  [[ -n "$COHORT" ]] && echo "Cohort: $COHORT"
  echo "Runs bundled: ${#finals[@]}   Files: $n_files (CSV/TSV: $n_csv, figures: $n_png)"
  echo "Source: the final/ (aggregate, allowlisted + small-cell-suppressed) outputs only."
  echo "intermediate/ (individual-level) content was NOT read."
  echo "Per-individual ID-column scan: ${#flagged[@]} file(s) flagged$([[ ${#flagged[@]} -gt 0 && "$ALLOW_FLAGGED" == 1 ]] && echo ' (packaged with --allow-flagged)')."
  echo "See MANIFEST.tsv for per-file sha256, size and row counts."
  [[ "$FORMAT" == "text" ]] && echo "Split this text bundle at lines matching '===== BEGIN FILE: <run>/<name> ====='."
} > "$DEST/README.txt"

# ---- Package into a SINGLE file ---------------------------------------------
mkdir -p "$(dirname "$OUT")"
if [[ "$FORMAT" == "tar" ]]; then
  tar -czf "$OUT" -C "$STAGE" "$ROOTNAME"
else
  {
    echo "########## $ROOTNAME (concatenated bundle) ##########"; echo
    for meta in README.txt MANIFEST.tsv; do
      echo "===== BEGIN FILE: $meta ====="; cat "$DEST/$meta"; echo "===== END FILE: $meta ====="; echo
    done
    # per-run tables
    ( cd "$DEST" && find . -type f \( -name '*.csv' -o -name '*.tsv' \) | sort ) | while IFS= read -r rf; do
      rf="${rf#./}"
      echo "===== BEGIN FILE: $rf ====="; cat "$DEST/$rf"; echo "===== END FILE: $rf ====="; echo
    done
  } > "$OUT"
fi

OUT_SHA="$(shasum -a 256 "$OUT" | awk '{print $1}')"
OUT_SZ="$(du -h "$OUT" | awk '{print $1}')"
echo ""
echo "✅ Wrote ONE export artifact:"
echo "   $OUT"
echo "   size: $OUT_SZ   sha256: $OUT_SHA"
echo "   ${#finals[@]} run(s), $n_files file(s). Export this single file from the sandbox."
