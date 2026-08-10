#!/usr/bin/env bash
# =============================================================================
# extract_results.sh — unpack a consolidated export back into the individual
# result files, and verify integrity.
#
# Companion to consolidate_results.sh. Accepts EITHER artifact it produces:
#   * a .tar.gz archive, or
#   * a single concatenated plain-text .txt bundle
# and restores the per-run/per-file tree:
#   <out-dir>/<run>/<file>.csv ...   plus MANIFEST.tsv and README.txt
# Then it re-checksums every extracted table against MANIFEST.tsv and reports
# any mismatch or missing file — so the coordinating centre knows the export
# arrived intact before running meta-analysis.
#
# Usage:
#   bash extract_results.sh --input=FILE [--out-dir=DIR] [--no-verify]
#
# Options:
#   --input=FILE   The consolidated .tar.gz or .txt to unpack (required).
#   --out-dir=DIR  Where to restore the files. Default: <input without extension>_extracted/
#   --no-verify    Skip the sha256 check against MANIFEST.tsv.
#   -h | --help
# =============================================================================
set -euo pipefail

INPUT=""; OUTDIR=""; VERIFY=1
die() { echo "ERROR: $*" >&2; exit 1; }

for arg in "$@"; do
  case "$arg" in
    --input=*)   INPUT="${arg#*=}" ;;
    --out-dir=*) OUTDIR="${arg#*=}" ;;
    --no-verify) VERIFY=0 ;;
    -h|--help)   sed -n '2,28p' "$0"; exit 0 ;;
    *) die "unknown option: $arg (see --help)" ;;
  esac
done

[[ -n "$INPUT" ]] || die "--input is required (see --help)."
[[ -f "$INPUT" ]] || die "input not found: $INPUT"
INPUT="$(cd "$(dirname "$INPUT")" && pwd)/$(basename "$INPUT")"

# ---- Detect format -----------------------------------------------------------
mode=""
case "$INPUT" in
  *.tar.gz|*.tgz) mode=tar ;;
  *.txt)          mode=text ;;
  *) if gzip -t "$INPUT" 2>/dev/null; then mode=tar; else mode=text; fi ;;
esac

# ---- Default output directory ------------------------------------------------
if [[ -z "$OUTDIR" ]]; then
  stem="${INPUT%.tar.gz}"; stem="${stem%.tgz}"; stem="${stem%.txt}"
  OUTDIR="${stem}_extracted"
fi
mkdir -p "$OUTDIR"

echo "Unpacking ($mode): $INPUT"
echo "         -> $OUTDIR"

# ---- Extract -----------------------------------------------------------------
if [[ "$mode" == "tar" ]]; then
  # --strip-components=1 drops the leading consolidated_prs_results/ so the layout
  # matches the text-bundle case (<out-dir>/<run>/<file>, MANIFEST.tsv at the root).
  tar -xzf "$INPUT" -C "$OUTDIR" --strip-components=1
else
  awk -v outdir="$OUTDIR" '
    function dirof(p,   i){ for(i=length(p); i>0; i--) if(substr(p,i,1)=="/") return substr(p,1,i-1); return "." }
    /^===== BEGIN FILE: .* =====$/ {
      path=$0; sub(/^===== BEGIN FILE: /,"",path); sub(/ =====$/,"",path)
      target=outdir "/" path
      system("mkdir -p \"" outdir "/" dirof(path) "\"")
      writing=1; next
    }
    /^===== END FILE: .* =====$/ { if(writing){ close(target); writing=0 } next }
    writing { print > target }
  ' "$INPUT"
fi

n_extracted=$(find "$OUTDIR" -type f \( -name '*.csv' -o -name '*.tsv' -o -name '*.png' -o -name '*.pdf' \) | wc -l | tr -d ' ')
echo "Restored $n_extracted result file(s)."

# ---- Verify against MANIFEST.tsv --------------------------------------------
MAN="$OUTDIR/MANIFEST.tsv"
if [[ "$VERIFY" == "1" ]]; then
  [[ -f "$MAN" ]] || die "MANIFEST.tsv not found in $OUTDIR — cannot verify. Re-run with --no-verify to skip."
  ok=0; bad=0; missing=0
  # skip header (run\tfile\tsha256\t...)
  while IFS=$'\t' read -r run file sha bytes rows; do
    [[ "$run" == "run" || -z "$run" ]] && continue
    p="$OUTDIR/$run/$file"
    if [[ ! -f "$p" ]]; then echo "  ✗ missing: $run/$file"; missing=$((missing+1)); continue; fi
    got="$(shasum -a 256 "$p" | awk '{print $1}')"
    if [[ "$got" == "$sha" ]]; then ok=$((ok+1)); else echo "  ✗ checksum mismatch: $run/$file"; bad=$((bad+1)); fi
  done < "$MAN"
  echo ""
  echo "Integrity check vs MANIFEST.tsv: $ok OK, $bad mismatched, $missing missing."
  [[ "$bad" -eq 0 && "$missing" -eq 0 ]] || die "extract verification FAILED — do not use these files."
  echo "✅ All $ok table(s) verified. Restored under: $OUTDIR"
else
  echo "(integrity check skipped; --no-verify)"
  echo "Restored under: $OUTDIR"
fi
