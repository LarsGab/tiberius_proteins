#!/usr/bin/env bash
#SBATCH --job-name=analyze
#SBATCH --array=0-11
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/05_analyze_%A_%a.log

# 12 tasks: 4 species × 3 exclusion levels.
# Sweeps filter rules over Diamond results + gffcompare labels,
# writes per-(species, level) precision/recall TSV tables to WORK_DIR/analysis/.

set -euo pipefail
source "${REPO_DIR}/config.sh"

TASKS=()
for entry in "${SPECIES_LIST[@]}"; do
    SPECIES="${entry%%:*}"
    for LEVEL in 1 2 3; do
        TASKS+=("${SPECIES}:${LEVEL}")
    done
done

TASK="${TASKS[$SLURM_ARRAY_TASK_ID]}"
SPECIES="${TASK%%:*}"
LEVEL="${TASK##*:}"
LEVEL_LABEL="${LEVEL_NAMES[$LEVEL]}"

DIAMOND_TSV="$WORK_DIR/diamond/results/${SPECIES}_excl_${LEVEL_LABEL}.tsv"
LABELS_TSV="$WORK_DIR/labels/${SPECIES}_labels.tsv"

echo "[analyze] Species=$SPECIES  Level=$LEVEL_LABEL"

for f in "$DIAMOND_TSV" "$LABELS_TSV"; do
    [[ -f "$f" ]] || { echo "ERROR: file not found: $f"; exit 1; }
done

python3 "$REPO_DIR/scripts/analyze_filter_rules.py" \
    --diamond  "$DIAMOND_TSV" \
    --labels   "$LABELS_TSV" \
    --species  "$SPECIES" \
    --level    "$LEVEL_LABEL" \
    --out-dir  "$WORK_DIR/analysis"

echo "[analyze] Done."
