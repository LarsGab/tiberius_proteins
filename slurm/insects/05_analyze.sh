#!/usr/bin/env bash
#SBATCH --job-name=ins_analyze
#SBATCH --array=0-9
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/ins_05_analyze_%A_%a.log

# 10 tasks: 10 insects × 1 level (order).
set -euo pipefail
source "${REPO_DIR}/config_insects.sh"

TASKS=()
for SPECIES in "${INSECTS_SPECIES_LIST[@]}"; do
    for LEVEL in "${INSECTS_LEVELS[@]}"; do
        TASKS+=("${SPECIES}:${LEVEL}")
    done
done

TASK="${TASKS[$SLURM_ARRAY_TASK_ID]}"
SPECIES="${TASK%%:*}"
LEVEL="${TASK##*:}"
LEVEL_LABEL="${LEVEL_NAMES[$LEVEL]}"

DIAMOND_TSV="$WORK_DIR/diamond/results/${SPECIES}_excl_${LEVEL_LABEL}.tsv"
LABELS_TSV="$WORK_DIR/labels/${SPECIES}_labels.tsv"

echo "[ins_analyze] Species=$SPECIES  Level=$LEVEL_LABEL"

for f in "$DIAMOND_TSV" "$LABELS_TSV"; do
    [[ -f "$f" ]] || { echo "ERROR: file not found: $f"; exit 1; }
done

python3 "$REPO_DIR/scripts/analyze_filter_rules.py" \
    --diamond  "$DIAMOND_TSV" \
    --labels   "$LABELS_TSV" \
    --species  "$SPECIES" \
    --level    "$LEVEL_LABEL" \
    --out-dir  "$WORK_DIR/analysis"

echo "[ins_analyze] Done."
