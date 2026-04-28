#!/usr/bin/env bash
#SBATCH --job-name=analyze_miniprot
#SBATCH --array=0-11
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/08_analyze_miniprot_%A_%a.log

# 12 tasks: 4 species × 3 exclusion levels.
# Sweeps miniprothint-based filter rules and writes PR-curve TSVs.

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

ANNOT_LABELS="$WORK_DIR/labels/${SPECIES}_labels.tsv"
MP_LABELS_ALL="$WORK_DIR/miniprot_labels/${SPECIES}_excl_${LEVEL_LABEL}_all_labels.tsv"
MP_LABELS_TRAINING="$WORK_DIR/miniprot_labels/${SPECIES}_excl_${LEVEL_LABEL}_training_labels.tsv"

echo "[analyze_miniprot] Species=$SPECIES  Level=$LEVEL_LABEL"

for f in "$ANNOT_LABELS" "$MP_LABELS_ALL" "$MP_LABELS_TRAINING"; do
    [[ -f "$f" ]] || { echo "ERROR: file not found: $f"; exit 1; }
done

python3 "$REPO_DIR/scripts/analyze_miniprot_filter.py" \
    --annot-labels              "$ANNOT_LABELS" \
    --miniprot-labels-all       "$MP_LABELS_ALL" \
    --miniprot-labels-training  "$MP_LABELS_TRAINING" \
    --species  "$SPECIES" \
    --level    "$LEVEL_LABEL" \
    --out-dir  "$WORK_DIR/miniprot_analysis"

echo "[analyze_miniprot] Done."
