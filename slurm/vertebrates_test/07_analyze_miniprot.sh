#!/usr/bin/env bash
#SBATCH --job-name=vt_analyze_miniprot
#SBATCH --array=0-5
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/vt_07_analyze_miniprot_%A_%a.log

# 6 tasks: sweep miniprothint-based filter rules and write PR-curve TSVs
# for each vertebrate test species at order-level exclusion.

set -euo pipefail
source "${REPO_DIR}/config_vertebrates_test.sh"

TASKS=()
for SPECIES in "${VERTEBRATES_SPECIES_LIST[@]}"; do
    for LEVEL in "${VERTEBRATES_LEVELS[@]}"; do
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

echo "[vt_analyze_miniprot] Species=$SPECIES  Level=$LEVEL_LABEL"

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

echo "[vt_analyze_miniprot] Done."
