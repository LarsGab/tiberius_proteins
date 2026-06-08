#!/usr/bin/env bash
#SBATCH --job-name=ins_analyze_miniprot
#SBATCH --array=0-9
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/ins_08_analyze_miniprot_%A_%a.log

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

ANNOT_LABELS="$WORK_DIR/labels/${SPECIES}_labels.tsv"
MP_LABELS_ALL="$WORK_DIR/miniprot_labels/${SPECIES}_excl_${LEVEL_LABEL}_all_labels.tsv"
MP_LABELS_TRAINING="$WORK_DIR/miniprot_labels/${SPECIES}_excl_${LEVEL_LABEL}_training_labels.tsv"

echo "[ins_analyze_miniprot] Species=$SPECIES  Level=$LEVEL_LABEL"

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

echo "[ins_analyze_miniprot] Done."
