#!/usr/bin/env bash
#SBATCH --job-name=ins_filter_odb
#SBATCH --array=0-9
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --time=03:00:00
#SBATCH --output=logs/ins_01_filter_odb_%A_%a.log

# 10 tasks: one per insect species (order level only).
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

ODB="${ODB_PARTITION[$SPECIES]}"
EXCL_TAXID="${EXCL_TAXON[${SPECIES}_${LEVEL}]}"

INPUT="$WORK_DIR/odb/raw/${ODB}.fa.gz"
OUTPUT="$WORK_DIR/odb/filtered/${SPECIES}_excl_${LEVEL_LABEL}.fa"
NODES="$WORK_DIR/odb/nodes.dmp"

echo "[ins_filter_odb] Species=$SPECIES  Level=$LEVEL_LABEL  ExclTaxID=$EXCL_TAXID"
echo "[ins_filter_odb] $INPUT  ->  $OUTPUT"

python3 "$REPO_DIR/scripts/filter_odb.py" \
    "$INPUT" "$NODES" "$EXCL_TAXID" "$OUTPUT"

echo "[ins_filter_odb] Done: $(wc -l < "$OUTPUT") lines written."
