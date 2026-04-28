#!/usr/bin/env bash
#SBATCH --job-name=filter_odb
#SBATCH --array=0-11
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --time=03:00:00
#SBATCH --output=logs/01_filter_odb_%A_%a.log

# 12 tasks: 4 species × 3 exclusion levels
# Task index: species_idx * 3 + (level - 1)

set -euo pipefail
source "${REPO_DIR}/config.sh"

# Build ordered (SPECIES, LEVEL) list matching the array index
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

ODB="${ODB_PARTITION[$SPECIES]}"
EXCL_TAXID="${EXCL_TAXON[${SPECIES}_${LEVEL}]}"

INPUT="$WORK_DIR/odb/raw/${ODB}.fa.gz"
OUTPUT="$WORK_DIR/odb/filtered/${SPECIES}_excl_${LEVEL_LABEL}.fa"
NODES="$WORK_DIR/odb/nodes.dmp"

echo "[filter_odb] Species=$SPECIES  Level=$LEVEL_LABEL  ExclTaxID=$EXCL_TAXID"
echo "[filter_odb] $INPUT  ->  $OUTPUT"

python3 "$REPO_DIR/scripts/filter_odb.py" \
    "$INPUT" "$NODES" "$EXCL_TAXID" "$OUTPUT"

echo "[filter_odb] Done: $(wc -l < "$OUTPUT") lines written."
