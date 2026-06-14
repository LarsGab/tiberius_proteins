#!/usr/bin/env bash
#SBATCH --job-name=vt_filter_odb
#SBATCH --array=0-8
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --time=03:00:00
#SBATCH --output=logs/vt_02_filter_odb_%A_%a.log

# 9 tasks: one per vertebrate test species (order level only).
# Uses the Vertebrata ODB partition + nodes.dmp downloaded by the original
# slurm/00_setup.sh — re-run that first if those are missing.

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

ODB="${ODB_PARTITION[$SPECIES]}"
EXCL_TAXID="${EXCL_TAXON[${SPECIES}_${LEVEL}]}"

INPUT="$WORK_DIR/odb/raw/${ODB}.fa.gz"
OUTPUT="$WORK_DIR/odb/filtered/${SPECIES}_excl_${LEVEL_LABEL}.fa"
NODES="$WORK_DIR/odb/nodes.dmp"

echo "[vt_filter_odb] Species=$SPECIES  Level=$LEVEL_LABEL  ExclTaxID=$EXCL_TAXID"
echo "[vt_filter_odb] $INPUT  ->  $OUTPUT"

for f in "$INPUT" "$NODES"; do
    [[ -f "$f" ]] || { echo "ERROR: file not found: $f (run slurm/00_setup.sh)"; exit 1; }
done

mkdir -p "$(dirname "$OUTPUT")"

python3 "$REPO_DIR/scripts/filter_odb.py" \
    "$INPUT" "$NODES" "$EXCL_TAXID" "$OUTPUT"

echo "[vt_filter_odb] Done: $(wc -l < "$OUTPUT") lines written."
