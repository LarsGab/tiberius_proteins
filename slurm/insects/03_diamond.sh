#!/usr/bin/env bash
#SBATCH --job-name=ins_diamond
#SBATCH --array=0-9
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/ins_03_diamond_%A_%a.log

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

QUERY="$WORK_DIR/peptides/${SPECIES}.pep.fa"
DB_FASTA="$WORK_DIR/odb/filtered/${SPECIES}_excl_${LEVEL_LABEL}.fa"
DB_PREFIX="$WORK_DIR/diamond/db/${SPECIES}_excl_${LEVEL_LABEL}"
OUT="$WORK_DIR/diamond/results/${SPECIES}_excl_${LEVEL_LABEL}.tsv"

echo "[ins_diamond] Species=$SPECIES  Level=$LEVEL_LABEL"
for f in "$QUERY" "$DB_FASTA"; do
    [[ -f "$f" ]] || { echo "ERROR: file not found: $f"; exit 1; }
done

mkdir -p "$(dirname "$DB_PREFIX")" "$(dirname "$OUT")"

echo "[ins_diamond] Building database ..."
diamond makedb \
    --in "$DB_FASTA" \
    --db "$DB_PREFIX" \
    --threads "$DIAMOND_THREADS"

echo "[ins_diamond] Running blastp ..."
diamond blastp \
    --query "$QUERY" \
    --db "${DB_PREFIX}.dmnd" \
    --out "$OUT" \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen \
               qstart qend sstart send evalue bitscore qlen slen \
    --max-target-seqs 5 \
    --evalue 1e-5 \
    --more-sensitive \
    --threads "$DIAMOND_THREADS"

N=$(wc -l < "$OUT" || echo 0)
echo "[ins_diamond] Done: $N hit lines written."
