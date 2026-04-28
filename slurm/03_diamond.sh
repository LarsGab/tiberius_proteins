#!/usr/bin/env bash
#SBATCH --job-name=diamond
#SBATCH --array=0-11
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/03_diamond_%A_%a.log

# 12 tasks: 4 species × 3 exclusion levels.
# Each task: (1) builds a Diamond DB from the filtered ODB FASTA,
#            (2) runs blastp of the species' Tiberius peptides against it.
# Memory budget: Vertebrata.fa can be several GB; 64 GB should be safe.

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

QUERY="$WORK_DIR/peptides/${SPECIES}.pep.fa"
DB_FASTA="$WORK_DIR/odb/filtered/${SPECIES}_excl_${LEVEL_LABEL}.fa"
DB_PREFIX="$WORK_DIR/diamond/db/${SPECIES}_excl_${LEVEL_LABEL}"
OUT="$WORK_DIR/diamond/results/${SPECIES}_excl_${LEVEL_LABEL}.tsv"

echo "[diamond] Species=$SPECIES  Level=$LEVEL_LABEL"
echo "[diamond] Query    : $QUERY"
echo "[diamond] DB FASTA : $DB_FASTA"
echo "[diamond] Output   : $OUT"

for f in "$QUERY" "$DB_FASTA"; do
    [[ -f "$f" ]] || { echo "ERROR: file not found: $f"; exit 1; }
done

# Build Diamond database
echo "[diamond] Building database ..."
diamond makedb \
    --in "$DB_FASTA" \
    --db "$DB_PREFIX" \
    --threads "$DIAMOND_THREADS"

# Run blastp
# Output format includes qlen (col 13) and slen (col 14) for coverage calculation
echo "[diamond] Running blastp ..."
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
echo "[diamond] Done: $N hit lines written."
