#!/usr/bin/env bash
#SBATCH --job-name=extract_pep
#SBATCH --array=0-3
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/02_extract_pep_%A_%a.log

# 4 tasks: one per species.
# Uses gffread to translate Tiberius CDS predictions to protein sequences.

set -euo pipefail
source "${REPO_DIR}/config.sh"

ENTRIES=("${SPECIES_LIST[@]}")
ENTRY="${ENTRIES[$SLURM_ARRAY_TASK_ID]}"
SPECIES="${ENTRY%%:*}"
CLADE="${ENTRY##*:}"

GENOME="$BENCH_DIR/$CLADE/$SPECIES/genome.fa"
GTF="$BENCH_DIR/paper/$CLADE/$SPECIES/results/predictions/tiberius/tiberius_seqlen.gtf"
OUT="$WORK_DIR/peptides/${SPECIES}.pep.fa"

echo "[extract_pep] Species=$SPECIES  Clade=$CLADE"
echo "[extract_pep] Genome  : $GENOME"
echo "[extract_pep] GTF     : $GTF"
echo "[extract_pep] Output  : $OUT"

for f in "$GENOME" "$GTF"; do
    [[ -f "$f" ]] || { echo "ERROR: file not found: $f"; exit 1; }
done

# -y  translate CDS to protein
# Strip terminal stop-codon asterisks that some tools add (Diamond may reject them)
gffread -y /dev/stdout -g "$GENOME" "$GTF" \
    | sed 's/\*//g' \
    > "$OUT"

N=$(grep -c '^>' "$OUT" || true)
echo "[extract_pep] Done: $N protein sequences written."
