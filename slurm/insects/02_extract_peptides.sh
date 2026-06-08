#!/usr/bin/env bash
#SBATCH --job-name=ins_extract_pep
#SBATCH --array=0-9
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/ins_02_extract_pep_%A_%a.log

# 10 tasks: one per insect species. Translates Tiberius CDS to proteins.
set -euo pipefail
source "${REPO_DIR}/config_insects.sh"

SPECIES="${INSECTS_SPECIES_LIST[$SLURM_ARRAY_TASK_ID]}"
CLADE="${SPECIES_CLADE[$SPECIES]}"

GENOME="$BENCH_DIR/$CLADE/$SPECIES/genome.fa"
GTF="$BENCH_DIR/paper/$CLADE/$SPECIES/results/predictions/tiberius/tiberius_seqlen.gtf"
OUT="$WORK_DIR/peptides/${SPECIES}.pep.fa"

echo "[ins_extract_pep] Species=$SPECIES  Clade=$CLADE"
echo "[ins_extract_pep] Genome  : $GENOME"
echo "[ins_extract_pep] GTF     : $GTF"
echo "[ins_extract_pep] Output  : $OUT"

for f in "$GENOME" "$GTF"; do
    [[ -f "$f" ]] || { echo "ERROR: file not found: $f"; exit 1; }
done

mkdir -p "$(dirname "$OUT")"

# -y  translate CDS to protein. Strip stop-codon * and gffread mask '.'
gffread -y /dev/stdout -g "$GENOME" "$GTF" \
    | sed -e 's/\*//g' -e '/^>/!s/\.//g' \
    > "$OUT"

N=$(grep -c '^>' "$OUT" || true)
echo "[ins_extract_pep] Done: $N protein sequences written."
