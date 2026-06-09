#!/usr/bin/env bash
#SBATCH --job-name=vt_extract_pep
#SBATCH --array=0-5
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/vt_03_extract_pep_%A_%a.log

# 6 tasks: one per vertebrate test species. Translates Tiberius CDS to proteins
# for the Diamond preprocessing step inside 05_miniprot.sh.

set -euo pipefail
source /etc/profile.d/modules.sh
module load singularity/3.11.3
source "${REPO_DIR}/config_vertebrates_test.sh"

SPECIES="${VERTEBRATES_SPECIES_LIST[$SLURM_ARRAY_TASK_ID]}"
CLADE="${SPECIES_CLADE[$SPECIES]}"

GENOME="$BENCH_DIR/$CLADE/$SPECIES/genome.fa"
GTF="$BENCH_DIR/paper/$CLADE/$SPECIES/results/predictions/tiberius/tiberius_seqlen.gtf"
OUT="$WORK_DIR/peptides/${SPECIES}.pep.fa"

echo "[vt_extract_pep] Species=$SPECIES  Clade=$CLADE"
echo "[vt_extract_pep] Genome  : $GENOME"
echo "[vt_extract_pep] GTF     : $GTF"
echo "[vt_extract_pep] Output  : $OUT"

for f in "$GENOME" "$GTF"; do
    [[ -f "$f" ]] || { echo "ERROR: file not found: $f"; exit 1; }
done

mkdir -p "$(dirname "$OUT")"

run_tool() {
    if [[ -n "${TIBERIUS_SIF:-}" ]]; then
        singularity exec "$TIBERIUS_SIF" "$@"
    else
        "$@"
    fi
}

# -y translate CDS to protein. Strip stop-codon * and gffread mask '.'
run_tool gffread -y /dev/stdout -g "$GENOME" "$GTF" \
    | sed -e 's/\*//g' -e '/^>/!s/\.//g' \
    > "$OUT"

N=$(grep -c '^>' "$OUT" || true)
echo "[vt_extract_pep] Done: $N protein sequences written."
