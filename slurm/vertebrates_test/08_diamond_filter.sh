#!/usr/bin/env bash
#SBATCH --job-name=vt_diamond_filter
#SBATCH --array=0-8
#SBATCH --partition=pinky,snowball
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=200
#SBATCH --mem=150gb
#SBATCH --time=72:00:00
#SBATCH --output=logs/vt_08_diamond_filter_%A_%a.log

# One task per vertebrate test species: extract peptides from Tiberius GTF,
# run Diamond blastp against the order-excluded ODB, and keep genes whose
# proteins have a hit with qcov>=75 and tcov>=75.
#
# Inputs (per species):
#   GENOME    $VERTEBRATES_ASSEMBLY_DIR/<species>/assembly/genome.fa
#   GTF       $BENCH_DIR/paper/Vertebrata/<species>/results/predictions/tiberius/tiberius_seqlen.gtf
#   PROTEINS  $WORK_DIR/odb/filtered/<species>_excl_order.fa   (produced by 02_filter_odb.sh)
#
# Output:
#   $WORK_DIR/diamond_filter/<species>/   (filtered.gtf, filter_report.json,
#                                          peptides.fa, diamond.tsv)

set -euo pipefail
source "${REPO_DIR}/config_vertebrates_test.sh"

SPECIES="${VERTEBRATES_SPECIES_LIST[$SLURM_ARRAY_TASK_ID]}"

GENOME="$VERTEBRATES_ASSEMBLY_DIR/$SPECIES/assembly/genome.fa"
GTF="$BENCH_DIR/paper/Vertebrata/$SPECIES/results/predictions/tiberius/tiberius_seqlen.gtf"
PROTEINS="$WORK_DIR/odb/filtered/${SPECIES}_excl_order.fa"
OUT_DIR="$WORK_DIR/diamond_filter/$SPECIES"

echo "[vt_diamond_filter] Species  : $SPECIES"
echo "[vt_diamond_filter] Genome   : $GENOME"
echo "[vt_diamond_filter] GTF      : $GTF"
echo "[vt_diamond_filter] Proteins : $PROTEINS"
echo "[vt_diamond_filter] Out dir  : $OUT_DIR"

for f in "$GENOME" "$GTF" "$PROTEINS"; do
    [[ -f "$f" ]] || { echo "ERROR: file not found: $f"; exit 1; }
done

mkdir -p "$OUT_DIR"

python3 "$REPO_DIR/scripts/filter_tiberius_by_diamond.py" \
    --gtf       "$GTF" \
    --genome    "$GENOME" \
    --proteins  "$PROTEINS" \
    --out-dir   "$OUT_DIR" \
    --threads   "${SLURM_CPUS_PER_TASK:-200}" \
    --qcov-min  75 \
    --tcov-min  75

echo "[vt_diamond_filter] Done."
