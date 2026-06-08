#!/usr/bin/env bash
#SBATCH --job-name=vt_tiberius
#SBATCH --array=0-5
#SBATCH --partition=vision
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --output=logs/vt_01_tiberius_%A_%a.log

# 6 tasks: one per vertebrate test species.
# Runs Tiberius ab initio gene prediction on each downloaded genome and
# writes the GTF to the layout that the rest of the pipeline expects:
#   $BENCH_DIR/paper/$CLADE/$SPECIES/results/predictions/tiberius/tiberius_seqlen.gtf
#
# Uses the Tiberius launcher (tiberius.py) inside Singularity.

set -euo pipefail
source /etc/profile.d/modules.sh
module load singularity/3.11.3
source "${REPO_DIR}/config_vertebrates_test.sh"

SPECIES="${VERTEBRATES_SPECIES_LIST[$SLURM_ARRAY_TASK_ID]}"
CLADE="${SPECIES_CLADE[$SPECIES]}"

GENOME="$BENCH_DIR/$CLADE/$SPECIES/genome.fa"
OUT_DIR="$BENCH_DIR/paper/$CLADE/$SPECIES/results/predictions/tiberius"
OUT_GTF="$OUT_DIR/tiberius_seqlen.gtf"

echo "[vt_tiberius] Species=$SPECIES  Clade=$CLADE"
echo "[vt_tiberius] Genome     : $GENOME"
echo "[vt_tiberius] Model cfg  : $TIBERIUS_MODEL_CFG"
echo "[vt_tiberius] Output GTF : $OUT_GTF"

[[ -f "$GENOME"          ]] || { echo "ERROR: genome not found: $GENOME"; exit 1; }
[[ -f "$TIBERIUS_LAUNCHER" ]] || { echo "ERROR: tiberius launcher not found: $TIBERIUS_LAUNCHER"; exit 1; }

if [[ -s "$OUT_GTF" ]]; then
    echo "[vt_tiberius] Output exists, skipping."
    exit 0
fi

mkdir -p "$OUT_DIR"

# Tiberius launcher invokes Singularity itself when given --singularity, so we
# do not wrap it in `singularity exec` here. The launcher needs Python with the
# tiberius package installed; activate the orffinder micromamba env (per memory).
eval "$(micromamba shell hook --shell bash)"
micromamba activate orffinder

python "$TIBERIUS_LAUNCHER" \
    --singularity \
    --genome     "$GENOME" \
    --model_cfg  "$TIBERIUS_MODEL_CFG" \
    --batch_size "$TIBERIUS_BATCH_SIZE" \
    --out        "$OUT_GTF"

[[ -s "$OUT_GTF" ]] || { echo "ERROR: Tiberius produced no GTF"; exit 2; }

N_TX=$(awk -F'\t' '$3=="transcript"' "$OUT_GTF" | wc -l || true)
echo "[vt_tiberius] Done. $N_TX transcripts predicted."
