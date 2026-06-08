#!/usr/bin/env bash
#SBATCH --job-name=vt_parse_labels
#SBATCH --array=0-5
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/vt_04_labels_%A_%a.log

# 6 tasks: one per vertebrate test species.
# Builds the ground-truth annotation labels TSV (Tiberius transcript class
# codes vs RefSeq annotation) needed by step 07.

set -euo pipefail
source "${REPO_DIR}/config_vertebrates_test.sh"

SPECIES="${VERTEBRATES_SPECIES_LIST[$SLURM_ARRAY_TASK_ID]}"
CLADE="${SPECIES_CLADE[$SPECIES]}"

GC_DIR="$BENCH_DIR/paper/bin/gffcompare_benchmark_output/$CLADE/$SPECIES/tiberius_seqlen"
ANNOTATED_GTF="$GC_DIR/tiberius_seqlen.annotated.gtf"
OUT="$WORK_DIR/labels/${SPECIES}_labels.tsv"

echo "[vt_labels] Species=$SPECIES  Clade=$CLADE"

if [[ ! -f "$ANNOTATED_GTF" ]]; then
    echo "[vt_labels] annotated.gtf not found – running gffcompare ..."
    ANNOT="$BENCH_DIR/$CLADE/$SPECIES/annot_cds.gff"
    GTF="$BENCH_DIR/paper/$CLADE/$SPECIES/results/predictions/tiberius/tiberius_seqlen.gtf"
    for f in "$ANNOT" "$GTF"; do
        [[ -f "$f" ]] || { echo "ERROR: file not found: $f"; exit 1; }
    done
    mkdir -p "$GC_DIR"
    gffcompare \
        --strict-match \
        -e 3 \
        -T \
        -r "$ANNOT" \
        "$GTF" \
        -o "$GC_DIR/tiberius_seqlen"
    echo "[vt_labels] gffcompare complete."
fi

[[ -f "$ANNOTATED_GTF" ]] || { echo "ERROR: annotated.gtf still missing"; exit 1; }

mkdir -p "$(dirname "$OUT")"
echo "[vt_labels] Parsing $ANNOTATED_GTF ..."
python3 "$REPO_DIR/scripts/extract_labels.py" "$ANNOTATED_GTF" "$OUT"
echo "[vt_labels] Done: $OUT"
