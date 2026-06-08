#!/usr/bin/env bash
#SBATCH --job-name=ins_parse_labels
#SBATCH --array=0-9
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/ins_04_labels_%A_%a.log

# 10 tasks: one per insect species.
set -euo pipefail
source "${REPO_DIR}/config_insects.sh"

SPECIES="${INSECTS_SPECIES_LIST[$SLURM_ARRAY_TASK_ID]}"
CLADE="${SPECIES_CLADE[$SPECIES]}"

GC_DIR="$BENCH_DIR/paper/bin/gffcompare_benchmark_output/$CLADE/$SPECIES/tiberius_seqlen"
ANNOTATED_GTF="$GC_DIR/tiberius_seqlen.annotated.gtf"
OUT="$WORK_DIR/labels/${SPECIES}_labels.tsv"

echo "[ins_labels] Species=$SPECIES  Clade=$CLADE"

if [[ ! -f "$ANNOTATED_GTF" ]]; then
    echo "[ins_labels] annotated.gtf not found – running gffcompare ..."
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
    echo "[ins_labels] gffcompare complete."
fi

[[ -f "$ANNOTATED_GTF" ]] || { echo "ERROR: annotated.gtf still missing"; exit 1; }

mkdir -p "$(dirname "$OUT")"
echo "[ins_labels] Parsing $ANNOTATED_GTF ..."
python3 "$REPO_DIR/scripts/extract_labels.py" "$ANNOTATED_GTF" "$OUT"
