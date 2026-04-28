#!/usr/bin/env bash
#SBATCH --job-name=parse_labels
#SBATCH --array=0-3
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/04_labels_%A_%a.log

# 4 tasks: one per species.
# Parses the gffcompare annotated.gtf for per-transcript class codes.
# If the annotated.gtf does not exist, reruns gffcompare first.

set -euo pipefail
source "${REPO_DIR}/config.sh"

ENTRIES=("${SPECIES_LIST[@]}")
ENTRY="${ENTRIES[$SLURM_ARRAY_TASK_ID]}"
SPECIES="${ENTRY%%:*}"
CLADE="${ENTRY##*:}"

GC_DIR="$BENCH_DIR/paper/bin/gffcompare_benchmark_output/$CLADE/$SPECIES/tiberius_seqlen"
ANNOTATED_GTF="$GC_DIR/tiberius_seqlen.annotated.gtf"
OUT="$WORK_DIR/labels/${SPECIES}_labels.tsv"

echo "[labels] Species=$SPECIES  Clade=$CLADE"

# Run gffcompare if the annotated GTF is missing
if [[ ! -f "$ANNOTATED_GTF" ]]; then
    echo "[labels] annotated.gtf not found – running gffcompare ..."
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
    echo "[labels] gffcompare complete."
fi

[[ -f "$ANNOTATED_GTF" ]] || { echo "ERROR: annotated.gtf still missing after gffcompare run"; exit 1; }

echo "[labels] Parsing $ANNOTATED_GTF ..."
python3 "$REPO_DIR/scripts/extract_labels.py" "$ANNOTATED_GTF" "$OUT"
