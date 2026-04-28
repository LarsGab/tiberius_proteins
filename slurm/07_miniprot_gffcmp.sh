#!/usr/bin/env bash
#SBATCH --job-name=miniprot_gffcmp
#SBATCH --array=0-11
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/07_miniprot_gffcmp_%A_%a.log

# 12 tasks: 4 species × 3 exclusion levels.
# For each task:
#   gffcompare -r miniprot.gtf tiberius.gtf  -> per-Tiberius-transcript class
#              codes relative to miniprothint gene models.
#   extract_labels.py on the resulting annotated.gtf -> TSV for analysis.
#
# Two comparisons are run per task:
#   (a) miniprot.gtf              (all miniprothint models)
#   (b) miniprot_trainingGenes.gff (high-confidence subset)

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
CLADE="${SPECIES_CLADE[$SPECIES]}"

TIBERIUS_GTF="$BENCH_DIR/paper/$CLADE/$SPECIES/results/predictions/tiberius/tiberius_seqlen.gtf"
MINIPROT_DIR="$WORK_DIR/miniprot/${SPECIES}_excl_${LEVEL_LABEL}"
OUT_DIR="$WORK_DIR/miniprot_labels"

echo "[gffcmp] Species=$SPECIES  Level=$LEVEL_LABEL"
echo "[gffcmp] Tiberius GTF : $TIBERIUS_GTF"
echo "[gffcmp] Miniprot dir : $MINIPROT_DIR"

for f in "$TIBERIUS_GTF" \
         "$MINIPROT_DIR/miniprot.gtf" \
         "$MINIPROT_DIR/miniprot_trainingGenes.gff"; do
    [[ -f "$f" ]] || { echo "ERROR: file not found: $f"; exit 1; }
done

mkdir -p "$OUT_DIR"

run_gffcmp_and_extract() {
    local ref="$1"
    local label="$2"   # "all" or "training"
    local gc_prefix="$MINIPROT_DIR/gffcmp_${label}"
    local out_tsv="$OUT_DIR/${SPECIES}_excl_${LEVEL_LABEL}_${label}_labels.tsv"

    echo "[gffcmp] Comparing Tiberius vs miniprothint ($label) ..."
    gffcompare \
        -e 3 \
        -T \
        -r "$ref" \
        "$TIBERIUS_GTF" \
        -o "$gc_prefix"

    local annotated="${gc_prefix}.$(basename "$TIBERIUS_GTF").tmap"
    # gffcompare names the annotated GTF after the query filename
    local annotated_gtf="${gc_prefix}.annotated.gtf"
    [[ -f "$annotated_gtf" ]] || {
        echo "ERROR: gffcompare did not produce annotated GTF at $annotated_gtf"
        exit 1
    }

    python3 "$REPO_DIR/scripts/extract_labels.py" \
        "$annotated_gtf" \
        "$out_tsv"

    echo "[gffcmp] Done ($label): $out_tsv"
}

run_gffcmp_and_extract "$MINIPROT_DIR/miniprot.gtf"             "all"
run_gffcmp_and_extract "$MINIPROT_DIR/miniprot_trainingGenes.gff" "training"

echo "[gffcmp] Step 07 complete."
