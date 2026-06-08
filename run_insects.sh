#!/usr/bin/env bash
# Submit the extended-insects pipeline to SLURM (10 species × order level only).
# Steps 0 (download ODB+taxonomy) is NOT re-run – Arthropoda + nodes.dmp from
# the original pipeline are reused.
#
# Run from the repo root on the cluster:
#   bash run_insects.sh

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export REPO_DIR

source "$REPO_DIR/config_insects.sh"

mkdir -p "$WORK_DIR/logs"

submit() {
    sbatch --parsable "$@"
}

PART="--partition=pinky,snowball"

echo "=== Step 1 (insects): Filter ODB Arthropoda at order level ==="
JID1=$(submit $PART \
    --output="$WORK_DIR/logs/ins_01_filter_odb_%A_%a.log" \
    "$REPO_DIR/slurm/insects/01_filter_odb.sh")
echo "  Job ID: $JID1"

echo "=== Step 2 (insects): Extract Tiberius peptides ==="
JID2=$(submit $PART \
    --output="$WORK_DIR/logs/ins_02_extract_pep_%A_%a.log" \
    "$REPO_DIR/slurm/insects/02_extract_peptides.sh")
echo "  Job ID: $JID2"

echo "=== Step 3 (insects): Diamond makedb + blastp ==="
JID3=$(submit $PART \
    --dependency=afterok:"$JID1",afterok:"$JID2" \
    --output="$WORK_DIR/logs/ins_03_diamond_%A_%a.log" \
    "$REPO_DIR/slurm/insects/03_diamond.sh")
echo "  Job ID: $JID3"

echo "=== Step 4 (insects): Parse gffcompare labels ==="
JID4=$(submit $PART \
    --output="$WORK_DIR/logs/ins_04_labels_%A_%a.log" \
    "$REPO_DIR/slurm/insects/04_parse_labels.sh")
echo "  Job ID: $JID4"

echo "=== Step 5 (insects): Sweep filter rules → PR-curve TSVs ==="
JID5=$(submit $PART \
    --dependency=afterok:"$JID3",afterok:"$JID4" \
    --output="$WORK_DIR/logs/ins_05_analyze_%A_%a.log" \
    "$REPO_DIR/slurm/insects/05_analyze.sh")
echo "  Job ID: $JID5"

echo "=== Step 6 (insects): miniprot + boundary scorer + miniprothint ==="
JID6=$(submit $PART \
    --dependency=afterok:"$JID2",afterok:"$JID1" \
    --output="$WORK_DIR/logs/ins_06_miniprot_%A_%a.log" \
    "$REPO_DIR/slurm/insects/06_miniprot.sh")
echo "  Job ID: $JID6"

echo "=== Step 7 (insects): gffcompare Tiberius vs miniprothint ==="
JID7=$(submit $PART \
    --dependency=afterok:"$JID6",afterok:"$JID4" \
    --output="$WORK_DIR/logs/ins_07_miniprot_gffcmp_%A_%a.log" \
    "$REPO_DIR/slurm/insects/07_miniprot_gffcmp.sh")
echo "  Job ID: $JID7"

echo "=== Step 8 (insects): Analyze miniprothint-based filtering ==="
JID8=$(submit $PART \
    --dependency=afterok:"$JID7" \
    --output="$WORK_DIR/logs/ins_08_analyze_miniprot_%A_%a.log" \
    "$REPO_DIR/slurm/insects/08_analyze_miniprot.sh")
echo "  Job ID: $JID8"

echo ""
echo "All insect-extension jobs submitted."
echo "Monitor with:  squeue -u \$USER"
echo ""
echo "When done, download:"
echo "  rsync -av brain:$WORK_DIR/analysis/         ./results/analysis/"
echo "  rsync -av brain:$WORK_DIR/miniprot_analysis/ ./results/miniprot_analysis/"
