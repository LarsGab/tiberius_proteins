#!/usr/bin/env bash
# Submit all pipeline jobs to SLURM with proper dependencies.
# Run this script from the repo root on the cluster:
#   bash run_all.sh

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export REPO_DIR

source "$REPO_DIR/config.sh"

mkdir -p "$WORK_DIR/logs"

# Helper: submit and capture job ID
submit() {
    sbatch --parsable "$@"
}

echo "=== Step 0: Download ODB partitions + NCBI taxonomy ==="
JID0=$(submit \
    --partition=pinky,snowball \
    --output="$WORK_DIR/logs/00_setup_%j.log" \
    "$REPO_DIR/slurm/00_setup.sh")
echo "  Job ID: $JID0"

echo "=== Step 1: Filter ODB at 3 exclusion levels (array 0-11) ==="
JID1=$(submit \
    --partition=pinky,snowball \
    --dependency=afterok:"$JID0" \
    --output="$WORK_DIR/logs/01_filter_odb_%A_%a.log" \
    "$REPO_DIR/slurm/01_filter_odb.sh")
echo "  Job ID: $JID1"

echo "=== Step 2: Extract Tiberius peptides via gffread (array 0-3) ==="
JID2=$(submit \
    --partition=pinky,snowball \
    --output="$WORK_DIR/logs/02_extract_pep_%A_%a.log" \
    "$REPO_DIR/slurm/02_extract_peptides.sh")
echo "  Job ID: $JID2"

echo "=== Step 3: Diamond makedb + blastp (array 0-11) ==="
JID3=$(submit \
    --partition=pinky,snowball \
    --dependency=afterok:"$JID1",afterok:"$JID2" \
    --output="$WORK_DIR/logs/03_diamond_%A_%a.log" \
    "$REPO_DIR/slurm/03_diamond.sh")
echo "  Job ID: $JID3"

echo "=== Step 4: Parse gffcompare labels (array 0-3) ==="
JID4=$(submit \
    --partition=pinky,snowball \
    --output="$WORK_DIR/logs/04_labels_%A_%a.log" \
    "$REPO_DIR/slurm/04_parse_labels.sh")
echo "  Job ID: $JID4"

echo "=== Step 5: Sweep filter rules → PR-curve TSVs (array 0-11) ==="
JID5=$(submit \
    --partition=pinky,snowball \
    --dependency=afterok:"$JID3",afterok:"$JID4" \
    --output="$WORK_DIR/logs/05_analyze_%A_%a.log" \
    "$REPO_DIR/slurm/05_analyze.sh")
echo "  Job ID: $JID5"

echo ""
echo "All jobs submitted."
echo "Monitor with:  squeue -u \$USER"
echo "Results will appear in:  $WORK_DIR/analysis/"
echo ""
echo "When done, download results:"
echo "  rsync -av brain:$WORK_DIR/analysis/ ./results/analysis/"
