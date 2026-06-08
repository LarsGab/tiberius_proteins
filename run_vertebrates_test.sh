#!/usr/bin/env bash
# Submit the vertebrates-test pipeline to SLURM (6 species × order level only).
#
# Pipeline:
#   00 download  -> NCBI datasets fetch genome + annot_cds.gff
#   01 tiberius  -> ab initio Tiberius predictions  (GPU)
#   02 filter_odb -> Vertebrata ODB minus the species' order subtree
#   03 extract_pep -> Tiberius CDS -> protein FASTA  (for Diamond pre-filter)
#   04 parse_labels -> gffcompare Tiberius vs RefSeq -> per-tx class codes
#   05 miniprot  -> Diamond pre-filter + miniprot + boundary scorer + miniprothint
#   06 gffcmp    -> Tiberius vs miniprothint -> per-tx support class codes
#   07 analyze   -> sweep filter rules -> PR-curve TSVs
#
# A vertebrates-specific setup job runs first to populate
#   $WORK_DIR/odb/raw/Vertebrata.fa.gz  and  $WORK_DIR/odb/nodes.dmp.
# It is idempotent: if both files already exist on disk it exits immediately.
#
# Run from the repo root on the cluster:
#   bash run_vertebrates_test.sh

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export REPO_DIR

source "$REPO_DIR/config_vertebrates_test.sh"

mkdir -p "$WORK_DIR/logs"

submit() {
    sbatch --parsable "$@"
}

# Default partition for CPU jobs. The Tiberius step overrides this with its own
# --partition=vision in its #SBATCH header.
PART="--partition=pinky,snowball"

echo "=== Setup (vertebrates_test): Download Vertebrata ODB + nodes.dmp ==="
JIDS=$(submit $PART \
    --output="$WORK_DIR/logs/vt_setup_%j.log" \
    "$REPO_DIR/slurm/vertebrates_test/setup.sh")
echo "  Job ID: $JIDS"

echo "=== Step 00 (vertebrates_test): Download RefSeq genomes ==="
JID0=$(submit $PART \
    --output="$WORK_DIR/logs/vt_00_download_%A_%a.log" \
    "$REPO_DIR/slurm/vertebrates_test/00_download_genomes.sh")
echo "  Job ID: $JID0"

echo "=== Step 01 (vertebrates_test): Tiberius ab initio predictions (GPU) ==="
JID1=$(submit \
    --dependency=afterok:"$JID0" \
    --output="$WORK_DIR/logs/vt_01_tiberius_%A_%a.log" \
    "$REPO_DIR/slurm/vertebrates_test/01_tiberius_predict.sh")
echo "  Job ID: $JID1"

echo "=== Step 02 (vertebrates_test): Filter ODB Vertebrata at order level ==="
JID2=$(submit $PART \
    --dependency=afterok:"$JIDS" \
    --output="$WORK_DIR/logs/vt_02_filter_odb_%A_%a.log" \
    "$REPO_DIR/slurm/vertebrates_test/02_filter_odb.sh")
echo "  Job ID: $JID2"

echo "=== Step 03 (vertebrates_test): Extract Tiberius peptides ==="
JID3=$(submit $PART \
    --dependency=afterok:"$JID1" \
    --output="$WORK_DIR/logs/vt_03_extract_pep_%A_%a.log" \
    "$REPO_DIR/slurm/vertebrates_test/03_extract_peptides.sh")
echo "  Job ID: $JID3"

echo "=== Step 04 (vertebrates_test): Parse gffcompare labels ==="
JID4=$(submit $PART \
    --dependency=afterok:"$JID1" \
    --output="$WORK_DIR/logs/vt_04_labels_%A_%a.log" \
    "$REPO_DIR/slurm/vertebrates_test/04_parse_labels.sh")
echo "  Job ID: $JID4"

echo "=== Step 05 (vertebrates_test): miniprot + boundary scorer + miniprothint ==="
JID5=$(submit $PART \
    --dependency=afterok:"$JID2",afterok:"$JID3" \
    --output="$WORK_DIR/logs/vt_05_miniprot_%A_%a.log" \
    "$REPO_DIR/slurm/vertebrates_test/05_miniprot.sh")
echo "  Job ID: $JID5"

echo "=== Step 06 (vertebrates_test): gffcompare Tiberius vs miniprothint ==="
JID6=$(submit $PART \
    --dependency=afterok:"$JID5",afterok:"$JID4" \
    --output="$WORK_DIR/logs/vt_06_miniprot_gffcmp_%A_%a.log" \
    "$REPO_DIR/slurm/vertebrates_test/06_miniprot_gffcmp.sh")
echo "  Job ID: $JID6"

echo "=== Step 07 (vertebrates_test): Analyze miniprothint-based filtering ==="
JID7=$(submit $PART \
    --dependency=afterok:"$JID6" \
    --output="$WORK_DIR/logs/vt_07_analyze_miniprot_%A_%a.log" \
    "$REPO_DIR/slurm/vertebrates_test/07_analyze_miniprot.sh")
echo "  Job ID: $JID7"

echo ""
echo "All vertebrates-test jobs submitted."
echo "Monitor with:  squeue -u \$USER"
echo ""
echo "When done, download:"
echo "  rsync -av brain:$WORK_DIR/miniprot_analysis/ ./results/miniprot_analysis/"
