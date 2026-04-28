#!/usr/bin/env bash
# Central configuration for the protein-filtering analysis pipeline.
# Source this file from all SLURM scripts: source "$(dirname "$0")/../config.sh"

# ── Paths ─────────────────────────────────────────────────────────────────────
# Root of this repo on the cluster (set by run_all.sh via REPO_DIR env var,
# or auto-detected from this file's location)
REPO_DIR="${REPO_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

# Output directory for all generated data on the cluster
WORK_DIR=/home/gabriell/tiberius_proteins_analysis

# Tiberius benchmarking data root
BENCH_DIR=/home/gabriell/tiberius_benchmarking

# ── Species list: "SPECIES:CLADE" ─────────────────────────────────────────────
SPECIES_LIST=(
    "Drosophila_melanogaster:Insecta"
    "Homo_sapiens:Mammalia"
    "Arabidopsis_thaliana:Embryophyta"
    "Phaeodactylum_tricornutum:Diatoms"
)

# ── OrthoDB partition per species ─────────────────────────────────────────────
declare -A ODB_PARTITION=(
    [Drosophila_melanogaster]="Arthropoda"
    [Homo_sapiens]="Vertebrata"
    [Arabidopsis_thaliana]="Viridiplantae"
    [Phaeodactylum_tricornutum]="Stramenopiles"
)

ODB_BASE_URL="https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12"
NCBI_TAXONOMY_URL="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"

# ── Exclusion levels ──────────────────────────────────────────────────────────
declare -A LEVEL_NAMES=(
    [1]="genus"
    [2]="order"
    [3]="subphylum"
)

# ── Exclusion NCBI taxon IDs (inclusive subtree is removed from the DB) ───────
# VERIFY these before running: https://www.ncbi.nlm.nih.gov/taxonomy
#
# Drosophila melanogaster (taxid 7227)
#   L1 genus      : Drosophila      7215
#   L2 order      : Diptera         7147
#   L3 subphylum  : Hexapoda        6960
#
# Homo sapiens (taxid 9606)
#   L1 genus      : Homo            9605
#   L2 order      : Primates        9443
#   L3 class*     : Mammalia       40674  (*subphylum Vertebrata = whole partition)
#
# Arabidopsis thaliana (taxid 3702)
#   L1 genus      : Arabidopsis     3701
#   L2 order      : Brassicales     3699
#   L3 subphylum  : Embryophyta     3193
#
# Phaeodactylum tricornutum (taxid 2850)
#   L1 genus      : Phaeodactylum  33849
#   L2 order      : Naviculales    33653
#   L3 class      : Bacillariophyta 33634
declare -A EXCL_TAXON=(
    [Drosophila_melanogaster_1]=7215
    [Drosophila_melanogaster_2]=7147
    [Drosophila_melanogaster_3]=6960
    [Homo_sapiens_1]=9605
    [Homo_sapiens_2]=9443
    [Homo_sapiens_3]=40674
    [Arabidopsis_thaliana_1]=3701
    [Arabidopsis_thaliana_2]=3699
    [Arabidopsis_thaliana_3]=3193
    [Phaeodactylum_tricornutum_1]=33849
    [Phaeodactylum_tricornutum_2]=33653
    [Phaeodactylum_tricornutum_3]=33634
)

# ── SLURM defaults (override per-job as needed) ───────────────────────────────
SLURM_ACCOUNT=""          # leave empty to use cluster default
SLURM_PARTITION=""        # leave empty to use cluster default
DIAMOND_THREADS=16
