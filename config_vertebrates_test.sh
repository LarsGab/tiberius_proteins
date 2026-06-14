#!/usr/bin/env bash
# Configuration for the vertebrates-test pipeline (order-level exclusion only).
# 6 species from tiberius_orf_finder/nextflow/conf/species_vertebrates_test.csv.
#
# Both genomes and Tiberius predictions are re-generated here (brain crash
# wiped the previous staging area), so this config drives prep steps 00 + 01
# in addition to the standard 02-07 miniprot analysis.
#
# Sources the main config first, then overrides species + adds new vars.

REPO_DIR="${REPO_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "$REPO_DIR/config.sh"

# ── Species list ──────────────────────────────────────────────────────────────
# All species use the Vertebrata ODB partition. Genomes are now staged by the
# tiberius_orf_finder pipeline at:
#   $VERTEBRATES_ASSEMBLY_DIR/<species>/assembly/genome.fa
VERTEBRATES_SPECIES_LIST=(
    "Archocentrus_centrarchus"
    "Betta_splendens"
    "Gallus_gallus"
    "Pristiophorus_japonicus"
    "Takifugu_rubripes"
    "Zootoca_vivipara"
    "Bos_taurus"
    "Delphinapterus_leucas"
    "Homo_sapiens"
)

# NCBI RefSeq accessions (from species_vertebrates_test.csv)
declare -A VERTEBRATES_ACCESSION=(
    [Archocentrus_centrarchus]="GCF_007364275.1"
    [Betta_splendens]="GCF_900634795.4"
    [Gallus_gallus]="GCF_016700215.2"
    [Pristiophorus_japonicus]="GCF_044704955.1"
    [Takifugu_rubripes]="GCF_901000725.2"
    [Zootoca_vivipara]="GCF_963506605.1"
    [Bos_taurus]="GCF_002263795.3"
    [Delphinapterus_leucas]="GCF_002288925.2"
    [Homo_sapiens]="GCF_000001405.40"
)

# Restructured genome layout (tiberius_orf_finder Nextflow output)
VERTEBRATES_ASSEMBLY_DIR="${VERTEBRATES_ASSEMBLY_DIR:-/home/gabriell/tiberius_orf_finder/results/vertebrates_test}"

# Add ODB partition + clade entries for each species
for SP in "${VERTEBRATES_SPECIES_LIST[@]}"; do
    ODB_PARTITION[$SP]="Vertebrata"
    SPECIES_CLADE[$SP]="Vertebrata"
done

# ── Order-level exclusion taxon IDs ───────────────────────────────────────────
# VERIFY EVERY ID before running: https://www.ncbi.nlm.nih.gov/taxonomy
# These are best-effort lookups; double-check that NCBI still uses the same
# order rank and that the taxid covers the intended subtree.
declare -A VERTEBRATES_ORDER_TAXON=(
    [Archocentrus_centrarchus]=1489872   # Cichliformes      VERIFY
    [Betta_splendens]=1489913            # Anabantiformes    VERIFY
    [Gallus_gallus]=8976                 # Galliformes       VERIFY
    [Pristiophorus_japonicus]=7779       # Pristiophoriformes (or Squaliformes) VERIFY
    [Takifugu_rubripes]=31022            # Tetraodontiformes VERIFY
    [Zootoca_vivipara]=8509              # Squamata          VERIFY
    [Bos_taurus]=91561                   # Artiodactyla (Cetartiodactyla) VERIFY
    [Delphinapterus_leucas]=91561        # Artiodactyla (Cetacea is infraorder under it) VERIFY
    [Homo_sapiens]=9443                  # Primates          VERIFY
)

# Populate EXCL_TAXON for the order level only (LEVEL 2)
for SP in "${VERTEBRATES_SPECIES_LIST[@]}"; do
    EXCL_TAXON[${SP}_2]="${VERTEBRATES_ORDER_TAXON[$SP]}"
done

# Convenience: only-order level loop
VERTEBRATES_LEVELS=(2)

# Number of (species, level) tasks for the array jobs
VERTEBRATES_N_TASKS=$(( ${#VERTEBRATES_SPECIES_LIST[@]} * ${#VERTEBRATES_LEVELS[@]} ))
VERTEBRATES_LAST_IDX=$(( VERTEBRATES_N_TASKS - 1 ))

# ── Tiberius prediction settings ──────────────────────────────────────────────
# Path to the Tiberius launcher on the cluster (cloned repo with tiberius.py).
TIBERIUS_LAUNCHER="${TIBERIUS_LAUNCHER:-$TIBERIUS_REPO/tiberius.py}"

# Pre-trained model config to use for these vertebrate species.
# See $TIBERIUS_REPO/model_cfg/README.md for available cfgs.
TIBERIUS_MODEL_CFG="${TIBERIUS_MODEL_CFG:-vertebrates}"

# Inference batch size (memory: brain prefers >= 200 for GPU throughput).
TIBERIUS_BATCH_SIZE="${TIBERIUS_BATCH_SIZE:-200}"
