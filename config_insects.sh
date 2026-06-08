#!/usr/bin/env bash
# Configuration for the extended-insects analysis (order-level exclusion only).
# Sources the main config first, then overrides the species-related arrays.

REPO_DIR="${REPO_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "$REPO_DIR/config.sh"

# ── Species list ──────────────────────────────────────────────────────────────
# All 10 species belong to clade Insecta and use the Arthropoda ODB partition.
# CLADE is set to Insecta to match $BENCH_DIR/Insecta/<species>/ layout.
INSECTS_SPECIES_LIST=(
    "Bombyx_mori"
    "Colias_croceus"
    "Nymphalis_io"
    "Tribolium_castaneum"
    "Zerene_cesonia"
    "Cataglyphis_hispanica"
    "Danaus_plexippus"
    "Leptidea_sinapis"
    "Osmia_bicornis"
    "Vanessa_cardui"
)

# Add ODB partition + clade entries for each new species
for SP in "${INSECTS_SPECIES_LIST[@]}"; do
    ODB_PARTITION[$SP]="Arthropoda"
    SPECIES_CLADE[$SP]="Insecta"
done

# ── Order-level exclusion taxon IDs ───────────────────────────────────────────
# VERIFY: https://www.ncbi.nlm.nih.gov/taxonomy
# Lepidoptera: 7088   Coleoptera: 7041   Hymenoptera: 7399
declare -A INSECTS_ORDER_TAXON=(
    [Bombyx_mori]=7088              # Lepidoptera
    [Colias_croceus]=7088           # Lepidoptera
    [Nymphalis_io]=7088             # Lepidoptera
    [Zerene_cesonia]=7088           # Lepidoptera
    [Danaus_plexippus]=7088         # Lepidoptera
    [Leptidea_sinapis]=7088         # Lepidoptera
    [Vanessa_cardui]=7088           # Lepidoptera
    [Tribolium_castaneum]=7041      # Coleoptera
    [Cataglyphis_hispanica]=7399    # Hymenoptera
    [Osmia_bicornis]=7399           # Hymenoptera
)

# Populate EXCL_TAXON for the order level only (LEVEL 2)
for SP in "${INSECTS_SPECIES_LIST[@]}"; do
    EXCL_TAXON[${SP}_2]="${INSECTS_ORDER_TAXON[$SP]}"
done

# Convenience: only-order level loop
INSECTS_LEVELS=(2)

# Number of (species, level) tasks for the array jobs
INSECTS_N_TASKS=$(( ${#INSECTS_SPECIES_LIST[@]} * ${#INSECTS_LEVELS[@]} ))
INSECTS_LAST_IDX=$(( INSECTS_N_TASKS - 1 ))
