#!/usr/bin/env python3
"""Filter an OrthoDB partitioned FASTA to remove all sequences from a given
NCBI taxon and all of its descendants.

OrthoDB12 sequence headers have the form:
    >TAXONID_N:GENEID ...
where TAXONID is the NCBI taxon ID of the source species (e.g. 7227_0:001234).

Usage:
    filter_odb.py INPUT_FASTA NODES_DMP EXCLUDE_TAXID OUTPUT_FASTA

Arguments:
    INPUT_FASTA    OrthoDB FASTA (plain or .gz)
    NODES_DMP      NCBI nodes.dmp from taxdump.tar.gz
    EXCLUDE_TAXID  NCBI taxon ID whose entire subtree will be excluded
    OUTPUT_FASTA   Output FASTA (plain text)
"""

from __future__ import annotations

import argparse
import gzip
import sys
from collections import defaultdict
from pathlib import Path


def load_nodes(nodes_dmp: Path) -> dict[int, int]:
    """Parse NCBI nodes.dmp, return {child_taxid: parent_taxid}."""
    parents: dict[int, int] = {}
    with open(nodes_dmp) as fh:
        for line in fh:
            parts = line.split("\t|\t")
            if len(parts) < 2:
                continue
            try:
                child = int(parts[0].strip())
                parent = int(parts[1].strip())
                parents[child] = parent
            except ValueError:
                continue
    return parents


def descendants(parents: dict[int, int], root: int) -> set[int]:
    """BFS – collect all taxon IDs rooted at *root* (inclusive)."""
    children: dict[int, list[int]] = defaultdict(list)
    for child, parent in parents.items():
        if child != parent:  # NCBI root node is its own parent
            children[parent].append(child)

    result: set[int] = set()
    queue: list[int] = [root]
    while queue:
        node = queue.pop()
        result.add(node)
        queue.extend(children.get(node, []))
    return result


def taxon_from_header(header: str) -> int | None:
    """Extract NCBI taxon ID from header '>7227_0:001234 ...' -> 7227."""
    name = header[1:].split()[0]   # strip '>' and trailing annotation
    taxon_str = name.split("_")[0]
    try:
        return int(taxon_str)
    except ValueError:
        return None


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("input",        help="OrthoDB FASTA (plain or .gz)")
    ap.add_argument("nodes_dmp",    help="NCBI nodes.dmp")
    ap.add_argument("exclude_taxid", type=int,
                    help="NCBI taxon ID whose subtree is excluded")
    ap.add_argument("output",       help="Output FASTA (plain text)")
    args = ap.parse_args(argv)

    print(f"Loading taxonomy from {args.nodes_dmp} ...", flush=True)
    parents = load_nodes(Path(args.nodes_dmp))

    print(f"Computing descendants of taxon {args.exclude_taxid} ...", flush=True)
    excl = descendants(parents, args.exclude_taxid)
    print(f"  -> {len(excl):,} taxa in exclusion set", flush=True)

    opener = gzip.open if str(args.input).endswith(".gz") else open
    n_total = n_excluded = 0
    keep = True

    with opener(args.input, "rt") as fin, open(args.output, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                n_total += 1
                taxon = taxon_from_header(line)
                keep = taxon is not None and taxon not in excl
                if not keep:
                    n_excluded += 1
            if keep:
                fout.write(line)

    kept = n_total - n_excluded
    print(
        f"Total: {n_total:,}  |  Excluded: {n_excluded:,}  |  Kept: {kept:,}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
