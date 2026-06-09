#!/usr/bin/env python3
"""Locus-style Diamond filter for Tiberius predictions.

Replicates the miniprot+miniprothint filter (keep a Tiberius gene if its
locus has any protein hint) with Diamond hits instead: a gene is kept iff
at least one of its transcripts has a Diamond hit above a minimum e-value.
No pident / qcov tuning -- it is just "has any hit".

Two input modes:
  (A) provide a protein DB FASTA (and genome) -> the script runs gffread
      + diamond makedb + diamond blastp itself
  (B) provide a Diamond blastp TSV (outfmt 6, at least qseqid + evalue)
      -> Diamond step is skipped

The query IDs in the Diamond TSV must match the transcript_id values in
the GTF (gffread emits transcript_id as the FASTA header, which is what
the existing pipeline uses).

Usage (mode A):
    filter_tiberius_by_diamond.py \\
        --gtf       tiberius.gtf \\
        --genome    genome.fa \\
        --proteins  protein_db.fa \\
        --out-dir   results/diamond_filter/ \\
        [--threads 8] [--evalue 1e-5]

Usage (mode B):
    filter_tiberius_by_diamond.py \\
        --gtf         tiberius.gtf \\
        --diamond-tsv diamond.tsv \\
        --out-dir     results/diamond_filter/ \\
        [--evalue 1e-5]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# Reuse helpers from the existing transcript-level filter script
sys.path.insert(0, str(Path(__file__).resolve().parent))
from filter_genes_by_proteins import (  # noqa: E402
    _GENE_ID_RE,
    _TRANSCRIPT_ID_RE,
    extract_peptides,
    run_diamond,
    write_filtered_gtf,
)


def transcripts_with_hit(diamond_tsv: Path, evalue_max: float) -> set[str]:
    """Return the set of query (transcript) IDs with any hit at evalue <= max."""
    keep: set[str] = set()
    with open(diamond_tsv) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 11:
                continue
            qseqid = cols[0]
            try:
                evalue = float(cols[10])
            except ValueError:
                continue
            if evalue <= evalue_max:
                keep.add(qseqid)
    return keep


def gene_to_transcripts(gtf: Path) -> dict[str, set[str]]:
    """Map gene_id -> set of transcript_ids found in the GTF."""
    mapping: dict[str, set[str]] = {}
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            tid = _TRANSCRIPT_ID_RE.search(line)
            gid = _GENE_ID_RE.search(line)
            if tid and gid:
                mapping.setdefault(gid.group(1), set()).add(tid.group(1))
    return mapping


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--gtf", required=True, type=Path,
                    help="Tiberius predictions (GTF)")
    ap.add_argument("--out-dir", required=True, type=Path)

    # Mode A: run Diamond ourselves
    ap.add_argument("--genome", type=Path,
                    help="Genome FASTA (required if --diamond-tsv not given)")
    ap.add_argument("--proteins", type=Path,
                    help="Protein DB FASTA (required if --diamond-tsv not given)")
    ap.add_argument("--peptides-fa", type=Path,
                    help="Skip gffread, use this peptide FASTA as Diamond query")

    # Mode B: pre-computed Diamond TSV
    ap.add_argument("--diamond-tsv", type=Path,
                    help="Pre-computed Diamond blastp TSV (outfmt 6). "
                         "If given, --genome/--proteins are ignored.")

    # Filtering
    ap.add_argument("--evalue", type=float, default=1e-5,
                    help="Max e-value for a hit to count as support (default 1e-5)")

    # Diamond knobs (only used in mode A)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--sensitivity",
                    choices=["--fast", "--mid-sensitive", "--sensitive",
                             "--more-sensitive", "--very-sensitive",
                             "--ultra-sensitive"],
                    default="--more-sensitive")
    ap.add_argument("--max-target-seqs", type=int, default=1,
                    help="Only need one hit per query (default 1)")
    ap.add_argument("--force", action="store_true",
                    help="Recompute peptides / Diamond even if outputs exist")

    return ap.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    if args.diamond_tsv is None:
        if args.genome is None or args.proteins is None:
            sys.exit("ERROR: provide either --diamond-tsv OR --genome + --proteins")

        peptides_fa = args.peptides_fa or extract_peptides(
            args.genome, args.gtf, args.out_dir / "peptides.fa", force=args.force,
        )
        diamond_tsv = run_diamond(
            peptides_fa, args.proteins, args.out_dir,
            threads=args.threads,
            evalue=args.evalue,
            sensitivity=args.sensitivity,
            max_target_seqs=args.max_target_seqs,
            force=args.force,
        )
    else:
        diamond_tsv = args.diamond_tsv
        print(f"Using pre-computed Diamond TSV: {diamond_tsv}", flush=True)

    print(f"Loading hits at evalue <= {args.evalue:g}", flush=True)
    supported_tids = transcripts_with_hit(diamond_tsv, args.evalue)

    gene_tids = gene_to_transcripts(args.gtf)
    n_genes_total = len(gene_tids)
    n_tids_total  = sum(len(v) for v in gene_tids.values())

    # Locus-style: a gene is kept if ANY of its transcripts has a hit.
    # All transcripts of a kept gene are then retained.
    kept_genes = {gid for gid, tids in gene_tids.items()
                  if tids & supported_tids}
    kept_tids: set[str] = set()
    for gid in kept_genes:
        kept_tids.update(gene_tids[gid])

    print(f"  transcripts with hit : {len(supported_tids):,} / {n_tids_total:,}",
          flush=True)
    print(f"  genes kept           : {len(kept_genes):,} / {n_genes_total:,} "
          f"({len(kept_genes)/max(1,n_genes_total):.1%})", flush=True)
    print(f"  transcripts kept     : {len(kept_tids):,} / {n_tids_total:,}",
          flush=True)

    out_gtf = args.out_dir / "filtered.gtf"
    n_kept_lines, n_total_lines = write_filtered_gtf(args.gtf, out_gtf, kept_tids)
    print(f"Wrote {out_gtf}  ({n_kept_lines:,} / {n_total_lines:,} lines)",
          flush=True)

    report = {
        "inputs": {
            "gtf":         str(args.gtf),
            "genome":      str(args.genome)   if args.genome   else None,
            "proteins":    str(args.proteins) if args.proteins else None,
            "diamond_tsv": str(diamond_tsv),
        },
        "evalue_max": args.evalue,
        "result": {
            "n_genes_total":         n_genes_total,
            "n_genes_kept":          len(kept_genes),
            "n_transcripts_total":   n_tids_total,
            "n_transcripts_with_hit":len(supported_tids),
            "n_transcripts_kept":    len(kept_tids),
            "n_lines_total":         n_total_lines,
            "n_lines_kept":          n_kept_lines,
        },
    }
    (args.out_dir / "filter_report.json").write_text(
        json.dumps(report, indent=2, default=str)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
