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


_SCI_NOTATION_HINT = ("e", "E")


def _looks_like_evalue(tok: str) -> bool:
    """A token is e-value-like if it parses as float AND uses sci notation."""
    if not any(c in tok for c in _SCI_NOTATION_HINT):
        return False
    try:
        float(tok)
    except ValueError:
        return False
    return True


def detect_evalue_col(diamond_tsv: Path, scan_lines: int = 50) -> int:
    """Return the 0-based column index that looks like the e-value column.

    Heuristic: find the column whose tokens consistently use scientific
    notation (e.g. '1e-5', '3.33e-239'). Diamond writes bitscore as a plain
    number, so this uniquely identifies e-value across common outfmt 6
    variants (default 12-col -> idx 10; custom layouts -> different idx).
    """
    counts: dict[int, int] = {}
    n_lines = 0
    with open(diamond_tsv) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            for i, tok in enumerate(cols):
                if _looks_like_evalue(tok):
                    counts[i] = counts.get(i, 0) + 1
            n_lines += 1
            if n_lines >= scan_lines:
                break
    if not counts:
        raise ValueError(
            f"Could not auto-detect e-value column in {diamond_tsv}. "
            "Pass --evalue-col explicitly (1-based)."
        )
    # Column with the most sci-notation tokens wins; tiebreak on lowest index.
    best_idx = max(counts, key=lambda i: (counts[i], -i))
    return best_idx


def transcripts_with_hit(
    diamond_tsv: Path,
    evalue_max: float,
    evalue_col: int,
) -> set[str]:
    """Return the set of query (transcript) IDs with any hit at evalue <= max.

    `evalue_col` is the 0-based column index of the e-value.
    """
    keep: set[str] = set()
    with open(diamond_tsv) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) <= evalue_col:
                continue
            try:
                evalue = float(cols[evalue_col])
            except ValueError:
                continue
            if evalue <= evalue_max:
                keep.add(cols[0])
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
    ap.add_argument("--evalue-col", type=int, default=None,
                    help="1-based column index of e-value in the Diamond TSV. "
                         "If omitted, auto-detected by scanning for sci notation "
                         "(works for the default 12-col outfmt 6 and most "
                         "custom layouts).")

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

    if args.evalue_col is None:
        evalue_col = detect_evalue_col(diamond_tsv)
        print(f"Auto-detected e-value column: {evalue_col + 1} (1-based)",
              flush=True)
    else:
        evalue_col = args.evalue_col - 1  # CLI is 1-based

    print(f"Loading hits at evalue <= {args.evalue:g}", flush=True)
    supported_tids = transcripts_with_hit(diamond_tsv, args.evalue, evalue_col)

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
        "evalue_col": evalue_col + 1,  # 1-based for readability
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
