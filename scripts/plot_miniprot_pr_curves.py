#!/usr/bin/env python3
"""Plot PR curves from the miniprothint-based filter analysis.

Reads *_miniprot_pr_curves.tsv files and produces:

  miniprot_pr_bars.pdf  —  grouped bar chart: precision and recall per rule,
                            grouped by exclusion level, one subplot per species.
  miniprot_pr_table.tsv —  flat summary table for all (species, level, rule).

Usage:
    python scripts/plot_miniprot_pr_curves.py \\
        --results-dir results/miniprot_analysis/ \\
        --out-dir     results/figures/
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

SPECIES_LABEL = {
    "Arabidopsis_thaliana":      "A. thaliana",
    "Drosophila_melanogaster":   "D. melanogaster",
    "Homo_sapiens":              "H. sapiens",
    "Phaeodactylum_tricornutum": "P. tricornutum",
}
LEVEL_ORDER = ["genus", "order", "subphylum"]
LEVEL_COLOR = {"genus": "#1f77b4", "order": "#ff7f0e", "subphylum": "#2ca02c"}

# Display order and short names for filter rules
RULE_ORDER = [
    "no_filter",
    "all_any_overlap",
    "all_junction_support",
    "all_compatible",
    "all_exact_match",
    "training_any_overlap",
    "training_junction_support",
    "training_compatible",
    "training_exact_match",
]
RULE_LABEL = {
    "no_filter":                 "no filter",
    "all_any_overlap":           "any overlap\n(all)",
    "all_junction_support":      "junc. support\n(all)",
    "all_compatible":            "compatible\n(all)",
    "all_exact_match":           "exact match\n(all)",
    "training_any_overlap":      "any overlap\n(training)",
    "training_junction_support": "junc. support\n(training)",
    "training_compatible":       "compatible\n(training)",
    "training_exact_match":      "exact match\n(training)",
}


def load_results(results_dir: Path) -> pd.DataFrame:
    dfs: list[pd.DataFrame] = []
    pattern = re.compile(r"^(.+)_(genus|order|subphylum)_miniprot_pr_curves\.tsv$")
    for tsv in sorted(results_dir.glob("*_miniprot_pr_curves.tsv")):
        m = pattern.match(tsv.name)
        if not m:
            print(f"  skip: {tsv.name}")
            continue
        df = pd.read_csv(tsv, sep="\t")
        df["species"] = m.group(1)
        df["level"]   = m.group(2)
        dfs.append(df)
    if not dfs:
        raise FileNotFoundError(f"No *_miniprot_pr_curves.tsv files in {results_dir}")
    combined = pd.concat(dfs, ignore_index=True)
    combined["precision"] = combined["precision"].astype(float)
    combined["recall"]    = combined["recall"].astype(float)
    combined["f1"]        = combined["f1"].astype(float)
    return combined


def plot_bars(data: pd.DataFrame, out_dir: Path) -> None:
    species_list = sorted(SPECIES_LABEL)
    rules = [r for r in RULE_ORDER if r in data["filter_rule"].unique()]
    n_rules = len(rules)
    bar_w = 0.25
    offsets = {lv: (i - 1) * bar_w for i, lv in enumerate(LEVEL_ORDER)}
    x = range(n_rules)

    fig, axes = plt.subplots(
        len(species_list), 2,
        figsize=(max(12, n_rules * 1.1), 4 * len(species_list)),
        sharey="col",
    )

    for row, species in enumerate(species_list):
        for col, metric in enumerate(("precision", "recall")):
            ax = axes[row][col]
            sdata = data[data["species"] == species]

            for level in LEVEL_ORDER:
                ldata = sdata[sdata["level"] == level]
                vals = []
                for rule in rules:
                    rdata = ldata[ldata["filter_rule"] == rule][metric]
                    vals.append(float(rdata.iloc[0]) if not rdata.empty else float("nan"))
                ax.bar(
                    [xi + offsets[level] for xi in x],
                    vals,
                    width=bar_w,
                    color=LEVEL_COLOR[level],
                    alpha=0.85,
                    label=level.capitalize() if col == 0 else "_",
                )

            ax.set_xticks(list(x))
            ax.set_xticklabels(
                [RULE_LABEL.get(r, r) for r in rules],
                rotation=35, ha="right", fontsize=7,
            )
            ax.set_ylim(0, 1.05)
            ax.set_ylabel(metric.capitalize(), fontsize=9)
            if col == 0:
                ax.set_title(
                    f"{SPECIES_LABEL.get(species, species)} — {metric.capitalize()}",
                    fontsize=10,
                )
            else:
                ax.set_title(f"{metric.capitalize()}", fontsize=10)
            ax.axhline(
                float(sdata[sdata["filter_rule"] == "no_filter"][metric].mean()),
                color="black", linestyle="--", linewidth=0.8, alpha=0.6,
            )
            ax.grid(axis="y", alpha=0.3, linewidth=0.5)

    level_handles = [
        mpatches.Patch(color=LEVEL_COLOR[lv], label=lv.capitalize())
        for lv in LEVEL_ORDER
    ]
    baseline_handle = plt.Line2D(
        [], [], color="black", linestyle="--", linewidth=0.8, label="Baseline (no filter)",
    )
    fig.legend(
        handles=level_handles + [baseline_handle],
        loc="lower center", ncol=4, fontsize=9,
        frameon=True, bbox_to_anchor=(0.5, 0.0),
    )
    fig.suptitle(
        "Miniprothint-based filtering: precision and recall per rule × exclusion level",
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0.06, 1, 1])

    for ext in ("pdf", "png"):
        path = out_dir / f"miniprot_pr_bars.{ext}"
        fig.savefig(path, dpi=150, bbox_inches="tight")
        print(f"Saved {path}")
    plt.close(fig)


def save_table(data: pd.DataFrame, out_dir: Path) -> None:
    cols = ["species", "level", "filter_rule", "precision", "recall", "f1",
            "n_kept", "n_total"]
    tbl = data[[c for c in cols if c in data.columns]].copy()
    tbl = tbl.sort_values(["species", "level", "filter_rule"])
    path = out_dir / "miniprot_pr_table.tsv"
    tbl.to_csv(path, sep="\t", index=False)
    print(f"Saved {path}")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", type=Path, default=Path("results/miniprot_analysis"),
                    help="Directory with *_miniprot_pr_curves.tsv files")
    ap.add_argument("--out-dir", type=Path, default=Path("results/figures"),
                    help="Output directory for figures")
    args = ap.parse_args(argv)

    args.out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading results from {args.results_dir} ...", flush=True)
    data = load_results(args.results_dir)
    n_combos = data.groupby(["species", "level"]).ngroups
    print(f"  {len(data)} rows across {n_combos} species × level combinations", flush=True)

    plot_bars(data, args.out_dir)
    save_table(data, args.out_dir)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
