#!/usr/bin/env python3
"""Plot precision-recall curves from the Diamond filter-rule analysis.

Reads all *_pr_curves.tsv files produced by analyze_filter_rules.py and
saves two figures:

  pr_curves_by_rule.pdf  —  6-panel PR curves (one per rule family)
  has_any_hit_summary.pdf —  precision / recall gain from the simplest
                             filter (any Diamond hit) across levels

Usage:
    python scripts/plot_pr_curves.py \
        --results-dir results/analysis/ \
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
import matplotlib.lines as mlines
import matplotlib.patches as mpatches

# ── Colours / styles ──────────────────────────────────────────────────────────
SPECIES_COLOR = {
    "Arabidopsis_thaliana":      "#2ca02c",
    "Drosophila_melanogaster":   "#ff7f0e",
    "Homo_sapiens":              "#1f77b4",
    "Phaeodactylum_tricornutum": "#9467bd",
}
SPECIES_LABEL = {
    "Arabidopsis_thaliana":      "A. thaliana",
    "Drosophila_melanogaster":   "D. melanogaster",
    "Homo_sapiens":              "H. sapiens",
    "Phaeodactylum_tricornutum": "P. tricornutum",
}
LEVEL_STYLE  = {"genus": "-",  "order": "--", "subphylum": ":"}
LEVEL_ORDER  = ["genus", "order", "subphylum"]

# Rule families → (title, x-axis label, whether higher threshold = stricter)
RULE_FAMILIES: list[tuple[str, str]] = [
    ("pident_ge",        "Min. % identity"),
    ("qcov_ge",          "Min. query coverage (%)"),
    ("bitscore_ge",      "Min. bitscore"),
    ("evalue_le_1e",     "Max. e-value  (10^threshold)"),
    ("pident_ge_qcov50", "Min. % identity  (qcov ≥ 50 %)"),
    ("pident_ge_qcov70", "Min. % identity  (qcov ≥ 70 %)"),
]


# ── I/O ───────────────────────────────────────────────────────────────────────
def load_results(results_dir: Path) -> pd.DataFrame:
    dfs: list[pd.DataFrame] = []
    pattern = re.compile(r"^(.+)_(genus|order|subphylum)_pr_curves\.tsv$")
    for tsv in sorted(results_dir.glob("*_pr_curves.tsv")):
        m = pattern.match(tsv.name)
        if not m:
            print(f"  skip: {tsv.name}")
            continue
        df = pd.read_csv(tsv, sep="\t")
        df["species"] = m.group(1)
        df["level"]   = m.group(2)
        dfs.append(df)
    if not dfs:
        raise FileNotFoundError(f"No *_pr_curves.tsv files in {results_dir}")
    combined = pd.concat(dfs, ignore_index=True)
    combined["precision"] = combined["precision"].astype(float)
    combined["recall"]    = combined["recall"].astype(float)
    combined["threshold"] = combined["threshold"].astype(float)
    return combined


# ── Figure 1: PR curves ───────────────────────────────────────────────────────
def _draw_panel(ax: plt.Axes, data: pd.DataFrame, rule: str, title: str) -> None:
    subset = data[data["filter_rule"] == rule]
    if subset.empty:
        ax.set_visible(False)
        return

    for species in sorted(SPECIES_COLOR):
        color = SPECIES_COLOR[species]
        for level in LEVEL_ORDER:
            rows = (
                subset[(subset["species"] == species) & (subset["level"] == level)]
                .sort_values("threshold")
            )
            if rows.empty:
                continue
            ax.plot(
                rows["recall"], rows["precision"],
                color=color, linestyle=LEVEL_STYLE[level],
                linewidth=1.6, alpha=0.85,
            )

    # Baseline cross for each species (same regardless of level)
    for species in sorted(SPECIES_COLOR):
        base = data[
            (data["filter_rule"] == "no_filter") & (data["species"] == species)
        ]
        if base.empty:
            continue
        ax.scatter(
            base["recall"].mean(), base["precision"].mean(),
            color=SPECIES_COLOR[species], marker="x", s=70,
            linewidths=2, zorder=5,
        )

    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(0.5, 1.02)
    ax.set_xlabel("Recall", fontsize=9)
    ax.set_ylabel("Precision", fontsize=9)
    ax.set_title(title, fontsize=10)
    ax.grid(True, alpha=0.3, linewidth=0.5)


def _make_legend(fig: plt.Figure) -> None:
    species_handles = [
        mpatches.Patch(color=c, label=SPECIES_LABEL[s])
        for s, c in sorted(SPECIES_COLOR.items())
    ]
    level_handles = [
        mlines.Line2D([], [], color="gray", linestyle=ls,
                      linewidth=1.6, label=lv.capitalize())
        for lv, ls in LEVEL_STYLE.items()
    ]
    baseline_handle = mlines.Line2D(
        [], [], color="gray", marker="x", linestyle="None",
        markersize=8, markeredgewidth=2, label="No filter (baseline)",
    )
    fig.legend(
        handles=species_handles + [plt.Line2D([], [], visible=False)]
               + level_handles + [baseline_handle],
        loc="lower center", ncol=5, fontsize=8,
        frameon=True, bbox_to_anchor=(0.5, 0.0),
    )


def plot_pr_figure(data: pd.DataFrame, out_dir: Path) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    fig.suptitle(
        "Precision – Recall under Diamond filter rules\n"
        "(lines = ODB exclusion level;  × = no-filter baseline)",
        fontsize=12,
    )
    for ax, (rule, title) in zip(axes.flat, RULE_FAMILIES):
        _draw_panel(ax, data, rule, title)
    _make_legend(fig)
    fig.tight_layout(rect=[0, 0.07, 1, 1])
    for ext in ("pdf", "png"):
        path = out_dir / f"pr_curves_by_rule.{ext}"
        fig.savefig(path, dpi=150, bbox_inches="tight")
        print(f"Saved {path}")
    plt.close(fig)


# ── Figure 2: has_any_hit summary ────────────────────────────────────────────
def plot_has_any_hit_figure(data: pd.DataFrame, out_dir: Path) -> None:
    """Grouped bar chart: baseline precision/recall vs. 'has any hit', per level."""
    species_list = sorted(SPECIES_COLOR)
    n = len(species_list)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(
        "Effect of requiring any Diamond hit (per species × exclusion level)\n"
        "Dark bars = filtered; light bars = baseline",
        fontsize=11,
    )

    bar_w = 0.18
    offsets = {"genus": -bar_w, "order": 0.0, "subphylum": bar_w}

    for ax, metric in zip(axes, ("precision", "recall")):
        for i, species in enumerate(species_list):
            color = SPECIES_COLOR[species]
            # baseline (level-independent; same for all levels, just use first)
            base_val = data[
                (data["filter_rule"] == "no_filter") & (data["species"] == species)
            ][metric].mean()
            ax.bar(i, base_val, width=bar_w * 3 + 0.04,
                   color=color, alpha=0.20, zorder=1)

            for level in LEVEL_ORDER:
                val = data[
                    (data["filter_rule"] == "has_any_hit") &
                    (data["species"] == species) &
                    (data["level"] == level)
                ][metric]
                if val.empty:
                    continue
                ax.bar(
                    i + offsets[level], float(val.iloc[0]),
                    width=bar_w, color=color, alpha=0.85,
                    hatch={"genus": "", "order": "//", "subphylum": "xx"}[level],
                    edgecolor="white", zorder=2,
                )

        ax.set_xticks(range(n))
        ax.set_xticklabels(
            [SPECIES_LABEL[s] for s in species_list],
            rotation=25, ha="right", fontsize=9,
        )
        ax.set_ylabel(metric.capitalize(), fontsize=10)
        ax.set_ylim(0.5, 1.02)
        ax.set_title(metric.capitalize(), fontsize=10)
        ax.grid(axis="y", alpha=0.3, linewidth=0.5)

    # legend for levels
    level_handles = [
        mpatches.Patch(facecolor="gray", alpha=0.85,
                       hatch={"genus": "", "order": "//", "subphylum": "xx"}[lv],
                       edgecolor="white", label=lv.capitalize())
        for lv in LEVEL_ORDER
    ]
    baseline_handle = mpatches.Patch(facecolor="gray", alpha=0.20, label="Baseline (no filter)")
    fig.legend(
        handles=level_handles + [baseline_handle],
        loc="lower center", ncol=4, fontsize=9,
        frameon=True, bbox_to_anchor=(0.5, 0.0),
    )
    fig.tight_layout(rect=[0, 0.08, 1, 1])
    for ext in ("pdf", "png"):
        path = out_dir / f"has_any_hit_summary.{ext}"
        fig.savefig(path, dpi=150, bbox_inches="tight")
        print(f"Saved {path}")
    plt.close(fig)


# ── Figure 3: F1 heatmap ──────────────────────────────────────────────────────
def plot_f1_heatmap(data: pd.DataFrame, out_dir: Path) -> None:
    """Best F1 score per (species, level, rule) as a heatmap."""
    records = []
    for (species, level, rule), grp in data.groupby(["species", "level", "filter_rule"]):
        if rule in ("no_filter",):
            continue
        best_f1 = grp["f1"].astype(float).max()
        records.append({"species": species, "level": level, "rule": rule, "best_f1": best_f1})
    if not records:
        return

    df = pd.DataFrame(records)
    pivot = df.pivot_table(index=["species", "level"], columns="rule", values="best_f1")

    fig, ax = plt.subplots(figsize=(14, 6))
    im = ax.imshow(pivot.values, aspect="auto", cmap="YlGn", vmin=0.7, vmax=1.0)
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=40, ha="right", fontsize=8)
    row_labels = [
        f"{SPECIES_LABEL[s]}  [{l}]"
        for s, l in pivot.index
    ]
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=8)
    plt.colorbar(im, ax=ax, label="Best F1", shrink=0.7)
    ax.set_title("Best F1 score per filter rule × species × exclusion level", fontsize=11)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        path = out_dir / f"f1_heatmap.{ext}"
        fig.savefig(path, dpi=150, bbox_inches="tight")
        print(f"Saved {path}")
    plt.close(fig)


# ── Main ──────────────────────────────────────────────────────────────────────
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results-dir", type=Path, default=Path("results/analysis"),
                    help="Directory with *_pr_curves.tsv files")
    ap.add_argument("--out-dir", type=Path, default=Path("results/figures"),
                    help="Output directory for figures")
    args = ap.parse_args(argv)

    args.out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading results from {args.results_dir} ...", flush=True)
    data = load_results(args.results_dir)
    n_combos = data.groupby(["species", "level"]).ngroups
    print(f"  {len(data)} rows across {n_combos} species × level combinations", flush=True)

    plot_pr_figure(data, args.out_dir)
    plot_has_any_hit_figure(data, args.out_dir)
    plot_f1_heatmap(data, args.out_dir)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
