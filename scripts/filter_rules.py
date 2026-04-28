"""Filter rule abstractions for protein-based gene filtering.

A `FilterRule` decides whether the best Diamond hit of a query transcript is
"good enough" to keep that transcript.  An `AdaptiveRule` picks one of several
`FilterRule`s based on the *DB-to-genome distance* estimated from the
best-hit %ID distribution and the fraction of queries with any hit.

To change the filter behaviour, edit the defaults below or pass custom rules
into `AdaptiveRule(...)`.

A note on semantics:
    * A transcript with NO Diamond hit is *dropped* by default
      (FilterRule.passes(None) -> False).
    * If the rule is `disabled=True`, *every* transcript is kept regardless
      (no filtering).  This is the safe choice when the DB is too small or
      too distant to be a reliable filter.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import NamedTuple, Optional


# ── Hit / stats ───────────────────────────────────────────────────────────────
class Hit(NamedTuple):
    """Best Diamond hit summary for a single query transcript."""
    pident:   float    # % identity, 0–100
    qcov:     float    # query coverage,  align_len / qlen * 100
    tcov:     float    # target coverage, align_len / slen * 100
    evalue:   float
    bitscore: float


@dataclass
class DBStats:
    """Aggregate statistics about a Diamond run, used by AdaptiveRule."""
    n_queries:     int      # total query transcripts evaluated
    n_with_hit:    int      # how many had at least one hit
    median_pident: float    # median best-hit %ID (over hit-bearing queries)

    @property
    def hit_fraction(self) -> float:
        return self.n_with_hit / self.n_queries if self.n_queries > 0 else 0.0


# ── Per-hit threshold rule ────────────────────────────────────────────────────
@dataclass
class FilterRule:
    """A simple AND-of-thresholds rule applied to one Diamond hit.

    Defaults are permissive (every threshold off).  Set whichever fields you
    want to constrain.  `disabled=True` means: skip filtering entirely
    (keep every transcript, regardless of whether it has a hit).
    """
    name:         str   = "fixed"
    pident_min:   float = 0.0       # %
    qcov_min:     float = 0.0       # %
    tcov_min:     float = 0.0       # %
    evalue_max:   float = float("inf")
    bitscore_min: float = 0.0
    disabled:     bool  = False     # if True, every transcript passes

    def passes(self, hit: Optional[Hit]) -> bool:
        if self.disabled:
            return True
        if hit is None:
            return False                # no hit -> drop
        return (
            hit.pident   >= self.pident_min
            and hit.qcov >= self.qcov_min
            and hit.tcov >= self.tcov_min
            and hit.evalue   <= self.evalue_max
            and hit.bitscore >= self.bitscore_min
        )

    def describe(self) -> str:
        if self.disabled:
            return f"{self.name}(DISABLED — keep everything)"
        parts: list[str] = []
        if self.pident_min   > 0:            parts.append(f"pident≥{self.pident_min:g}")
        if self.qcov_min     > 0:            parts.append(f"qcov≥{self.qcov_min:g}")
        if self.tcov_min     > 0:            parts.append(f"tcov≥{self.tcov_min:g}")
        if self.evalue_max   < float("inf"): parts.append(f"evalue≤{self.evalue_max:g}")
        if self.bitscore_min > 0:            parts.append(f"bitscore≥{self.bitscore_min:g}")
        body = ", ".join(parts) if parts else "any-hit"
        return f"{self.name}({body})"


# ── Adaptive rule selector ────────────────────────────────────────────────────
@dataclass
class AdaptiveRule:
    """Pick a FilterRule based on the DB's apparent distance to the genome.

    Decision tree, top to bottom:

        1. hit_fraction < min_hit_fraction
           -> rule_low_coverage       (DB too small/sparse — don't trust it)
        2. median_pident ≥ close_threshold
           -> rule_close              (close DB — has_any_hit is enough)
        3. median_pident ≥ medium_threshold
           -> rule_medium             (moderate DB — combined thresholds)
        4. otherwise
           -> rule_distant            (distant DB — filtering hurts; keep all)
    """
    close_threshold:  float = 70.0
    medium_threshold: float = 40.0
    min_hit_fraction: float = 0.3

    rule_close: FilterRule = field(default_factory=lambda: FilterRule(
        name="close",
        pident_min=50.0,
        qcov_min=90.0,
        evalue_max=1e-5,
    ))
    rule_medium: FilterRule = field(default_factory=lambda: FilterRule(
        name="medium",
        pident_min=30.0,
        qcov_min=50.0,
        evalue_max=1e-5,
    ))
    rule_distant: FilterRule = field(default_factory=lambda: FilterRule(
        name="distant",
        disabled=True,
    ))
    rule_low_coverage: FilterRule = field(default_factory=lambda: FilterRule(
        name="low_coverage",
        disabled=True,
    ))

    def select(self, stats: DBStats) -> FilterRule:
        if stats.hit_fraction < self.min_hit_fraction:
            return self.rule_low_coverage
        if stats.median_pident >= self.close_threshold:
            return self.rule_close
        if stats.median_pident >= self.medium_threshold:
            return self.rule_medium
        return self.rule_distant

    def describe(self) -> str:
        return (
            f"AdaptiveRule(close≥{self.close_threshold}, "
            f"medium≥{self.medium_threshold}, "
            f"min_hit_frac={self.min_hit_fraction})\n"
            f"  close      : {self.rule_close.describe()}\n"
            f"  medium     : {self.rule_medium.describe()}\n"
            f"  distant    : {self.rule_distant.describe()}\n"
            f"  low-cov    : {self.rule_low_coverage.describe()}"
        )


# Convenient default singleton — import from elsewhere as
#     from filter_rules import DEFAULT_ADAPTIVE_RULE
DEFAULT_ADAPTIVE_RULE = AdaptiveRule()
