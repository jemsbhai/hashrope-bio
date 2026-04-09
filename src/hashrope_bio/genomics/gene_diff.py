"""Gene-level change detection via hash comparison.

Given two ropes (e.g., reference and patient genome) and a list of
gene annotations, reports which genes and exons differ — in O(N · log w)
total time, where N is the number of regions and w is the rope weight.

This implements a two-pass architecture:
  Pass 1 (screen): O(log w) hash comparison per region — identifies changes
  Pass 2 (detail): only changed regions need expensive byte-level analysis

Typical usage:
    ref_rope, h, _ = load_fasta_to_rope("reference.fa")
    sample_rope, _, _ = load_fasta_to_rope("patient.fa")
    genes = load_gene_regions("chr22_genes.tsv")
    report = diff_genes(ref_rope, sample_rope, genes, h)
    for g in report.changed_genes:
        print(f"{g.name}: gene body changed, {len(g.changed_exons)}/{g.exon_count} exons differ")
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

from hashrope import PolynomialHash, rope_substr_hash


# ---------------------------------------------------------------------------
# Data types
# ---------------------------------------------------------------------------

@dataclass
class Exon:
    """A single exon with genomic coordinates (0-based, half-open)."""
    index: int       # 1-based exon number
    start: int       # 0-based genomic start
    end: int         # 0-based genomic end (exclusive)

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass
class GeneRegion:
    """A gene with its transcript boundaries and exon structure."""
    name: str
    transcript: str
    chrom: str
    strand: str
    tx_start: int    # 0-based genomic start
    tx_end: int      # 0-based genomic end (exclusive)
    exons: list[Exon] = field(default_factory=list)

    @property
    def tx_length(self) -> int:
        return self.tx_end - self.tx_start

    @property
    def exon_count(self) -> int:
        return len(self.exons)


@dataclass
class ExonDiff:
    """Result of comparing one exon between two ropes."""
    exon: Exon
    ref_hash: int
    sample_hash: int

    @property
    def changed(self) -> bool:
        return self.ref_hash != self.sample_hash


@dataclass
class GeneDiff:
    """Result of comparing one gene (body + exons) between two ropes."""
    gene: GeneRegion
    body_ref_hash: int
    body_sample_hash: int
    exon_diffs: list[ExonDiff] = field(default_factory=list)

    @property
    def body_changed(self) -> bool:
        return self.body_ref_hash != self.body_sample_hash

    @property
    def changed_exons(self) -> list[ExonDiff]:
        return [d for d in self.exon_diffs if d.changed]

    @property
    def unchanged_exons(self) -> list[ExonDiff]:
        return [d for d in self.exon_diffs if not d.changed]


@dataclass
class GeneDiffReport:
    """Full report from diffing a set of genes between two ropes."""
    gene_diffs: list[GeneDiff]

    @property
    def total_genes(self) -> int:
        return len(self.gene_diffs)

    @property
    def changed_genes(self) -> list[GeneDiff]:
        return [d for d in self.gene_diffs if d.body_changed]

    @property
    def unchanged_genes(self) -> list[GeneDiff]:
        return [d for d in self.gene_diffs if not d.body_changed]

    @property
    def total_exons(self) -> int:
        return sum(len(d.exon_diffs) for d in self.gene_diffs)

    @property
    def changed_exon_count(self) -> int:
        return sum(len(d.changed_exons) for d in self.gene_diffs)

    def summary(self) -> str:
        """Human-readable summary."""
        lines = [
            f"Gene diff: {len(self.changed_genes)}/{self.total_genes} genes changed, "
            f"{self.changed_exon_count}/{self.total_exons} exons changed",
        ]
        for gd in self.changed_genes:
            ce = len(gd.changed_exons)
            lines.append(
                f"  {gd.gene.name} ({gd.gene.tx_length:,} bp): "
                f"body changed, {ce}/{gd.gene.exon_count} exons differ"
            )
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Gene loading
# ---------------------------------------------------------------------------

def load_gene_regions(
    path: str | Path,
    chrom_filter: str | None = None,
) -> list[GeneRegion]:
    """Load gene annotations from TSV (chr22_genes.tsv format).

    Expected columns: gene_name, transcript, strand, tx_start, tx_end,
    tx_len, cds_start, cds_end, exon_count, exon_starts, exon_ends, exon_lengths

    Args:
        path: Path to TSV file.
        chrom_filter: If set, only load genes on this chromosome.

    Returns:
        List of GeneRegion, sorted by tx_start.
    """
    genes = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            name = row["gene_name"]
            transcript = row["transcript"]
            strand = row["strand"]
            tx_start = int(row["tx_start"])
            tx_end = int(row["tx_end"])

            exon_starts = [int(x) for x in row["exon_starts"].split(",") if x]
            exon_ends = [int(x) for x in row["exon_ends"].split(",") if x]

            exons = [
                Exon(index=i + 1, start=s, end=e)
                for i, (s, e) in enumerate(zip(exon_starts, exon_ends))
            ]

            genes.append(GeneRegion(
                name=name, transcript=transcript, chrom="chr22",
                strand=strand, tx_start=tx_start, tx_end=tx_end,
                exons=exons,
            ))

    genes.sort(key=lambda g: g.tx_start)
    return genes


# ---------------------------------------------------------------------------
# Core: gene-level diff
# ---------------------------------------------------------------------------

def diff_genes(
    ref_rope,
    sample_rope,
    genes: list[GeneRegion],
    h: PolynomialHash,
    *,
    include_exons: bool = True,
    seq_len: int | None = None,
) -> GeneDiffReport:
    """Compare gene regions between two ropes.

    For each gene, compares the gene body hash. If include_exons is True,
    also compares each exon individually (useful for pinpointing which
    exons changed within a changed gene).

    The total cost is O((G + E) · log w) where G is the number of genes,
    E is the number of exons, and w is the rope weight. This is dramatically
    cheaper than O(G · L_avg) byte-by-byte comparison for large genes.

    Args:
        ref_rope: Reference genome rope.
        sample_rope: Sample/patient genome rope (same coordinate system).
        genes: List of GeneRegion to compare.
        h: Shared PolynomialHash instance.
        include_exons: Also diff individual exons (default True).
        seq_len: If known, skip genes beyond this position.

    Returns:
        GeneDiffReport with per-gene and per-exon results.
    """
    gene_diffs = []

    for gene in genes:
        # Skip genes outside sequence bounds
        if seq_len is not None and gene.tx_end > seq_len:
            continue

        # Gene body comparison: O(log w) per rope
        body_ref = rope_substr_hash(ref_rope, gene.tx_start, gene.tx_length, h)
        body_sample = rope_substr_hash(sample_rope, gene.tx_start, gene.tx_length, h)

        exon_diffs = []
        if include_exons:
            for exon in gene.exons:
                if seq_len is not None and exon.end > seq_len:
                    continue
                if exon.length == 0:
                    continue
                ex_ref = rope_substr_hash(ref_rope, exon.start, exon.length, h)
                ex_sample = rope_substr_hash(sample_rope, exon.start, exon.length, h)
                exon_diffs.append(ExonDiff(exon=exon, ref_hash=ex_ref, sample_hash=ex_sample))

        gene_diffs.append(GeneDiff(
            gene=gene,
            body_ref_hash=body_ref,
            body_sample_hash=body_sample,
            exon_diffs=exon_diffs,
        ))

    return GeneDiffReport(gene_diffs=gene_diffs)


def screen_genes(
    ref_rope,
    sample_rope,
    genes: list[GeneRegion],
    h: PolynomialHash,
    *,
    seq_len: int | None = None,
) -> list[str]:
    """Quick screen: return names of genes whose bodies differ.

    Gene-body-only comparison (no exon detail). Fastest possible
    screening pass — O(G · log w) total.

    Args:
        ref_rope: Reference genome rope.
        sample_rope: Sample genome rope.
        genes: List of GeneRegion.
        h: Shared PolynomialHash instance.
        seq_len: If known, skip genes beyond this position.

    Returns:
        List of gene names where the body hash differs.
    """
    changed = []
    for gene in genes:
        if seq_len is not None and gene.tx_end > seq_len:
            continue
        h_ref = rope_substr_hash(ref_rope, gene.tx_start, gene.tx_length, h)
        h_sample = rope_substr_hash(sample_rope, gene.tx_start, gene.tx_length, h)
        if h_ref != h_sample:
            changed.append(gene.name)
    return changed
