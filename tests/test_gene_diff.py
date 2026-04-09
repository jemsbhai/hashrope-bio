"""Tests for gene_diff module."""

import tempfile
from pathlib import Path

import pytest
from hashrope import PolynomialHash, Leaf, rope_concat

from hashrope_bio.genomics.gene_diff import (
    Exon, GeneRegion, diff_genes, screen_genes, load_gene_regions,
)


def _build_rope(data: bytes, chunk_size: int, h: PolynomialHash):
    """Build a chunked rope from bytes."""
    rope = None
    for i in range(0, len(data), chunk_size):
        chunk = data[i:i + chunk_size]
        leaf = Leaf(chunk, h)
        rope = rope_concat(rope, leaf, h)
    return rope


# --- Fixtures ---

@pytest.fixture
def h():
    return PolynomialHash()


@pytest.fixture
def ref_seq():
    """1000-byte reference sequence."""
    return b"ACGT" * 250


@pytest.fixture
def sample_same(ref_seq):
    """Sample identical to reference."""
    return ref_seq[:]


@pytest.fixture
def sample_mutated(ref_seq):
    """Sample with mutation at position 500."""
    s = bytearray(ref_seq)
    s[500] = ord('T') if s[500] != ord('T') else ord('C')
    return bytes(s)


@pytest.fixture
def genes():
    """Three test genes spanning the 1000-byte sequence."""
    return [
        GeneRegion(
            name="GENE_A", transcript="NM_001", chrom="chr1",
            strand="+", tx_start=0, tx_end=300,
            exons=[Exon(1, 0, 100), Exon(2, 150, 300)],
        ),
        GeneRegion(
            name="GENE_B", transcript="NM_002", chrom="chr1",
            strand="+", tx_start=400, tx_end=700,
            exons=[Exon(1, 400, 550), Exon(2, 600, 700)],
        ),
        GeneRegion(
            name="GENE_C", transcript="NM_003", chrom="chr1",
            strand="-", tx_start=750, tx_end=1000,
            exons=[Exon(1, 750, 900), Exon(2, 920, 1000)],
        ),
    ]


# --- Tests ---

class TestDiffGenes:
    def test_identical_sequences_no_changes(self, h, ref_seq, sample_same, genes):
        ref_rope = _build_rope(ref_seq, 64, h)
        sample_rope = _build_rope(sample_same, 64, h)
        report = diff_genes(ref_rope, sample_rope, genes, h)

        assert report.total_genes == 3
        assert len(report.changed_genes) == 0
        assert len(report.unchanged_genes) == 3
        assert report.changed_exon_count == 0

    def test_mutation_detected_in_correct_gene(self, h, ref_seq, sample_mutated, genes):
        """Mutation at pos 500 is inside GENE_B (400-700)."""
        ref_rope = _build_rope(ref_seq, 64, h)
        sample_rope = _build_rope(sample_mutated, 64, h)
        report = diff_genes(ref_rope, sample_rope, genes, h)

        changed_names = [gd.gene.name for gd in report.changed_genes]
        assert "GENE_B" in changed_names
        assert "GENE_A" not in changed_names
        assert "GENE_C" not in changed_names

    def test_mutation_detected_in_correct_exon(self, h, ref_seq, sample_mutated, genes):
        """Mutation at pos 500 is inside GENE_B exon 1 (400-550)."""
        ref_rope = _build_rope(ref_seq, 64, h)
        sample_rope = _build_rope(sample_mutated, 64, h)
        report = diff_genes(ref_rope, sample_rope, genes, h)

        gene_b_diff = [gd for gd in report.gene_diffs if gd.gene.name == "GENE_B"][0]
        changed_exon_indices = [ed.exon.index for ed in gene_b_diff.changed_exons]
        assert 1 in changed_exon_indices  # exon 1 (400-550) contains pos 500
        assert 2 not in changed_exon_indices  # exon 2 (600-700) is unchanged

    def test_no_exons_when_disabled(self, h, ref_seq, sample_mutated, genes):
        ref_rope = _build_rope(ref_seq, 64, h)
        sample_rope = _build_rope(sample_mutated, 64, h)
        report = diff_genes(ref_rope, sample_rope, genes, h, include_exons=False)

        assert report.total_exons == 0
        assert len(report.changed_genes) == 1  # gene body still detected

    def test_seq_len_filter(self, h, ref_seq, sample_same, genes):
        ref_rope = _build_rope(ref_seq, 64, h)
        sample_rope = _build_rope(sample_same, 64, h)
        # Only genes within first 500 bp
        report = diff_genes(ref_rope, sample_rope, genes, h, seq_len=500)

        # GENE_A (0-300) included, GENE_B (400-700) excluded, GENE_C (750-1000) excluded
        assert report.total_genes == 1
        assert report.gene_diffs[0].gene.name == "GENE_A"


class TestScreenGenes:
    def test_screen_identical(self, h, ref_seq, sample_same, genes):
        ref_rope = _build_rope(ref_seq, 64, h)
        sample_rope = _build_rope(sample_same, 64, h)
        changed = screen_genes(ref_rope, sample_rope, genes, h)
        assert changed == []

    def test_screen_mutation(self, h, ref_seq, sample_mutated, genes):
        ref_rope = _build_rope(ref_seq, 64, h)
        sample_rope = _build_rope(sample_mutated, 64, h)
        changed = screen_genes(ref_rope, sample_rope, genes, h)
        assert changed == ["GENE_B"]


class TestSummary:
    def test_summary_format(self, h, ref_seq, sample_mutated, genes):
        ref_rope = _build_rope(ref_seq, 64, h)
        sample_rope = _build_rope(sample_mutated, 64, h)
        report = diff_genes(ref_rope, sample_rope, genes, h)
        summary = report.summary()

        assert "1/3 genes changed" in summary
        assert "GENE_B" in summary


class TestLoadGeneRegions:
    def test_load_from_tsv(self, tmp_path):
        tsv = tmp_path / "test_genes.tsv"
        tsv.write_text(
            "gene_name\ttranscript\tstrand\ttx_start\ttx_end\ttx_len\t"
            "cds_start\tcds_end\texon_count\texon_starts\texon_ends\texon_lengths\n"
            "BRCA1\tNM_001\t+\t100\t500\t400\t150\t450\t2\t100,300\t200,500\t100,200\n"
            "TP53\tNM_002\t-\t1000\t2000\t1000\t1100\t1900\t3\t1000,1300,1700\t1200,1500,2000\t200,200,300\n"
        )
        genes = load_gene_regions(tsv)
        assert len(genes) == 2
        assert genes[0].name == "BRCA1"
        assert genes[0].exon_count == 2
        assert genes[0].exons[0].start == 100
        assert genes[0].exons[0].end == 200
        assert genes[1].name == "TP53"
        assert genes[1].exon_count == 3
