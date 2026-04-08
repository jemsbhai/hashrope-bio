"""Tests for cheminformatics module — synthetic data."""

import pytest
from hashrope import PolynomialHash, Leaf, rope_substr_hash

from hashrope_bio.cheminformatics.resistance import (
    check_resistance_panel, ResistanceSite, HIV_RT_NNRTI_PANEL,
)


@pytest.fixture
def h():
    return PolynomialHash()


class TestResistancePanel:
    def test_wildtype_no_mutations(self, h):
        """Identical ref and sample → no mutations detected."""
        gene = b"A" * 300  # synthetic gene, 100 codons
        ref_rope = Leaf(gene, h)
        sample_rope = Leaf(gene, h)

        panel = [
            ResistanceSite("TEST", 50, "X", "Y", "TestDrug", "X50Y"),
            ResistanceSite("TEST", 100, "X", "Y", "TestDrug", "X100Y"),
        ]

        results = check_resistance_panel(ref_rope, sample_rope, panel, h)
        assert all(not r.is_mutant for r in results)

    def test_single_mutation_detected(self, h):
        """One codon differs → one mutation flagged."""
        ref_gene = b"AAAAAAAAA"   # 3 codons: AAA AAA AAA
        mut_gene = b"AAAGGGAAA"   # 3 codons: AAA GGG AAA (codon 2 mutated)
        ref_rope = Leaf(ref_gene, h)
        mut_rope = Leaf(mut_gene, h)

        panel = [
            ResistanceSite("TEST", 1, "K", "K", "TestDrug", "no change"),
            ResistanceSite("TEST", 2, "K", "G", "TestDrug", "K2G"),
            ResistanceSite("TEST", 3, "K", "K", "TestDrug", "no change"),
        ]

        results = check_resistance_panel(ref_rope, mut_rope, panel, h)
        assert not results[0].is_mutant
        assert results[1].is_mutant
        assert not results[2].is_mutant

    def test_hiv_panel_exists(self):
        """Verify the pre-defined HIV panel has expected entries."""
        assert len(HIV_RT_NNRTI_PANEL) > 0
        assert any(s.position == 103 for s in HIV_RT_NNRTI_PANEL)  # K103N
