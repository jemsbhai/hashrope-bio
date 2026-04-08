"""Tests for genomics module — synthetic data, no downloads needed."""

import pytest
from hashrope import PolynomialHash, Leaf, rope_concat, rope_hash, rope_substr_hash

from hashrope_bio.genomics.region_query import region_hash, regions_identical
from hashrope_bio.genomics.mutation import localize_mutation, localize_all_mutations
from hashrope_bio.genomics.repeats import build_repeat_node, detect_tandem_repeat


@pytest.fixture
def h():
    return PolynomialHash()


class TestRegionQuery:
    def test_region_hash_matches_direct(self, h):
        data = b"ACGTACGTACGTACGT"
        rope = Leaf(data, h)
        assert region_hash(rope, 4, 8, h) == h.hash(b"ACGTACGT")

    def test_regions_identical_same_data(self, h):
        data = b"ACGTACGTACGTACGT"
        rope_a = Leaf(data, h)
        rope_b = Leaf(data, h)
        assert regions_identical(rope_a, rope_b, 0, 16, h)

    def test_regions_identical_detects_snp(self, h):
        ref = b"ACGTACGTACGTACGT"
        mut = b"ACGTACCTACGTACGT"  # G->C at position 6
        rope_ref = Leaf(ref, h)
        rope_mut = Leaf(mut, h)
        assert not regions_identical(rope_ref, rope_mut, 0, 16, h)
        # But flanking regions are still identical
        assert regions_identical(rope_ref, rope_mut, 0, 6, h)
        assert regions_identical(rope_ref, rope_mut, 7, 9, h)

    def test_multi_leaf_rope(self, h):
        """Region query across internal node boundaries."""
        l1 = Leaf(b"ACGT", h)
        l2 = Leaf(b"TGCA", h)
        rope = rope_concat(l1, l2, h)
        # Spanning query
        assert region_hash(rope, 2, 4, h) == h.hash(b"GTTG")


class TestMutationLocalization:
    def test_single_snp(self, h):
        ref = b"ACGTACGTACGTACGT"
        mut = bytearray(ref)
        mut[10] = ord("T")  # G->T at position 10
        mut = bytes(mut)

        rope_ref = Leaf(ref, h)
        rope_mut = Leaf(mut, h)

        pos = localize_mutation(rope_ref, rope_mut, 0, 16, h)
        assert pos == 10

    def test_no_mutation_returns_none(self, h):
        data = b"ACGTACGT"
        rope_a = Leaf(data, h)
        rope_b = Leaf(data, h)
        assert localize_mutation(rope_a, rope_b, 0, 8, h) is None

    def test_multiple_mutations(self, h):
        ref = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"  # 30 A's
        mut = bytearray(ref)
        mut[5] = ord("T")
        mut[20] = ord("C")
        mut = bytes(mut)

        rope_ref = Leaf(ref, h)
        rope_mut = Leaf(mut, h)

        regions = localize_all_mutations(rope_ref, rope_mut, 0, 30, h)
        positions = [r[0] for r in regions]
        assert 5 in positions
        assert 20 in positions

    def test_localize_on_multi_leaf_rope(self, h):
        """Mutation localization on a rope built from multiple leaves."""
        ref_chunks = [b"ACGTACGT", b"TGCATGCA", b"AAAACCCC"]
        mut_chunks = [b"ACGTACGT", b"TGCATGCA", b"AAAATCCC"]  # C->T at position 20

        ref_rope = None
        for chunk in ref_chunks:
            leaf = Leaf(chunk, h)
            ref_rope = rope_concat(ref_rope, leaf, h) if ref_rope else leaf

        mut_rope = None
        for chunk in mut_chunks:
            leaf = Leaf(chunk, h)
            mut_rope = rope_concat(mut_rope, leaf, h) if mut_rope else leaf

        pos = localize_mutation(ref_rope, mut_rope, 0, 24, h)
        assert pos == 20


class TestRepeats:
    def test_build_repeat_node_hash(self, h):
        motif = b"CAG"
        count = 100
        node = build_repeat_node(motif, count, h)
        expected = h.hash(motif * count)
        assert rope_hash(node) == expected

    def test_detect_tandem_repeat(self):
        seq = b"CAGCAGCAGCAGCAGXYZ"
        assert detect_tandem_repeat(seq, b"CAG") == 5

    def test_detect_no_repeat(self):
        seq = b"ACGTACGT"
        assert detect_tandem_repeat(seq, b"CAG") is None

    def test_detect_single_repeat(self):
        seq = b"CAGXXX"
        assert detect_tandem_repeat(seq, b"CAG") == 1
