"""Tests for genomics module — synthetic data, no downloads needed."""

import tempfile
from pathlib import Path

import pytest
from hashrope import PolynomialHash, Leaf, rope_concat, rope_hash, rope_substr_hash, rope_len

from hashrope_bio.genomics.fasta import load_fasta_to_rope, load_fasta_chunked, load_fasta_bytes
from hashrope_bio.genomics.region_query import region_hash, regions_identical
from hashrope_bio.genomics.mutation import localize_mutation, localize_all_mutations
from hashrope_bio.genomics.repeats import build_repeat_node, detect_tandem_repeat


@pytest.fixture
def h():
    return PolynomialHash()


@pytest.fixture
def sample_fasta(tmp_path):
    """Create a small synthetic FASTA file for testing."""
    seq = "ACGTACGTACGTACGT" * 100  # 1600 bp
    fasta = tmp_path / "test.fa"
    fasta.write_text(f">chr_test synthetic test sequence\n{seq}\n")
    return fasta, seq.encode("ascii")


@pytest.fixture
def multi_line_fasta(tmp_path):
    """FASTA with sequence split across multiple lines (standard format)."""
    seq = "ACGTACGTACGTACGT" * 100  # 1600 bp
    lines = [seq[i:i+80] for i in range(0, len(seq), 80)]  # 80 chars per line
    content = ">chr_multi multi-line test\n" + "\n".join(lines) + "\n"
    fasta = tmp_path / "multiline.fa"
    fasta.write_text(content)
    return fasta, seq.encode("ascii")


@pytest.fixture
def two_seq_fasta(tmp_path):
    """FASTA with two sequences — only the first should be loaded."""
    seq1 = "AAAA" * 100
    seq2 = "TTTT" * 100
    content = f">seq1 first\n{seq1}\n>seq2 second\n{seq2}\n"
    fasta = tmp_path / "two_seq.fa"
    fasta.write_text(content)
    return fasta, seq1.encode("ascii")


class TestFastaLoading:
    def test_load_basic(self, sample_fasta, h):
        fasta_path, expected_bytes = sample_fasta
        rope, hasher, meta = load_fasta_to_rope(fasta_path, h=h)

        assert rope is not None
        assert meta["seq_name"] == "chr_test"
        assert meta["seq_len"] == 1600
        assert meta["backend"] == "pure_python"
        assert meta["chunk_count"] > 0
        assert meta["construction_time_s"] >= 0

        # Hash must match direct hash of the sequence
        assert rope_hash(rope) == h.hash(expected_bytes)

    def test_load_multi_line(self, multi_line_fasta, h):
        fasta_path, expected_bytes = multi_line_fasta
        rope, hasher, meta = load_fasta_to_rope(fasta_path, h=h)

        assert meta["seq_len"] == 1600
        assert rope_hash(rope) == h.hash(expected_bytes)

    def test_load_two_sequences_takes_first(self, two_seq_fasta, h):
        fasta_path, expected_bytes = two_seq_fasta
        rope, hasher, meta = load_fasta_to_rope(fasta_path, h=h)

        assert meta["seq_name"] == "seq1"
        assert meta["seq_len"] == 400
        assert rope_hash(rope) == h.hash(expected_bytes)

    def test_chunk_sizes(self, sample_fasta, h):
        fasta_path, expected_bytes = sample_fasta
        expected_hash = h.hash(expected_bytes)

        for cs in [64, 256, 1024, 4096]:
            rope, _, meta = load_fasta_to_rope(fasta_path, chunk_size=cs, h=h)
            assert rope_hash(rope) == expected_hash, f"Hash mismatch at chunk_size={cs}"
            assert meta["chunk_count"] == -(-1600 // cs)  # ceil division

    def test_chunked_iterator(self, sample_fasta):
        fasta_path, expected_bytes = sample_fasta
        chunks = list(load_fasta_chunked(fasta_path, chunk_size=256))

        assert len(chunks) == -(-1600 // 256)  # 7 chunks (6×256 + 1×64... wait, 1600/256=6.25)
        reassembled = b"".join(chunks)
        assert reassembled == expected_bytes

    def test_load_fasta_bytes(self, sample_fasta):
        fasta_path, expected_bytes = sample_fasta
        seq_bytes, name = load_fasta_bytes(fasta_path)
        assert seq_bytes == expected_bytes
        assert name == "chr_test"

    def test_uppercase_normalization(self, tmp_path, h):
        """Lowercase FASTA should be uppercased by default."""
        seq = "acgtacgt"
        fasta = tmp_path / "lower.fa"
        fasta.write_text(f">test\n{seq}\n")

        rope, _, meta = load_fasta_to_rope(fasta, h=h)
        assert rope_hash(rope) == h.hash(b"ACGTACGT")

    def test_rope_len_matches_seq_len(self, sample_fasta, h):
        fasta_path, _ = sample_fasta
        rope, _, meta = load_fasta_to_rope(fasta_path, h=h)
        assert rope_len(rope) == meta["seq_len"]

    def test_substr_hash_on_loaded_rope(self, sample_fasta, h):
        """Verify substr_hash works on a FASTA-loaded rope."""
        fasta_path, expected_bytes = sample_fasta
        rope, hasher, _ = load_fasta_to_rope(fasta_path, chunk_size=256, h=h)

        # Check various substrings
        for start, length in [(0, 16), (100, 50), (800, 200), (1500, 100)]:
            expected = h.hash(expected_bytes[start:start + length])
            got = rope_substr_hash(rope, start, length, h)
            assert got == expected, f"substr_hash mismatch at start={start}, length={length}"


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
        assert regions_identical(rope_ref, rope_mut, 0, 6, h)
        assert regions_identical(rope_ref, rope_mut, 7, 9, h)

    def test_multi_leaf_rope(self, h):
        l1 = Leaf(b"ACGT", h)
        l2 = Leaf(b"TGCA", h)
        rope = rope_concat(l1, l2, h)
        assert region_hash(rope, 2, 4, h) == h.hash(b"GTTG")


class TestMutationLocalization:
    def test_single_snp(self, h):
        ref = b"ACGTACGTACGTACGT"
        mut = bytearray(ref)
        mut[10] = ord("T")
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
        ref = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
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
        ref_chunks = [b"ACGTACGT", b"TGCATGCA", b"AAAACCCC"]
        mut_chunks = [b"ACGTACGT", b"TGCATGCA", b"AAAATCCC"]

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

    def test_mutation_on_fasta_loaded_rope(self, h, sample_fasta):
        """End-to-end: load FASTA, introduce mutation, localize it."""
        fasta_path, ref_bytes = sample_fasta
        ref_rope, _, _ = load_fasta_to_rope(fasta_path, chunk_size=256, h=h)

        # Mutate position 500
        mut_bytes = bytearray(ref_bytes)
        mut_bytes[500] = ord("N")
        mut_bytes = bytes(mut_bytes)

        # Build mutant rope from bytes directly (simulating a second FASTA)
        mut_rope = None
        for i in range(0, len(mut_bytes), 256):
            chunk = mut_bytes[i:i + 256]
            leaf = Leaf(chunk, h)
            mut_rope = rope_concat(mut_rope, leaf, h) if mut_rope else leaf

        pos = localize_mutation(ref_rope, mut_rope, 0, 1600, h)
        assert pos == 500


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
