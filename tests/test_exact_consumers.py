"""Regression checks for decisions that previously trusted matching hashes."""
import struct

import pytest
from hashrope import PolynomialHash, Leaf, rope_concat, rope_hash
from hashrope_bio.genomics.region_query import regions_identical
from hashrope_bio.genomics.mutation import localize_mutation, localize_all_mutations
from hashrope_bio.genomics.gene_diff import Exon, GeneRegion, diff_genes, screen_genes
from hashrope_bio.cheminformatics.resistance import ResistanceSite, check_resistance_panel
from hashrope_bio.cheminformatics.compounds import build_compound_index, lookup_compound
from hashrope_bio.proteomics.trajectory import frames_identical


def test_colliding_regions_mutations_genes_exons_and_codons():
    h = PolynomialHash()
    a, b = Leaf(b"\x00\x83*", h), Leaf(b"\x01\x00*", h)
    assert rope_hash(a) == rope_hash(b)
    assert not regions_identical(a, b, 0, 3, h)
    assert localize_mutation(a, b, 0, 3, h) == 0
    assert localize_all_mutations(a, b, 0, 3, h) == [(0, 1), (1, 1)]
    gene = GeneRegion("test", "test", "chr1", "+", 0, 3, [Exon(1, 0, 3)])
    report = diff_genes(a, b, [gene], h)
    result = report.gene_diffs[0]
    assert result.body_ref_hash == result.body_sample_hash
    assert result.body_changed
    assert result.changed_exons[0].changed
    assert screen_genes(a, b, [gene], h) == ["test"]
    site = ResistanceSite("test", 1, "X", "Y", "test", "byte regression")
    panel = check_resistance_panel(a, b, [site], h)
    assert panel[0].ref_hash == panel[0].sample_hash
    assert panel[0].is_mutant


def test_colliding_finite_float_frames_are_not_identical():
    # Distinct finite float32 coordinates with the base-131 carry collision in
    # the first two bytes; all other coordinates are identical.
    a = bytes.fromhex("0083003f") + struct.pack("<ff", 1.0, 2.0)
    b = bytes.fromhex("0100003f") + struct.pack("<ff", 1.0, 2.0)
    assert struct.unpack("<fff", a) != struct.unpack("<fff", b)
    h = PolynomialHash()
    assert h.hash(a) == h.hash(b)
    rope = rope_concat(Leaf(a, h), Leaf(b, h), h)
    assert not frames_identical(rope, 0, 1, 12, h)
    assert frames_identical(rope, 0, 0, 12, h)


def test_compound_buckets_preserve_collisions_and_duplicate_records(tmp_path):
    h = PolynomialHash(base=2)
    assert h.hash(b"AC") == h.hash(b"BA")
    source = tmp_path / "compounds.tsv"
    source.write_bytes(b"AC\tfirst\nBA\tsecond\nAC\tduplicate\n")
    index = build_compound_index(source, h)
    assert index[h.hash(b"AC")] == [0, 1, 2]
    assert lookup_compound(b"AC", index, h) == [0, 2]
    assert lookup_compound(b"BA", index, h) == [1]
    assert lookup_compound(b"CC", index, h) == []
    with pytest.raises(ValueError, match="retained bytes"):
        lookup_compound(b"AC", dict(index), h)
    with pytest.raises(ValueError, match="profiles differ"):
        lookup_compound(b"AC", index, PolynomialHash())


def test_invalid_regions_and_mutation_resolution_are_rejected():
    h = PolynomialHash()
    a = Leaf(b"ACG", h)
    with pytest.raises(ValueError):
        regions_identical(a, a, 1, 3, h)
    with pytest.raises(ValueError):
        frames_identical(a, 0, 1, 3, h)
    for fn in [localize_mutation, localize_all_mutations]:
        with pytest.raises(ValueError, match="min_block"):
            fn(a, a, 0, 3, h, min_block=0)


def test_existing_benchmark_consumers_also_verify_colliding_dna(monkeypatch):
    from pathlib import Path
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1] / "benchmarks"))
    import bench_cross_species
    import bench_resistance_realdata

    h = PolynomialHash(base=2)
    assert h.hash(b"AGAA") == h.hash(b"CCAA")
    compared = bench_cross_species.compare_gene("AGAA", "CCAA", h)
    assert not compared["is_conserved"]
    assert compared["first_divergence_pos"] == 0
    monkeypatch.setattr(bench_resistance_realdata, "PANEL_POSITIONS", [1])
    panel = bench_resistance_realdata.validate_sequence(b"AGA", "CCA", (1, 1), h)
    assert panel[1]["hashrope_mismatch"]
    assert panel[1]["hashrope_agrees_with_bytes"]
