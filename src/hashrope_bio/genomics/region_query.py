"""Region fingerprints and exact byte identity queries on genome-scale ropes."""

from __future__ import annotations

from hashrope import PolynomialHash, rope_substr_hash, rope_range_equal


def region_hash(rope, start: int, length: int, h: PolynomialHash) -> int:
    """Compute H(sequence[start:start+length]) in O(log w).

    Args:
        rope: Pre-built hash rope root node.
        start: 0-based byte offset into the sequence.
        length: Number of bytes to hash.
        h: PolynomialHash instance used to build the rope.

    Returns:
        Polynomial hash of the subsequence.
    """
    return rope_substr_hash(rope, start, length, h)


def regions_identical(
    rope_a,
    rope_b,
    start: int,
    length: int,
    h: PolynomialHash,
) -> bool:
    """Compare the actual bytes in [start, start+length).

    Shared structure can be skipped; independently stored equal regions may
    require reading every byte. ``h`` is retained for API compatibility and
    does not determine equality. Invalid ranges raise ValueError.
    """
    return rope_range_equal(rope_a, start, rope_b, start, length)
