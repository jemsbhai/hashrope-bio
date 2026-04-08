"""O(log w) region identity queries on genome-scale ropes.

Core operation: substr_hash on a pre-built rope answers
"is region [start, start+L) identical between two ropes?"
without reading the underlying bytes.
"""

from __future__ import annotations

from hashrope import PolynomialHash, rope_substr_hash


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
    """Test whether two ropes have identical content at [start, start+length).

    O(log w) per rope — total O(log w) since both ropes have similar depth.
    No bytes are read; this is pure hash comparison.

    Args:
        rope_a: Reference rope.
        rope_b: Sample rope.
        start: 0-based byte offset.
        length: Region size in bytes.
        h: Shared PolynomialHash instance.

    Returns:
        True if H(A[start:start+length]) == H(B[start:start+length]).

    Note:
        This is exact identity, not similarity. A single base difference
        causes a mismatch (no false negatives). False positives are
        theoretically possible but astronomically unlikely with Mersenne-61
        (probability ≈ 1/2^61 per comparison).
    """
    h_a = rope_substr_hash(rope_a, start, length, h)
    h_b = rope_substr_hash(rope_b, start, length, h)
    return h_a == h_b
