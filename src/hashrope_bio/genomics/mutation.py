"""Exact mutation localization using byte comparisons and shared structure.

Equal regions may require a full byte scan; fingerprints never prune changes.
"""

from __future__ import annotations

from hashrope import PolynomialHash, rope_range_equal, rope_range_lcp


def localize_mutation(
    ref_rope,
    sample_rope,
    start: int,
    length: int,
    h: PolynomialHash,
    min_block: int = 1,
) -> int | None:
    """Return the first divergent position, or its binary-search block start.

    ``min_block=1`` gives an exact byte position. Larger values preserve the
    prior binary partition boundaries. No differing-hash precondition applies.
    """
    if min_block < 1:
        raise ValueError("min_block must be at least 1")
    common = rope_range_lcp(ref_rope, start, sample_rope, start, length)
    if common == length:
        return None
    position = start + common
    lo, hi = start, start + length
    while hi - lo > min_block:
        mid = (lo + hi) // 2
        if position < mid:
            hi = mid
        else:
            lo = mid
    return lo


def localize_all_mutations(
    ref_rope,
    sample_rope,
    start: int,
    length: int,
    h: PolynomialHash,
    min_block: int = 1,
) -> list[tuple[int, int]]:
    """Recursively find all divergent regions between ref and sample.

    Returns a list of (start, length) tuples for each divergent block
    at min_block resolution.

    Args:
        ref_rope: Reference sequence rope.
        sample_rope: Sample sequence rope.
        start: Start of the search region (0-based).
        length: Size of the search region in bytes.
        h: Shared PolynomialHash instance.
        min_block: Resolution of the search (smallest reported block size).

    Returns:
        List of (start, length) tuples identifying divergent regions.
    """
    if min_block < 1:
        raise ValueError("min_block must be at least 1")
    results: list[tuple[int, int]] = []
    _find_divergent(ref_rope, sample_rope, start, length, h, min_block, results)
    return results


def _find_divergent(
    ref_rope,
    sample_rope,
    start: int,
    length: int,
    h: PolynomialHash,
    min_block: int,
    results: list[tuple[int, int]],
) -> None:
    """Recursive helper for localize_all_mutations."""
    if rope_range_equal(ref_rope, start, sample_rope, start, length):
        return  # Verified equal, including shared immutable ranges.

    if length <= min_block:
        results.append((start, length))
        return

    mid = length // 2
    _find_divergent(ref_rope, sample_rope, start, mid, h, min_block, results)
    _find_divergent(ref_rope, sample_rope, start + mid, length - mid, h, min_block, results)
