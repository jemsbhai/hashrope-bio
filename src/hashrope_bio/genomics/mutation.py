"""Binary search mutation localization via substr_hash.

Given a reference rope and a sample rope known to differ,
localize the mutation(s) in O(log N) hash comparisons where
N is the region size — each comparison is O(log w).

Total: O(log N · log w) for single mutation in N-bp region.
"""

from __future__ import annotations

from hashrope import PolynomialHash, rope_substr_hash


def localize_mutation(
    ref_rope,
    sample_rope,
    start: int,
    length: int,
    h: PolynomialHash,
    min_block: int = 1,
) -> int | None:
    """Binary search for a single mutation position.

    Precondition: the region [start, start+length) is known to differ
    between ref and sample (i.e., their substr_hashes are unequal).

    Args:
        ref_rope: Reference sequence rope.
        sample_rope: Sample sequence rope (contains mutation).
        start: Start of the search region (0-based).
        length: Size of the search region in bytes.
        h: Shared PolynomialHash instance.
        min_block: Stop when region is this small (default: 1 = exact position).

    Returns:
        Position of the mutation (0-based), or None if regions are identical.
    """
    lo, hi = start, start + length

    # Confirm the region actually differs
    h_ref = rope_substr_hash(ref_rope, lo, hi - lo, h)
    h_sam = rope_substr_hash(sample_rope, lo, hi - lo, h)
    if h_ref == h_sam:
        return None

    while hi - lo > min_block:
        mid = (lo + hi) // 2
        h_ref_left = rope_substr_hash(ref_rope, lo, mid - lo, h)
        h_sam_left = rope_substr_hash(sample_rope, lo, mid - lo, h)
        if h_ref_left != h_sam_left:
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
    h_ref = rope_substr_hash(ref_rope, start, length, h)
    h_sam = rope_substr_hash(sample_rope, start, length, h)

    if h_ref == h_sam:
        return  # Identical — prune entire subtree

    if length <= min_block:
        results.append((start, length))
        return

    mid = length // 2
    _find_divergent(ref_rope, sample_rope, start, mid, h, min_block, results)
    _find_divergent(ref_rope, sample_rope, start + mid, length - mid, h, min_block, results)
