"""Tandem repeat handling via RepeatNode.

RepeatNode represents "motif × q" in O(1) space (228 bytes) regardless
of q, with hash computed in O(log q). Relevant for STR/VNTR analysis
and repeat expansion disorders (Huntington's, Fragile X, etc.).
"""

from __future__ import annotations

from hashrope import PolynomialHash, Leaf, rope_repeat, rope_hash, rope_len


def build_repeat_node(
    motif: bytes,
    count: int,
    h: PolynomialHash,
):
    """Create a RepeatNode for motif × count.

    Args:
        motif: The repeat unit (e.g., b"CAG" for Huntington's).
        count: Number of repeats.
        h: PolynomialHash instance.

    Returns:
        Rope node representing the repeated sequence.
        Hash is computed in O(log count), memory is O(1).
    """
    leaf = Leaf(motif, h)
    return rope_repeat(leaf, count, h)


def detect_tandem_repeat(
    sequence: bytes,
    motif: bytes,
) -> int | None:
    """Count exact tandem repeats of motif at the start of sequence.

    Simple scan — this is the baseline, not the hashrope-accelerated version.

    Args:
        sequence: Byte sequence to scan.
        motif: Repeat unit to search for.

    Returns:
        Number of complete repeats found, or None if motif not found at start.
    """
    d = len(motif)
    if d == 0 or len(sequence) < d:
        return None
    if sequence[:d] != motif:
        return None

    count = 0
    for i in range(0, len(sequence) - d + 1, d):
        if sequence[i:i + d] == motif:
            count += 1
        else:
            break
    return count


# --- Clinically relevant repeat loci ---

CLINICAL_REPEATS = {
    "HTT": {"gene": "Huntingtin", "motif": b"CAG", "normal": (10, 35), "pathogenic": (36, None),
             "disease": "Huntington's disease"},
    "FMR1": {"gene": "FMRP", "motif": b"CGG", "normal": (5, 44), "pathogenic": (200, None),
              "disease": "Fragile X syndrome"},
    "ATXN1": {"gene": "Ataxin-1", "motif": b"CAG", "normal": (6, 35), "pathogenic": (39, None),
               "disease": "Spinocerebellar ataxia type 1"},
    "ATXN3": {"gene": "Ataxin-3", "motif": b"CAG", "normal": (12, 40), "pathogenic": (55, None),
               "disease": "Machado-Joseph disease"},
    "DMPK": {"gene": "DMPK", "motif": b"CTG", "normal": (5, 34), "pathogenic": (50, None),
              "disease": "Myotonic dystrophy type 1"},
}
