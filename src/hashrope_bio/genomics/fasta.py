"""Load FASTA sequences into hash ropes.

Handles both indexed (pysam) and unindexed (biopython) FASTA files.
Chunks sequences into fixed-size blocks for rope construction.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterator

from hashrope import PolynomialHash, Leaf, rope_concat


def load_fasta_to_rope(
    path: str | Path,
    region: str | None = None,
    chunk_size: int = 4096,
) -> tuple:
    """Load a FASTA file (or region) into a hash rope.

    Args:
        path: Path to .fa/.fasta file. Must have .fai index for pysam.
        region: Optional samtools-style region string, e.g. "chr22:1000-2000".
                If None, loads the entire first sequence.
        chunk_size: Bytes per leaf node. Affects tree depth and query speed.
                    Recommended: 4096 for genomes, 256 for small genes.

    Returns:
        (rope_root, hasher, metadata) where metadata is a dict with keys:
            - seq_name: str
            - seq_len: int
            - chunk_count: int
            - construction_time_s: float
    """
    raise NotImplementedError("TODO: implement FASTA loading")


def load_fasta_chunked(
    path: str | Path,
    chunk_size: int = 4096,
) -> Iterator[bytes]:
    """Yield fixed-size byte chunks from a FASTA file (headers stripped).

    Last chunk may be shorter than chunk_size.
    """
    raise NotImplementedError("TODO: implement chunked FASTA reader")


def _load_with_pysam(path: Path, region: str | None) -> bytes:
    """Extract sequence bytes via pysam (requires .fai index)."""
    raise NotImplementedError


def _load_with_biopython(path: Path) -> bytes:
    """Extract sequence bytes via biopython (no index needed, slower)."""
    raise NotImplementedError
