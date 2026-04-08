"""Load FASTA sequences into hash ropes.

Supports three backends (in priority order):
  1. pysam — fast indexed access, requires .fai index
  2. biopython — no index needed, slower
  3. pure Python — zero dependencies, handles standard FASTA

Chunks sequences into fixed-size blocks for rope construction.
"""

from __future__ import annotations

import gzip
import time
from pathlib import Path
from typing import Iterator

from hashrope import PolynomialHash, Leaf, rope_concat


def load_fasta_to_rope(
    path: str | Path,
    region: str | None = None,
    chunk_size: int = 4096,
    h: PolynomialHash | None = None,
    uppercase: bool = True,
) -> tuple:
    """Load a FASTA file (or region) into a hash rope.

    Args:
        path: Path to .fa/.fasta/.fa.gz file.
        region: Optional samtools-style region string, e.g. "chr22:1000-2000".
                If None, loads the entire first sequence.
                Only supported with pysam backend.
        chunk_size: Bytes per leaf node. Affects tree depth and query speed.
                    Recommended: 4096 for genomes, 256 for small genes.
        h: PolynomialHash instance. Uses default (base=131, p=2^61-1) if None.
        uppercase: Convert sequence to uppercase (recommended — FASTA case varies).

    Returns:
        (rope_root, hasher, metadata) where metadata is a dict with keys:
            - seq_name: str
            - seq_len: int
            - chunk_count: int
            - construction_time_s: float
            - backend: str ("pysam" | "biopython" | "pure_python")
    """
    path = Path(path)
    if h is None:
        h = PolynomialHash()

    # Choose backend
    if region is not None:
        # Region queries require pysam
        seq_bytes, seq_name, backend = _load_with_pysam(path, region)
    else:
        # Try pysam first (fastest), fall back through biopython to pure Python
        for loader, name in [(_try_pysam, "pysam"), (_try_biopython, "biopython"), (_try_pure, "pure_python")]:
            result = loader(path)
            if result is not None:
                seq_bytes, seq_name = result
                backend = name
                break
        else:
            raise RuntimeError(f"Could not load FASTA from {path} with any backend")

    if uppercase:
        seq_bytes = seq_bytes.upper()

    # Build rope from chunks
    t0 = time.perf_counter()
    rope_root, chunk_count = _build_rope_from_bytes(seq_bytes, chunk_size, h)
    t1 = time.perf_counter()

    metadata = {
        "seq_name": seq_name,
        "seq_len": len(seq_bytes),
        "chunk_count": chunk_count,
        "construction_time_s": t1 - t0,
        "backend": backend,
    }

    return rope_root, h, metadata


def load_fasta_chunked(
    path: str | Path,
    chunk_size: int = 4096,
    uppercase: bool = True,
) -> Iterator[bytes]:
    """Yield fixed-size byte chunks from a FASTA file (headers stripped).

    Last chunk may be shorter than chunk_size.
    """
    seq_bytes, _ = _load_first_sequence(Path(path))
    if uppercase:
        seq_bytes = seq_bytes.upper()

    for i in range(0, len(seq_bytes), chunk_size):
        yield seq_bytes[i:i + chunk_size]


def load_fasta_bytes(
    path: str | Path,
    region: str | None = None,
    uppercase: bool = True,
) -> tuple[bytes, str]:
    """Load raw sequence bytes from a FASTA file (no rope construction).

    Useful for baselines and verification.

    Returns:
        (sequence_bytes, sequence_name)
    """
    path = Path(path)

    if region is not None:
        seq_bytes, seq_name, _ = _load_with_pysam(path, region)
    else:
        seq_bytes, seq_name = _load_first_sequence(path)

    if uppercase:
        seq_bytes = seq_bytes.upper()

    return seq_bytes, seq_name


# ---------------------------------------------------------------------------
# Rope construction
# ---------------------------------------------------------------------------


def _build_rope_from_bytes(
    data: bytes,
    chunk_size: int,
    h: PolynomialHash,
) -> tuple:
    """Build a rope from contiguous bytes via chunked concatenation.

    Returns (rope_root, chunk_count).
    """
    if len(data) == 0:
        return None, 0

    rope = None
    chunk_count = 0

    for i in range(0, len(data), chunk_size):
        chunk = data[i:i + chunk_size]
        leaf = Leaf(chunk, h)
        rope = rope_concat(rope, leaf, h)
        chunk_count += 1

    return rope, chunk_count


# ---------------------------------------------------------------------------
# Backend: pysam (indexed FASTA)
# ---------------------------------------------------------------------------


def _load_with_pysam(path: Path, region: str) -> tuple[bytes, str, str]:
    """Extract sequence bytes via pysam (requires .fai index).

    Returns (bytes, seq_name, "pysam").
    """
    try:
        import pysam
    except ImportError:
        raise ImportError(
            "pysam is required for region queries. "
            "Install with: pip install pysam"
        )

    fai_path = Path(str(path) + ".fai")
    if not fai_path.exists():
        raise FileNotFoundError(
            f"FASTA index not found: {fai_path}. "
            f"Create it with: samtools faidx {path}"
        )

    with pysam.FastaFile(str(path)) as fasta:
        seq = fasta.fetch(region=region)
        # Parse seq_name from region string
        seq_name = region.split(":")[0] if ":" in region else region
        return seq.encode("ascii"), seq_name, "pysam"


def _try_pysam(path: Path) -> tuple[bytes, str] | None:
    """Try loading first sequence via pysam. Returns None if unavailable."""
    try:
        import pysam
    except ImportError:
        return None

    fai_path = Path(str(path) + ".fai")
    if not fai_path.exists():
        return None

    try:
        with pysam.FastaFile(str(path)) as fasta:
            name = fasta.references[0]
            seq = fasta.fetch(name)
            return seq.encode("ascii"), name
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Backend: biopython
# ---------------------------------------------------------------------------


def _try_biopython(path: Path) -> tuple[bytes, str] | None:
    """Try loading first sequence via biopython. Returns None if unavailable."""
    try:
        from Bio import SeqIO
    except ImportError:
        return None

    try:
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt") as f:
            record = next(SeqIO.parse(f, "fasta"))
            return bytes(record.seq), record.id
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Backend: pure Python (zero dependencies)
# ---------------------------------------------------------------------------


def _try_pure(path: Path) -> tuple[bytes, str] | None:
    """Load first sequence from a FASTA file using pure Python."""
    try:
        return _load_first_sequence(path)
    except Exception:
        return None


def _load_first_sequence(path: Path) -> tuple[bytes, str]:
    """Parse the first sequence from a FASTA file (pure Python).

    Handles .fa, .fasta, .fa.gz, .fasta.gz.
    Strips header lines (starting with '>') and whitespace.
    """
    opener = gzip.open if path.suffix == ".gz" else open
    seq_name = ""
    seq_parts: list[bytes] = []
    found_header = False

    with opener(path, "rt") as f:
        for line in f:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if found_header:
                    # We've hit the second sequence — stop
                    break
                seq_name = line[1:].split()[0]  # First word after '>'
                found_header = True
                continue
            if found_header:
                # Strip whitespace and append
                cleaned = line.replace(" ", "").replace("\t", "")
                if cleaned:
                    seq_parts.append(cleaned.encode("ascii"))

    if not seq_parts:
        raise ValueError(f"No sequence data found in {path}")

    return b"".join(seq_parts), seq_name
