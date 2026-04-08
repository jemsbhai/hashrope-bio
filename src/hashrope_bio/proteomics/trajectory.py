"""MD trajectory loading and frame-level comparison via hash ropes.

Serializes atomic coordinates (float32 × 3 × N_atoms) per frame
into byte arrays, builds a rope of frames, and enables O(log w)
frame identity testing via substr_hash.

Requires: MDAnalysis, numpy
"""

from __future__ import annotations

from pathlib import Path

from hashrope import PolynomialHash, rope_substr_hash


def load_trajectory_to_rope(
    topology: str | Path,
    trajectory: str | Path,
    chunk_frames: int = 1,
    selection: str = "all",
) -> tuple:
    """Load an MD trajectory into a hash rope.

    Each frame's atomic positions are serialized to bytes (float32, little-endian)
    and stored as a leaf. Frames can be grouped into chunks for efficiency.

    Args:
        topology: Path to topology file (.pdb, .gro, .psf).
        trajectory: Path to trajectory file (.xtc, .trr, .dcd).
        chunk_frames: Number of frames per leaf node.
        selection: MDAnalysis atom selection string.

    Returns:
        (rope_root, hasher, metadata) where metadata contains:
            - n_frames: int
            - n_atoms: int
            - frame_size_bytes: int (per frame)
            - construction_time_s: float
    """
    raise NotImplementedError("TODO: implement trajectory loading")


def frames_identical(
    rope,
    frame_i: int,
    frame_j: int,
    frame_size: int,
    h: PolynomialHash,
) -> bool:
    """Test bitwise identity of two frames in O(log w).

    Args:
        rope: Pre-built trajectory rope.
        frame_i: Index of first frame (0-based).
        frame_j: Index of second frame (0-based).
        frame_size: Byte size of one frame.
        h: PolynomialHash instance.

    Returns:
        True if frames are bitwise identical.

    Note:
        This tests bitwise identity of serialized coordinates, not
        structural similarity (RMSD). Floating-point noise means
        "same conformation" may not produce identical bytes.
    """
    h_i = rope_substr_hash(rope, frame_i * frame_size, frame_size, h)
    h_j = rope_substr_hash(rope, frame_j * frame_size, frame_size, h)
    return h_i == h_j
