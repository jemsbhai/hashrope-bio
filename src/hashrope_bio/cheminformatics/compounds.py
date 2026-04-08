"""Chemical compound exact lookup via polynomial hashing.

Builds a hash index over canonical SMILES/InChI strings for O(1)
exact compound identity testing. Not a similarity tool — use RDKit
fingerprints + Tanimoto for similarity search.

Requires: rdkit (optional, for canonicalization)
"""

from __future__ import annotations

from pathlib import Path

from hashrope import PolynomialHash


def build_compound_index(
    smiles_file: str | Path,
    h: PolynomialHash | None = None,
) -> dict[int, list[int]]:
    """Build a hash-based index from a SMILES file.

    Args:
        smiles_file: Path to a file with one SMILES per line
                     (optionally tab-separated: SMILES\\tID).
        h: PolynomialHash instance. Uses default if None.

    Returns:
        Dict mapping polynomial_hash(SMILES) -> list of line indices.
        Multiple entries per hash indicates a collision (expected: ~0 for <10M compounds).
    """
    raise NotImplementedError("TODO: implement compound index builder")


def lookup_compound(
    query_smiles: bytes,
    index: dict[int, list[int]],
    h: PolynomialHash,
) -> list[int]:
    """Look up a compound by exact SMILES hash match.

    Args:
        query_smiles: Canonical SMILES as bytes.
        index: Pre-built hash index from build_compound_index.
        h: Same PolynomialHash used to build the index.

    Returns:
        List of matching line indices (usually 0 or 1 entries).
    """
    query_hash = h.hash(query_smiles)
    return index.get(query_hash, [])


def canonicalize_smiles(smiles: str) -> str | None:
    """Canonicalize a SMILES string via RDKit.

    Returns None if SMILES is invalid.
    """
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol, canonical=True)
    except ImportError:
        # No RDKit — return as-is (caller responsible for canonical form)
        return smiles
