"""Exact lookup of serialized canonical compounds using verified hash buckets.

Fingerprints choose candidate buckets; retained bytes establish identity. This
compares canonical encodings, not chemical equivalence or structural similarity.
"""
from __future__ import annotations

from pathlib import Path
from hashrope import PolynomialHash


class CompoundIndex(dict[int, list[int]]):
    """Candidate buckets plus owned bytes needed for exact verification."""

    def __init__(self, h: PolynomialHash) -> None:
        super().__init__()
        self.profile = (h.prime, h.base)
        self.records: dict[int, bytes] = {}


def build_compound_index(
    smiles_file: str | Path,
    h: PolynomialHash | None = None,
) -> CompoundIndex:
    """Index one canonical encoding per line, optionally followed by a tab/ID.

    All candidates and duplicate records are retained. Returned IDs are zero-
    based physical line indices; empty lines are skipped. Input canonicalization
    remains the caller's responsibility. No text decoding changes the bytes.
    """
    h = h if h is not None else PolynomialHash()
    index = CompoundIndex(h)
    with open(smiles_file, "rb") as stream:
        for line_number, line in enumerate(stream):
            encoding = line.rstrip(b"\r\n").split(b"\t", 1)[0]
            if not encoding:
                continue
            index.records[line_number] = encoding
            index.setdefault(h.hash(encoding), []).append(line_number)
    return index


def lookup_compound(
    query_smiles: bytes,
    index: CompoundIndex,
    h: PolynomialHash,
) -> list[int]:
    """Return only candidates whose retained canonical encoding equals query.

    Legacy hash-only dictionaries cannot establish identity and are rejected;
    rebuild them with ``build_compound_index`` to retain verification bytes.
    """
    if not isinstance(index, CompoundIndex):
        raise ValueError("Exact lookup requires retained bytes; rebuild with build_compound_index")
    if index.profile != (h.prime, h.base):
        raise ValueError("Compound index and query hash profiles differ")
    query = bytes(query_smiles)
    return [row for row in index.get(h.hash(query), []) if index.records[row] == query]


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
