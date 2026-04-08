"""Drug resistance mutation panel checking via substr_hash.

Given a reference gene rope and a sample gene rope, check N known
resistance mutation sites in O(N · log w) total time. Each site
is a codon (3 bp) compared via substr_hash.
"""

from __future__ import annotations

from dataclasses import dataclass

from hashrope import PolynomialHash, rope_substr_hash


@dataclass
class ResistanceSite:
    """A known drug resistance mutation site."""
    gene: str
    position: int          # amino acid position (1-based)
    wildtype_aa: str       # e.g., "K" for lysine
    mutant_aa: str         # e.g., "N" for asparagine
    drug_class: str        # e.g., "NNRTI"
    annotation: str        # e.g., "K103N — efavirenz resistance"


@dataclass
class PanelResult:
    """Result of checking one resistance site."""
    site: ResistanceSite
    ref_hash: int
    sample_hash: int
    is_mutant: bool        # True if hashes differ (potential resistance)


def check_resistance_panel(
    ref_rope,
    sample_rope,
    panel: list[ResistanceSite],
    h: PolynomialHash,
    codon_size: int = 3,
) -> list[PanelResult]:
    """Check all sites in a resistance panel.

    For each site, compares the codon hash at that position between
    reference and sample. Hash mismatch = potential resistance mutation.

    Args:
        ref_rope: Reference gene sequence rope.
        sample_rope: Patient sample gene sequence rope.
        panel: List of ResistanceSite to check.
        h: Shared PolynomialHash instance.
        codon_size: Bytes per codon (default 3 for DNA).

    Returns:
        List of PanelResult, one per site.
    """
    results = []
    for site in panel:
        # Convert 1-based amino acid position to 0-based nucleotide offset
        nt_offset = (site.position - 1) * codon_size

        ref_hash = rope_substr_hash(ref_rope, nt_offset, codon_size, h)
        sample_hash = rope_substr_hash(sample_rope, nt_offset, codon_size, h)

        results.append(PanelResult(
            site=site,
            ref_hash=ref_hash,
            sample_hash=sample_hash,
            is_mutant=(ref_hash != sample_hash),
        ))

    return results


# --- Pre-defined panels ---

# HIV-1 RT NNRTI resistance (Stanford HIVDB, abridged)
HIV_RT_NNRTI_PANEL = [
    ResistanceSite("RT", 100, "L", "I", "NNRTI", "L100I"),
    ResistanceSite("RT", 101, "K", "E", "NNRTI", "K101E"),
    ResistanceSite("RT", 103, "K", "N", "NNRTI", "K103N — primary efavirenz resistance"),
    ResistanceSite("RT", 106, "V", "A", "NNRTI", "V106A"),
    ResistanceSite("RT", 108, "V", "I", "NNRTI", "V108I"),
    ResistanceSite("RT", 138, "E", "K", "NNRTI", "E138K — rilpivirine resistance"),
    ResistanceSite("RT", 181, "Y", "C", "NNRTI", "Y181C — nevirapine resistance"),
    ResistanceSite("RT", 188, "Y", "L", "NNRTI", "Y188L"),
    ResistanceSite("RT", 190, "G", "A", "NNRTI", "G190A"),
    ResistanceSite("RT", 225, "P", "H", "NNRTI", "P225H"),
    ResistanceSite("RT", 227, "F", "L", "NNRTI", "F227L"),
    ResistanceSite("RT", 230, "M", "L", "NNRTI", "M230L"),
]

# HIV-1 RT NRTI resistance (Stanford HIVDB, abridged)
HIV_RT_NRTI_PANEL = [
    ResistanceSite("RT", 41, "M", "L", "NRTI", "M41L — TAM"),
    ResistanceSite("RT", 65, "K", "R", "NRTI", "K65R — tenofovir resistance"),
    ResistanceSite("RT", 67, "D", "N", "NRTI", "D67N — TAM"),
    ResistanceSite("RT", 70, "K", "R", "NRTI", "K70R — TAM"),
    ResistanceSite("RT", 74, "L", "V", "NRTI", "L74V — didanosine resistance"),
    ResistanceSite("RT", 115, "Y", "F", "NRTI", "Y115F"),
    ResistanceSite("RT", 151, "Q", "M", "NRTI", "Q151M — multi-NRTI resistance"),
    ResistanceSite("RT", 184, "M", "V", "NRTI", "M184V — lamivudine/emtricitabine resistance"),
    ResistanceSite("RT", 210, "L", "W", "NRTI", "L210W — TAM"),
    ResistanceSite("RT", 215, "T", "Y", "NRTI", "T215Y — TAM"),
    ResistanceSite("RT", 219, "K", "Q", "NRTI", "K219Q — TAM"),
]
