//! Drug resistance mutation panel checking via substr_hash.
//!
//! Checks N known resistance sites by comparing codon hashes between
//! a reference and sample rope. Each check is a single `substr_hash`
//! call at O(log w).

use hashrope::Arena;
use hashrope::rope::Node;

/// A resistance mutation site.
#[derive(Debug, Clone)]
pub struct ResistanceSite {
    pub gene: &'static str,
    pub position: usize,       // 1-based amino acid position
    pub wildtype_aa: &'static str,
    pub mutant_aa: &'static str,
    pub drug_class: &'static str,
    pub annotation: &'static str,
}

/// Result of checking one site.
#[derive(Debug)]
pub struct PanelResult {
    pub position: usize,
    pub ref_hash: u64,
    pub sample_hash: u64,
    pub is_mutant: bool,
}

/// Check all sites in a resistance panel.
///
/// For each site, compares the 3-byte codon hash at that position
/// between reference and sample ropes.
pub fn check_resistance_panel(
    arena_ref: &mut Arena,
    rope_ref: Node,
    arena_sample: &mut Arena,
    rope_sample: Node,
    panel: &[ResistanceSite],
) -> Vec<PanelResult> {
    let mut results = Vec::with_capacity(panel.len());
    for site in panel {
        let offset = ((site.position - 1) * 3) as u64;
        let ref_hash = arena_ref.substr_hash(rope_ref, offset, 3);
        let sample_hash = arena_sample.substr_hash(rope_sample, offset, 3);
        results.push(PanelResult {
            position: site.position,
            ref_hash,
            sample_hash,
            is_mutant: ref_hash != sample_hash,
        });
    }
    results
}

/// Check panel using a single shared arena (both ropes in same arena).
pub fn check_panel_single_arena(
    arena: &mut Arena,
    rope_ref: Node,
    rope_sample: Node,
    panel: &[ResistanceSite],
) -> Vec<PanelResult> {
    let mut results = Vec::with_capacity(panel.len());
    for site in panel {
        let offset = ((site.position - 1) * 3) as u64;
        let ref_hash = arena.substr_hash(rope_ref, offset, 3);
        let sample_hash = arena.substr_hash(rope_sample, offset, 3);
        results.push(PanelResult {
            position: site.position,
            ref_hash,
            sample_hash,
            is_mutant: ref_hash != sample_hash,
        });
    }
    results
}

/// Byte-slice baseline: extract codons and compare directly.
pub fn check_panel_byte_baseline(
    ref_bytes: &[u8],
    sample_bytes: &[u8],
    panel: &[ResistanceSite],
) -> Vec<bool> {
    panel.iter().map(|site| {
        let offset = (site.position - 1) * 3;
        ref_bytes[offset..offset + 3] != sample_bytes[offset..offset + 3]
    }).collect()
}

// --- Pre-defined HIV-1 RT panels (matching Python) ---

pub const HIV_RT_NNRTI_PANEL: &[ResistanceSite] = &[
    ResistanceSite { gene: "RT", position: 100, wildtype_aa: "L", mutant_aa: "I", drug_class: "NNRTI", annotation: "L100I" },
    ResistanceSite { gene: "RT", position: 101, wildtype_aa: "K", mutant_aa: "E", drug_class: "NNRTI", annotation: "K101E" },
    ResistanceSite { gene: "RT", position: 103, wildtype_aa: "K", mutant_aa: "N", drug_class: "NNRTI", annotation: "K103N" },
    ResistanceSite { gene: "RT", position: 106, wildtype_aa: "V", mutant_aa: "A", drug_class: "NNRTI", annotation: "V106A" },
    ResistanceSite { gene: "RT", position: 108, wildtype_aa: "V", mutant_aa: "I", drug_class: "NNRTI", annotation: "V108I" },
    ResistanceSite { gene: "RT", position: 138, wildtype_aa: "E", mutant_aa: "K", drug_class: "NNRTI", annotation: "E138K" },
    ResistanceSite { gene: "RT", position: 181, wildtype_aa: "Y", mutant_aa: "C", drug_class: "NNRTI", annotation: "Y181C" },
    ResistanceSite { gene: "RT", position: 188, wildtype_aa: "Y", mutant_aa: "L", drug_class: "NNRTI", annotation: "Y188L" },
    ResistanceSite { gene: "RT", position: 190, wildtype_aa: "G", mutant_aa: "A", drug_class: "NNRTI", annotation: "G190A" },
    ResistanceSite { gene: "RT", position: 225, wildtype_aa: "P", mutant_aa: "H", drug_class: "NNRTI", annotation: "P225H" },
    ResistanceSite { gene: "RT", position: 227, wildtype_aa: "F", mutant_aa: "L", drug_class: "NNRTI", annotation: "F227L" },
    ResistanceSite { gene: "RT", position: 230, wildtype_aa: "M", mutant_aa: "L", drug_class: "NNRTI", annotation: "M230L" },
];

pub const HIV_RT_NRTI_PANEL: &[ResistanceSite] = &[
    ResistanceSite { gene: "RT", position: 41, wildtype_aa: "M", mutant_aa: "L", drug_class: "NRTI", annotation: "M41L" },
    ResistanceSite { gene: "RT", position: 65, wildtype_aa: "K", mutant_aa: "R", drug_class: "NRTI", annotation: "K65R" },
    ResistanceSite { gene: "RT", position: 67, wildtype_aa: "D", mutant_aa: "N", drug_class: "NRTI", annotation: "D67N" },
    ResistanceSite { gene: "RT", position: 70, wildtype_aa: "K", mutant_aa: "R", drug_class: "NRTI", annotation: "K70R" },
    ResistanceSite { gene: "RT", position: 74, wildtype_aa: "L", mutant_aa: "V", drug_class: "NRTI", annotation: "L74V" },
    ResistanceSite { gene: "RT", position: 115, wildtype_aa: "Y", mutant_aa: "F", drug_class: "NRTI", annotation: "Y115F" },
    ResistanceSite { gene: "RT", position: 151, wildtype_aa: "Q", mutant_aa: "M", drug_class: "NRTI", annotation: "Q151M" },
    ResistanceSite { gene: "RT", position: 184, wildtype_aa: "M", mutant_aa: "V", drug_class: "NRTI", annotation: "M184V" },
    ResistanceSite { gene: "RT", position: 210, wildtype_aa: "L", mutant_aa: "W", drug_class: "NRTI", annotation: "L210W" },
    ResistanceSite { gene: "RT", position: 215, wildtype_aa: "T", mutant_aa: "Y", drug_class: "NRTI", annotation: "T215Y" },
    ResistanceSite { gene: "RT", position: 219, wildtype_aa: "K", mutant_aa: "Q", drug_class: "NRTI", annotation: "K219Q" },
];

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gene::{make_synthetic_gene, mutate_codon};

    #[test]
    fn test_panel_detects_mutations() {
        let ref_gene = make_synthetic_gene(1680);
        let mut sample_gene = ref_gene.clone();
        mutate_codon(&mut sample_gene, 103);
        mutate_codon(&mut sample_gene, 184);

        let mut arena = Arena::new();
        let ref_rope = arena.from_bytes(&ref_gene);
        let sample_rope = arena.from_bytes(&sample_gene);

        let panel: Vec<ResistanceSite> = HIV_RT_NNRTI_PANEL.iter()
            .chain(HIV_RT_NRTI_PANEL.iter())
            .cloned()
            .collect();

        let results = check_panel_single_arena(&mut arena, ref_rope, sample_rope, &panel);

        let mutant_positions: Vec<usize> = results.iter()
            .filter(|r| r.is_mutant)
            .map(|r| r.position)
            .collect();

        assert!(mutant_positions.contains(&103));
        assert!(mutant_positions.contains(&184));
        assert_eq!(mutant_positions.len(), 2);
    }

    #[test]
    fn test_no_mutations_detected_on_identical() {
        let gene = make_synthetic_gene(1680);
        let mut arena = Arena::new();
        let ref_rope = arena.from_bytes(&gene);
        let sample_rope = arena.from_bytes(&gene);

        let results = check_panel_single_arena(
            &mut arena, ref_rope, sample_rope, HIV_RT_NNRTI_PANEL,
        );

        assert!(results.iter().all(|r| !r.is_mutant));
    }

    #[test]
    fn test_byte_baseline_agrees_with_hashrope() {
        let ref_gene = make_synthetic_gene(1680);
        let mut sample_gene = ref_gene.clone();
        mutate_codon(&mut sample_gene, 103);
        mutate_codon(&mut sample_gene, 184);

        let mut arena = Arena::new();
        let ref_rope = arena.from_bytes(&ref_gene);
        let sample_rope = arena.from_bytes(&sample_gene);

        let panel: Vec<ResistanceSite> = HIV_RT_NNRTI_PANEL.iter()
            .chain(HIV_RT_NRTI_PANEL.iter())
            .cloned()
            .collect();

        let hashrope_results = check_panel_single_arena(&mut arena, ref_rope, sample_rope, &panel);
        let baseline_results = check_panel_byte_baseline(&ref_gene, &sample_gene, &panel);

        for (hr, bl) in hashrope_results.iter().zip(baseline_results.iter()) {
            assert_eq!(hr.is_mutant, *bl, "Disagree at position {}", hr.position);
        }
    }
}
