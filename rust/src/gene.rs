//! Synthetic gene generation and sequence utilities.

/// Generate a deterministic pseudo-random gene sequence.
///
/// Uses the same LCG as the Python benchmark for cross-validation.
/// Seed = 42, matching `make_synthetic_gene()` in bench_resistance.py.
pub fn make_synthetic_gene(length: usize) -> Vec<u8> {
    let nucs = b"ACGT";
    let mut state: u64 = 42;
    let mut gene = Vec::with_capacity(length);
    for _ in 0..length {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        gene.push(nucs[((state >> 33) & 3) as usize]);
    }
    gene
}

/// Introduce a mutation at a codon position (1-based amino acid position).
///
/// Flips the first nucleotide of the codon: if it's G, change to A; otherwise change to G.
/// This matches the Python benchmark's mutation strategy.
pub fn mutate_codon(gene: &mut [u8], aa_position: usize) {
    let offset = (aa_position - 1) * 3;
    if gene[offset] == b'G' {
        gene[offset] = b'A';
    } else {
        gene[offset] = b'G';
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_synthetic_gene_deterministic() {
        let g1 = make_synthetic_gene(100);
        let g2 = make_synthetic_gene(100);
        assert_eq!(g1, g2);
        assert_eq!(g1.len(), 100);
        assert!(g1.iter().all(|&b| b == b'A' || b == b'C' || b == b'G' || b == b'T'));
    }

    #[test]
    fn test_mutate_codon() {
        let mut gene = b"AAACCCGGGTTT".to_vec();
        mutate_codon(&mut gene, 1);
        assert_eq!(gene[0], b'G');
        mutate_codon(&mut gene, 3);
        assert_eq!(gene[6], b'A');
    }
}
