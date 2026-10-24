use hashrope::{Arena, MERSENNE_61};
use hashrope_bio::gene_diff::{diff_genes, diff_genes_single_arena, screen_genes_single_arena, Exon, GeneRegion};
use hashrope_bio::resistance::{check_panel_single_arena, check_resistance_panel, ResistanceSite};

fn gene() -> GeneRegion {
    GeneRegion { name: "test".into(), transcript: "test".into(), chrom: "chr1".into(),
        strand: '+', tx_start: 0, tx_end: 3,
        exons: vec![Exon { index: 1, start: 0, end: 3 }] }
}

fn site() -> ResistanceSite {
    ResistanceSite { gene: "test", position: 1, wildtype_aa: "X", mutant_aa: "Y",
        drug_class: "test", annotation: "byte regression" }
}

#[test]
fn matching_hashes_do_not_hide_gene_exon_or_codon_changes() {
    let mut arena = Arena::new();
    let a = arena.from_bytes(&[0, 131, 42]);
    let b = arena.from_bytes(&[1, 0, 42]);
    assert_eq!(arena.hash(a), arena.hash(b));
    let report = diff_genes_single_arena(&mut arena, a, b, &[gene()], None, true);
    assert!(report.gene_diffs[0].body_changed());
    assert!(report.gene_diffs[0].exon_diffs[0].changed());
    assert_eq!(screen_genes_single_arena(&mut arena, a, b, &[gene()], None), vec!["test"]);
    let panel = check_panel_single_arena(&mut arena, a, b, &[site()]);
    assert_eq!(panel[0].ref_hash, panel[0].sample_hash);
    assert!(panel[0].is_mutant);
}

#[test]
fn different_profiles_do_not_create_spurious_mutations() {
    let mut a = Arena::new_lazy();
    let mut b = Arena::with_hash(MERSENNE_61, 257);
    let left = a.from_bytes(b"ACG");
    let right = b.from_bytes(b"ACG");
    let report = diff_genes(&mut a, left, &mut b, right, &[gene()], None);
    assert_ne!(report.gene_diffs[0].body_ref_hash, report.gene_diffs[0].body_sample_hash);
    assert!(!report.gene_diffs[0].body_changed());
    assert!(!report.gene_diffs[0].exon_diffs[0].changed());
    assert!(!check_resistance_panel(&mut a, left, &mut b, right, &[site()])[0].is_mutant);
}
