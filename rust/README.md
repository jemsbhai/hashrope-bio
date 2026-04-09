# hashrope-bio

Computational biology applications of [hashrope](https://crates.io/crates/hashrope) —
O(log w) substring hashing for gene-level change detection, drug resistance screening,
genomic region comparison, and tandem repeat analysis.

## Features

- **Gene-level change detection**: Screen hundreds of genes between two genomes in O(G · log w) total time
- **Drug resistance panels**: Check N known mutation sites in O(N · log w) via codon hash comparison
- **Pre-defined HIV-1 RT panels**: 23 NNRTI + NRTI resistance sites from Stanford HIVDB
- **Gene annotation loading**: Parse UCSC refGene TSV for biologically meaningful queries

## Quick start

### Gene-level change detection

```rust
use hashrope::Arena;
use hashrope_bio::gene_diff::*;

// Load gene annotations
let genes = load_gene_regions("chr22_genes.tsv").unwrap();

// Build ropes for reference and sample genomes
let mut arena = Arena::new();
let ref_rope = arena.from_bytes(&ref_seq);
let sample_rope = arena.from_bytes(&sample_seq);

// Screen all genes — O(G · log w) total
let changed = screen_genes_single_arena(
    &mut arena, ref_rope, sample_rope, &genes, None,
);
println!("{} genes differ", changed.len());

// Detailed diff with exon-level resolution
let report = diff_genes_single_arena(
    &mut arena, ref_rope, sample_rope, &genes, None, true,
);
println!("{}", report.summary());
// Gene diff: 12/643 genes changed, 47/5127 exons changed
//   LARGE1 (650,315 bp): body changed, 3/14 exons differ
```

### Drug resistance panels

```rust
use hashrope::Arena;
use hashrope_bio::resistance::*;
use hashrope_bio::gene::{make_synthetic_gene, mutate_codon};

let ref_gene = make_synthetic_gene(1680);
let mut sample = ref_gene.clone();
mutate_codon(&mut sample, 103); // K103N resistance mutation

let mut arena = Arena::new();
let ref_rope = arena.from_bytes(&ref_gene);
let sample_rope = arena.from_bytes(&sample);

let panel: Vec<_> = HIV_RT_NNRTI_PANEL.iter()
    .chain(HIV_RT_NRTI_PANEL.iter())
    .cloned().collect();

let results = check_panel_single_arena(&mut arena, ref_rope, sample_rope, &panel);
let mutant_sites: Vec<_> = results.iter()
    .filter(|r| r.is_mutant)
    .map(|r| r.position)
    .collect();
assert!(mutant_sites.contains(&103));
```

## Validated results

All benchmarks use real biological data, not synthetic:

| Experiment | Data | Key result |
|:-----------|:-----|:-----------|
| Gene change detection | 643 chr22 genes (UCSC refGene) | **5,770/5,770 correct**, 25–411× speedup |
| Resistance screening | 100 clinical HIV-1 sequences vs Stanford HIVDB | **103/103 DRMs detected (100% sensitivity)** |
| Mutation localization | 3,422 ClinVar pathogenic SNVs on chr22 | **3,422/3,422 correct (100%)** |
| Region queries | chr22 at L=1 Mbp, chunk=256 | **936× speedup** vs linear |

## Benchmarks

Criterion benchmarks included:

```bash
cargo bench
```

## License

MIT
