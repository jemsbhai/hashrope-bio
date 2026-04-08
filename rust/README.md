# hashrope-bio

Computational biology applications of [hashrope](https://crates.io/crates/hashrope) —
O(log w) substring hashing for drug resistance screening, genomic region comparison,
and tandem repeat analysis.

## Features

- **Drug resistance panels**: Check N known mutation sites in O(N · log w) total time via `substr_hash`
- **Synthetic gene generation**: Deterministic pseudo-random sequences for reproducible benchmarks
- **Pre-defined HIV-1 RT panels**: 23 NNRTI + NRTI resistance sites from Stanford HIVDB

## Quick start

```rust
use hashrope::Arena;
use hashrope_bio::gene::{make_synthetic_gene, mutate_codon};
use hashrope_bio::resistance::*;

let ref_gene = make_synthetic_gene(1680); // 560-codon gene
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

## Benchmarks

Criterion benchmarks included. Run with:

```bash
cargo bench
```

### Results (Rust, --release, 23-site panel on 1680 bp gene)

| Scenario | Total (ns) | Per site (ns) |
|---|---|---|
| hashrope single-leaf | 2,160 | 94 |
| hashrope chunked/256 | 1,917 | 83 |
| byte-slice baseline | 38 | 1.7 |

Rust hashrope is **17.6× faster** than the equivalent Python hashrope implementation
and matches Python's byte-slice baseline performance.

## License

MIT
