//! E-D2: Drug Resistance Mutation Panel — Rust Benchmark
//!
//! Mirrors the Python bench_resistance.py for direct cross-language comparison.
//! Tests: single-leaf rope, chunked ropes (64/256), byte-slice baseline.

use criterion::{black_box, criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use hashrope::{Arena, Node};
use hashrope_bio::gene::{make_synthetic_gene, mutate_codon};
use hashrope_bio::resistance::*;

/// Build a single-leaf rope from gene bytes.
fn build_single_leaf(arena: &mut Arena, gene: &[u8]) -> Node {
    arena.from_bytes(gene)
}

/// Build a chunked rope from gene bytes.
fn build_chunked(arena: &mut Arena, gene: &[u8], chunk_size: usize) -> Node {
    let mut rope: Node = None;
    for chunk in gene.chunks(chunk_size) {
        let leaf = arena.from_bytes(chunk);
        rope = arena.concat(rope, leaf);
    }
    rope
}

fn bench_resistance_panel(c: &mut Criterion) {
    let panel: Vec<ResistanceSite> = HIV_RT_NNRTI_PANEL.iter()
        .chain(HIV_RT_NRTI_PANEL.iter())
        .cloned()
        .collect();

    let ref_gene = make_synthetic_gene(1680);
    let mut sample_gene = ref_gene.clone();
    mutate_codon(&mut sample_gene, 103);
    mutate_codon(&mut sample_gene, 184);

    let mut group = c.benchmark_group("resistance_panel");

    // --- Single-leaf rope ---
    group.bench_function("hashrope_single_leaf", |b| {
        b.iter_batched(
            || {
                let mut arena = Arena::new();
                let ref_rope = build_single_leaf(&mut arena, &ref_gene);
                let sample_rope = build_single_leaf(&mut arena, &sample_gene);
                (arena, ref_rope, sample_rope)
            },
            |(mut arena, ref_rope, sample_rope)| {
                let results = check_panel_single_arena(
                    &mut arena, ref_rope, sample_rope, black_box(&panel),
                );
                black_box(&results);
                arena // return so drop is outside timing
            },
            BatchSize::SmallInput,
        );
    });

    // --- Chunked ropes ---
    for chunk_size in [64, 256] {
        group.bench_with_input(
            BenchmarkId::new("hashrope_chunked", chunk_size),
            &chunk_size,
            |b, &cs| {
                b.iter_batched(
                    || {
                        let mut arena = Arena::new();
                        let ref_rope = build_chunked(&mut arena, &ref_gene, cs);
                        let sample_rope = build_chunked(&mut arena, &sample_gene, cs);
                        (arena, ref_rope, sample_rope)
                    },
                    |(mut arena, ref_rope, sample_rope)| {
                        let results = check_panel_single_arena(
                            &mut arena, ref_rope, sample_rope, black_box(&panel),
                        );
                        black_box(&results);
                        arena
                    },
                    BatchSize::SmallInput,
                );
            },
        );
    }

    // --- Byte-slice baseline ---
    group.bench_function("byte_slice_baseline", |b| {
        b.iter(|| {
            let results = check_panel_byte_baseline(
                black_box(&ref_gene),
                black_box(&sample_gene),
                black_box(&panel),
            );
            black_box(&results);
        });
    });

    group.finish();
}

fn bench_per_site_substr_hash(c: &mut Criterion) {
    let ref_gene = make_synthetic_gene(1680);

    let mut group = c.benchmark_group("per_site_substr_hash");

    // Single-leaf: how fast is one substr_hash(offset, 3)?
    group.bench_function("single_leaf_3byte", |b| {
        b.iter_batched(
            || {
                let mut arena = Arena::new();
                let rope = build_single_leaf(&mut arena, &ref_gene);
                (arena, rope)
            },
            |(mut arena, rope)| {
                // Check position 103 (offset 306)
                let h = arena.substr_hash(rope, black_box(306), black_box(3));
                black_box(h);
                arena
            },
            BatchSize::SmallInput,
        );
    });

    // Chunked (64-byte): one substr_hash across internal nodes
    group.bench_function("chunked_64_3byte", |b| {
        b.iter_batched(
            || {
                let mut arena = Arena::new();
                let rope = build_chunked(&mut arena, &ref_gene, 64);
                (arena, rope)
            },
            |(mut arena, rope)| {
                let h = arena.substr_hash(rope, black_box(306), black_box(3));
                black_box(h);
                arena
            },
            BatchSize::SmallInput,
        );
    });

    // Chunked (256-byte): one substr_hash
    group.bench_function("chunked_256_3byte", |b| {
        b.iter_batched(
            || {
                let mut arena = Arena::new();
                let rope = build_chunked(&mut arena, &ref_gene, 256);
                (arena, rope)
            },
            |(mut arena, rope)| {
                let h = arena.substr_hash(rope, black_box(306), black_box(3));
                black_box(h);
                arena
            },
            BatchSize::SmallInput,
        );
    });

    // Byte-slice baseline: one 3-byte comparison
    group.bench_function("byte_slice_3byte", |b| {
        let sample_gene = ref_gene.clone();
        b.iter(|| {
            let eq = black_box(&ref_gene[306..309]) == black_box(&sample_gene[306..309]);
            black_box(eq);
        });
    });

    group.finish();
}

criterion_group!(benches, bench_resistance_panel, bench_per_site_substr_hash);
criterion_main!(benches);
