# PyO3 Bindings — Deferred

**Status**: Deferred — to be implemented after pure-Rust benchmarks establish baseline.

## Summary

The Python E-D2 resistance panel benchmark showed hashrope **18× slower** than Python byte-slice baseline on a 1,680 bp gene (38 µs vs 2.1 µs). Analysis shows this is Python interpreter overhead, not algorithmic — the pure Python `rope_substr_hash` costs ~1,648 ns/site due to isinstance checks, dict lookups, and method dispatch, while the underlying math (3 × mersenne_mul ≈ 3 ns) is negligible.

## PyO3 Architecture Decision

When implemented, use **coarse-grained bindings** (batch operations cross the FFI boundary once) rather than fine-grained (per-site FFI calls). 

Estimated per-call FFI overhead: ~35–70 ns (CPython dispatch + PyO3 handle_panic + arg extraction + return conversion).

**Fine-grained** (23 FFI calls): ~1,300 ns total — still 23× faster than pure Python hashrope
**Coarse-grained** (1 FFI call): ~165 ns total — 330× faster than pure Python hashrope, 12× faster than Python byte-slice baseline

Key design: wrap `Arena` as opaque `#[pyclass(frozen)]`, pass `u32` NodeIds as primitives, never expose per-node Python objects.

## Priority

After pure-Rust benchmarks confirm the algorithm's true constant factors, then implement PyO3 bindings to close the Python performance gap.
