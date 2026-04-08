"""E-G4: Rope Construction Cost and Amortization.

Measures the one-time cost of building a genome-scale rope
and computes the amortization point vs per-query savings.

Usage:
    python benchmarks/bench_construction.py --fasta data/chr22.fa
    python benchmarks/bench_construction.py --fasta data/chr22.fa --chunk-sizes 256,1024,4096,16384 --output ../../results/
"""

from __future__ import annotations

import argparse
import sys
import time
import tracemalloc
from dataclasses import dataclass
from pathlib import Path

from hashrope import (
    PolynomialHash,
    Leaf,
    Internal,
    RepeatNode,
    rope_concat,
    rope_hash,
    rope_height,
    rope_len,
    rope_substr_hash,
)

from hashrope_bio.genomics.fasta import load_fasta_bytes
from hashrope_bio.result_output import save_results


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

DEFAULT_CHUNK_SIZES = [256, 1024, 4096, 16384]
AMORT_QUERY_LEN = 10_000  # Length of region query used for amortization estimate
AMORT_QUERY_ITERS = 500   # Timing iterations for per-query cost
AMORT_WARMUP = 50


# ---------------------------------------------------------------------------
# Tree introspection
# ---------------------------------------------------------------------------


def count_nodes(node) -> dict[str, int]:
    """Count nodes by type in the rope tree.

    Returns dict with keys: leaves, internals, repeats, total.
    """
    counts = {"leaves": 0, "internals": 0, "repeats": 0}
    _count_recursive(node, counts)
    counts["total"] = counts["leaves"] + counts["internals"] + counts["repeats"]
    return counts


def _count_recursive(node, counts: dict[str, int]) -> None:
    if node is None:
        return
    if isinstance(node, Leaf):
        counts["leaves"] += 1
    elif isinstance(node, Internal):
        counts["internals"] += 1
        _count_recursive(node.left, counts)
        _count_recursive(node.right, counts)
    elif isinstance(node, RepeatNode):
        counts["repeats"] += 1
        _count_recursive(node.child, counts)


# ---------------------------------------------------------------------------
# Core benchmark
# ---------------------------------------------------------------------------


@dataclass
class ConstructionResult:
    """Result of one rope construction benchmark."""
    chunk_size: int
    seq_len: int
    leaves: int
    internals: int
    total_nodes: int
    height: int
    construction_time_s: float
    peak_memory_mb: float
    throughput_mbs: float  # MB/s of sequence data


@dataclass
class AmortizationResult:
    """Amortization analysis for a given chunk size."""
    chunk_size: int
    construction_time_s: float
    query_len: int
    hashrope_query_ns: float
    baseline_hash_ns: float
    saving_per_query_ns: float
    amortization_queries: int  # Queries needed to pay off construction


def build_rope_chunked(seq: bytes, chunk_size: int, h: PolynomialHash):
    """Build a rope from a byte sequence using fixed-size chunks.

    Returns the rope root node.
    """
    rope = None
    for i in range(0, len(seq), chunk_size):
        chunk = seq[i:i + chunk_size]
        leaf = Leaf(chunk, h)
        rope = rope_concat(rope, leaf, h)
    return rope


def bench_construction(
    seq: bytes,
    chunk_sizes: list[int],
    h: PolynomialHash,
) -> list[ConstructionResult]:
    """Benchmark rope construction at various chunk sizes."""
    results = []

    for cs in chunk_sizes:
        print(f"\n  chunk_size={cs} ...", end=" ", flush=True)

        # Measure peak memory with tracemalloc
        tracemalloc.start()
        t0 = time.perf_counter()
        rope = build_rope_chunked(seq, cs, h)
        t1 = time.perf_counter()
        _, peak_bytes = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        elapsed = t1 - t0
        height = rope_height(rope)
        nodes = count_nodes(rope)
        throughput = (len(seq) / (1024 * 1024)) / elapsed if elapsed > 0 else 0

        # Verify hash correctness against full-sequence hash
        full_hash = h.hash(seq)
        rope_h = rope_hash(rope)
        hash_ok = rope_h == full_hash
        if not hash_ok:
            print(f"HASH MISMATCH!", file=sys.stderr)
        else:
            print(f"OK ({elapsed:.3f}s, height={height}, {nodes['leaves']} leaves)")

        results.append(ConstructionResult(
            chunk_size=cs,
            seq_len=len(seq),
            leaves=nodes["leaves"],
            internals=nodes["internals"],
            total_nodes=nodes["total"],
            height=height,
            construction_time_s=elapsed,
            peak_memory_mb=peak_bytes / (1024 * 1024),
            throughput_mbs=throughput,
        ))

    return results


def bench_amortization(
    seq: bytes,
    chunk_sizes: list[int],
    construction_results: list[ConstructionResult],
    h: PolynomialHash,
) -> list[AmortizationResult]:
    """Estimate amortization point: how many queries pay off construction?

    Compares:
      - hashrope: rope_substr_hash (O(log w))
      - baseline: slice sequence + hash (O(L))
    for a fixed query length L.
    """
    results = []
    query_len = min(AMORT_QUERY_LEN, len(seq) // 2)

    for cs, cr in zip(chunk_sizes, construction_results):
        print(f"\n  Amortization for chunk_size={cs} (query_len={query_len}) ...", end=" ", flush=True)

        # Rebuild rope (construction was in a tracemalloc context)
        rope = build_rope_chunked(seq, cs, h)

        # Pick a query position in the middle of the sequence
        start = len(seq) // 4

        # --- hashrope query timing ---
        for _ in range(AMORT_WARMUP):
            rope_substr_hash(rope, start, query_len, h)

        t0 = time.perf_counter_ns()
        for _ in range(AMORT_QUERY_ITERS):
            rope_substr_hash(rope, start, query_len, h)
        t1 = time.perf_counter_ns()
        hashrope_ns = (t1 - t0) / AMORT_QUERY_ITERS

        # --- baseline: slice + hash ---
        for _ in range(AMORT_WARMUP):
            h.hash(seq[start:start + query_len])

        t0 = time.perf_counter_ns()
        for _ in range(AMORT_QUERY_ITERS):
            h.hash(seq[start:start + query_len])
        t1 = time.perf_counter_ns()
        baseline_ns = (t1 - t0) / AMORT_QUERY_ITERS

        saving_ns = baseline_ns - hashrope_ns
        if saving_ns > 0:
            construction_ns = cr.construction_time_s * 1e9
            amort_queries = int(construction_ns / saving_ns) + 1
        else:
            amort_queries = -1  # hashrope is slower — no amortization

        print(f"hashrope={hashrope_ns:.0f}ns, baseline={baseline_ns:.0f}ns, "
              f"saving={saving_ns:.0f}ns/query, amort={amort_queries} queries")

        results.append(AmortizationResult(
            chunk_size=cs,
            construction_time_s=cr.construction_time_s,
            query_len=query_len,
            hashrope_query_ns=hashrope_ns,
            baseline_hash_ns=baseline_ns,
            saving_per_query_ns=saving_ns,
            amortization_queries=amort_queries,
        ))

    return results


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------


def print_construction_results(results: list[ConstructionResult]) -> None:
    print("\n" + "=" * 100)
    print("E-G4: Rope Construction Cost")
    print("=" * 100)

    header = (f"{'chunk':>8} {'seq_len':>12} {'leaves':>8} {'internal':>8} "
              f"{'total':>8} {'height':>8} {'time_s':>10} {'mem_MB':>10} {'MB/s':>10}")
    print(header)
    print("-" * len(header))

    for r in results:
        print(f"{r.chunk_size:>8} {r.seq_len:>12,} {r.leaves:>8,} {r.internals:>8,} "
              f"{r.total_nodes:>8,} {r.height:>8} {r.construction_time_s:>10.3f} "
              f"{r.peak_memory_mb:>10.1f} {r.throughput_mbs:>10.1f}")


def print_amortization_results(results: list[AmortizationResult]) -> None:
    print("\n" + "=" * 100)
    print(f"E-G4: Amortization Analysis (query_len={results[0].query_len} bp)")
    print("=" * 100)

    header = (f"{'chunk':>8} {'construct_s':>12} {'rope_ns':>10} {'base_ns':>10} "
              f"{'save_ns':>10} {'amort_q':>10}")
    print(header)
    print("-" * len(header))

    for r in results:
        amort_str = f"{r.amortization_queries:>10,}" if r.amortization_queries > 0 else f"{'N/A':>10}"
        print(f"{r.chunk_size:>8} {r.construction_time_s:>12.3f} "
              f"{r.hashrope_query_ns:>10.0f} {r.baseline_hash_ns:>10.0f} "
              f"{r.saving_per_query_ns:>10.0f} {amort_str}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description="E-G4: Rope Construction Benchmark")
    parser.add_argument("--fasta", type=Path, required=True, help="Path to FASTA file")
    parser.add_argument("--chunk-sizes", type=str, default=",".join(str(x) for x in DEFAULT_CHUNK_SIZES),
                        help="Comma-separated chunk sizes (default: 256,1024,4096,16384)")
    parser.add_argument("--output", type=Path, default=None,
                        help="Directory for JSON output (e.g., ../../results/)")
    parser.add_argument("--skip-amortization", action="store_true",
                        help="Skip amortization benchmark (faster)")
    args = parser.parse_args()

    chunk_sizes = [int(x) for x in args.chunk_sizes.split(",")]

    print("hashrope-bio E-G4: Rope Construction Cost and Amortization")
    print(f"Python {sys.version}")
    print(f"FASTA: {args.fasta}")
    print(f"Chunk sizes: {chunk_sizes}")

    h = PolynomialHash()
    print(f"Hash: base={h.base}, prime=2^61-1")

    # Load sequence
    print(f"\nLoading sequence from {args.fasta} ...")
    t0 = time.perf_counter()
    seq, seq_name = load_fasta_bytes(args.fasta, uppercase=True)
    t1 = time.perf_counter()
    print(f"  Loaded: {seq_name}, {len(seq):,} bp ({len(seq) / (1024*1024):.1f} MB) in {t1-t0:.2f}s")

    # Construction benchmark
    print("\n--- Construction Benchmark ---")
    construction_results = bench_construction(seq, chunk_sizes, h)
    print_construction_results(construction_results)

    # Amortization benchmark
    amort_results = []
    if not args.skip_amortization:
        print("\n--- Amortization Benchmark ---")
        amort_results = bench_amortization(seq, chunk_sizes, construction_results, h)
        print_amortization_results(amort_results)

    # Verify all constructions produced same hash
    print("\n--- Hash Verification ---")
    full_hash = h.hash(seq)
    print(f"  Full sequence hash: {full_hash}")
    for cs in chunk_sizes:
        rope = build_rope_chunked(seq, cs, h)
        rh = rope_hash(rope)
        status = "OK" if rh == full_hash else "FAIL"
        print(f"  chunk_size={cs}: {status}")

    # Save results
    if args.output:
        config = {
            "fasta": str(args.fasta),
            "seq_name": seq_name,
            "seq_len": len(seq),
            "chunk_sizes": chunk_sizes,
            "amort_query_len": AMORT_QUERY_LEN,
            "amort_query_iters": AMORT_QUERY_ITERS,
            "hash_base": h.base,
            "hash_prime": str(h.prime),
        }

        save_results("construction", args.output, construction_results, config=config)
        if amort_results:
            save_results("amortization", args.output, amort_results, config=config)


if __name__ == "__main__":
    main()
