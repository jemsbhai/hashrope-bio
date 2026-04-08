"""E-G2: Tandem Repeat Compression via RepeatNode.

Measures RepeatNode memory savings and hash computation speed
on clinically relevant repeat expansion loci.

Runs both synthetic benchmarks (controlled) and real data validation
(actual motifs and clinically documented repeat counts).

Usage:
    python benchmarks/bench_repeats.py
    python benchmarks/bench_repeats.py --output ../../results/
"""

from __future__ import annotations

import argparse
import sys
import time
from dataclasses import dataclass
from pathlib import Path

from hashrope import PolynomialHash, Leaf, rope_repeat, rope_hash, rope_len

from hashrope_bio.genomics.repeats import CLINICAL_REPEATS, build_repeat_node
from hashrope_bio.result_output import save_results


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Repeat counts to benchmark at — spans normal, borderline, pathogenic, extreme
REPEAT_COUNTS = [5, 10, 20, 35, 50, 100, 200, 500, 1000, 2000, 5000, 10000]

# Number of timing iterations for stable measurements
WARMUP_ITERS = 100
TIMING_ITERS = 10_000
COMPARISON_ITERS = 100_000

# For large q materialization becomes expensive — skip naive baseline above this
NAIVE_MAX_Q = 100_000


@dataclass
class RepeatResult:
    """Result of one repeat benchmark measurement."""
    motif: str
    motif_len: int
    q: int
    repeat_ns: float
    naive_ns: float | None
    repeat_bytes_est: int
    naive_bytes: int
    compression_ratio: float
    hash_match: bool
    category: str  # "normal", "borderline", "pathogenic", "extreme"


@dataclass
class ComparisonResult:
    """Result of repeat-count comparison benchmark."""
    motif: str
    q1: int
    q2: int
    hash_compare_ns: float
    byte_compare_ns: float
    identical: bool


# ---------------------------------------------------------------------------
# Benchmarks
# ---------------------------------------------------------------------------


def bench_repeat_construction(h: PolynomialHash) -> list[RepeatResult]:
    """Benchmark RepeatNode construction and hash computation."""
    results = []

    for locus_name, info in CLINICAL_REPEATS.items():
        motif = info["motif"]
        normal_max = info["normal"][1]
        pathogenic_min = info["pathogenic"][0]

        for q in REPEAT_COUNTS:
            if q <= info["normal"][1]:
                category = "normal"
            elif pathogenic_min is not None and q < pathogenic_min:
                category = "borderline"
            elif pathogenic_min is not None and q < pathogenic_min * 5:
                category = "pathogenic"
            else:
                category = "extreme"

            # --- RepeatNode ---
            for _ in range(WARMUP_ITERS):
                node = build_repeat_node(motif, q, h)
                _ = rope_hash(node)

            t0 = time.perf_counter_ns()
            for _ in range(TIMING_ITERS):
                node = build_repeat_node(motif, q, h)
                _ = rope_hash(node)
            t1 = time.perf_counter_ns()
            repeat_ns = (t1 - t0) / TIMING_ITERS
            repeat_hash = rope_hash(node)

            # --- Naive baseline ---
            naive_ns = None
            naive_hash = None
            total_len = len(motif) * q

            if q <= NAIVE_MAX_Q:
                materialized = motif * q
                for _ in range(min(WARMUP_ITERS, 10)):
                    _ = h.hash(materialized)
                iters = max(1, TIMING_ITERS // max(1, q // 100))
                t0 = time.perf_counter_ns()
                for _ in range(iters):
                    _ = h.hash(materialized)
                t1 = time.perf_counter_ns()
                naive_ns = (t1 - t0) / iters
                naive_hash = h.hash(materialized)

            hash_match = True
            if naive_hash is not None:
                hash_match = repeat_hash == naive_hash
                if not hash_match:
                    print(f"  HASH MISMATCH: {locus_name} motif={motif} q={q}", file=sys.stderr)

            repeat_bytes_est = sys.getsizeof(motif) + 200
            naive_bytes = total_len
            compression_ratio = naive_bytes / repeat_bytes_est if repeat_bytes_est > 0 else 0

            results.append(RepeatResult(
                motif=motif.decode(), motif_len=len(motif), q=q,
                repeat_ns=repeat_ns, naive_ns=naive_ns,
                repeat_bytes_est=repeat_bytes_est, naive_bytes=naive_bytes,
                compression_ratio=compression_ratio,
                hash_match=hash_match, category=category,
            ))

    return results


def bench_repeat_comparison(h: PolynomialHash) -> list[ComparisonResult]:
    """Benchmark repeat-count comparison via hash vs byte comparison."""
    results = []
    motif = b"CAG"

    pairs = [
        (35, 35), (35, 36), (35, 40), (36, 36), (36, 70),
        (10, 10), (100, 100), (100, 2000), (2000, 2000),
    ]

    for q1, q2 in pairs:
        node1 = build_repeat_node(motif, q1, h)
        node2 = build_repeat_node(motif, q2, h)
        h1, h2 = rope_hash(node1), rope_hash(node2)

        for _ in range(WARMUP_ITERS):
            _ = h1 == h2

        t0 = time.perf_counter_ns()
        for _ in range(COMPARISON_ITERS):
            _ = h1 == h2
        t1 = time.perf_counter_ns()
        hash_ns = (t1 - t0) / COMPARISON_ITERS

        mat1, mat2 = motif * q1, motif * q2
        for _ in range(WARMUP_ITERS):
            _ = mat1 == mat2

        byte_iters = min(COMPARISON_ITERS, 10_000)
        t0 = time.perf_counter_ns()
        for _ in range(byte_iters):
            _ = mat1 == mat2
        t1 = time.perf_counter_ns()
        byte_ns = (t1 - t0) / byte_iters

        identical = h1 == h2
        assert identical == (mat1 == mat2), f"Disagree at q1={q1}, q2={q2}"

        results.append(ComparisonResult(
            motif=motif.decode(), q1=q1, q2=q2,
            hash_compare_ns=hash_ns, byte_compare_ns=byte_ns,
            identical=identical,
        ))

    return results


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------


def print_construction_results(results: list[RepeatResult]) -> None:
    print("\n" + "=" * 100)
    print("E-G2: RepeatNode Construction & Hash Computation")
    print("=" * 100)

    motifs_seen = []
    for r in results:
        if r.motif not in motifs_seen:
            motifs_seen.append(r.motif)

    for motif in motifs_seen:
        subset = [r for r in results if r.motif == motif]
        locus = [name for name, info in CLINICAL_REPEATS.items()
                 if info["motif"].decode() == motif]
        locus_str = ", ".join(locus) if locus else "unknown"

        print(f"\n--- Motif: {motif} (len={len(motif)}) — Loci: {locus_str} ---")
        header = f"{'q':>8} {'category':>12} {'repeat_ns':>12} {'naive_ns':>12} {'speedup':>10} {'naive_KB':>10} {'ratio':>8} {'hash_ok':>8}"
        print(header)
        print("-" * len(header))

        for r in subset:
            naive_str = f"{r.naive_ns:>12.0f}" if r.naive_ns is not None else f"{'—':>12}"
            speedup = f"{r.naive_ns / r.repeat_ns:>10.1f}×" if r.naive_ns and r.repeat_ns > 0 else f"{'—':>10}"
            print(f"{r.q:>8} {r.category:>12} {r.repeat_ns:>12.0f} {naive_str} {speedup} {r.naive_bytes / 1024:>10.1f} {r.compression_ratio:>8.1f} {'OK' if r.hash_match else 'FAIL':>8}")


def print_comparison_results(results: list[ComparisonResult]) -> None:
    print("\n" + "=" * 100)
    print("E-G2: Repeat Count Comparison (CAG — Huntington's Disease)")
    print("=" * 100)
    header = f"{'q1':>8} {'q2':>8} {'hash_ns':>10} {'byte_ns':>10} {'speedup':>10} {'identical':>10}"
    print(header)
    print("-" * len(header))
    for r in results:
        speedup = f"{r.byte_compare_ns / r.hash_compare_ns:.1f}×" if r.hash_compare_ns > 0 else "—"
        print(f"{r.q1:>8} {r.q2:>8} {r.hash_compare_ns:>10.1f} {r.byte_compare_ns:>10.1f} {speedup:>10} {str(r.identical):>10}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description="E-G2: Tandem Repeat Benchmarks")
    parser.add_argument("--output", type=Path, default=None,
                        help="Directory for JSON output (e.g., ../../results/)")
    args = parser.parse_args()

    print("hashrope-bio E-G2: Tandem Repeat Compression via RepeatNode")
    print(f"Python {sys.version}")
    print(f"Timing iterations: {TIMING_ITERS} (construction), {COMPARISON_ITERS} (comparison)")

    h = PolynomialHash()
    print(f"Hash: base={h.base}, prime=2^61-1")

    print("\n--- Clinical loci under test ---")
    for name, info in CLINICAL_REPEATS.items():
        print(f"  {name}: motif={info['motif'].decode()}, normal={info['normal']}, "
              f"pathogenic={info['pathogenic']}, disease={info['disease']}")

    construction_results = bench_repeat_construction(h)
    comparison_results = bench_repeat_comparison(h)

    print_construction_results(construction_results)
    print_comparison_results(comparison_results)

    mismatches = [r for r in construction_results if not r.hash_match]
    if mismatches:
        print(f"\nFAILED: {len(mismatches)} hash mismatches!", file=sys.stderr)
        sys.exit(1)
    else:
        print(f"\nAll {len(construction_results)} hash verifications passed.")

    if args.output:
        config = {
            "repeat_counts": REPEAT_COUNTS,
            "warmup_iters": WARMUP_ITERS,
            "timing_iters": TIMING_ITERS,
            "comparison_iters": COMPARISON_ITERS,
            "naive_max_q": NAIVE_MAX_Q,
            "hash_base": h.base,
            "hash_prime": str(h.prime),
            "clinical_loci": {name: {
                "motif": info["motif"].decode(),
                "normal_range": info["normal"],
                "pathogenic_threshold": info["pathogenic"],
                "disease": info["disease"],
            } for name, info in CLINICAL_REPEATS.items()},
        }

        save_results("repeat_construction", args.output, construction_results, config=config)
        save_results("repeat_comparison", args.output, comparison_results, config=config)


if __name__ == "__main__":
    main()
