"""E-G2: Tandem Repeat Compression via RepeatNode.

Measures RepeatNode memory savings and hash computation speed
on known STR/VNTR loci and clinically relevant repeat expansions.

Usage:
    python benchmarks/bench_repeats.py
"""

from __future__ import annotations

import sys
import time

from hashrope import PolynomialHash, Leaf, rope_repeat, rope_hash

from hashrope_bio.genomics.repeats import CLINICAL_REPEATS, build_repeat_node


def bench_repeat_compression() -> None:
    """Benchmark RepeatNode vs materialized for clinical repeat loci."""
    h = PolynomialHash()

    header = f"{'locus':>8} {'motif':>6} {'q':>8} {'repeat_ns':>12} {'naive_ns':>12} {'repeat_bytes':>14} {'naive_bytes':>14} {'ratio':>8}"
    print(header)
    print("-" * len(header))

    for name, info in CLINICAL_REPEATS.items():
        motif = info["motif"]
        # Test at normal range, pathogenic threshold, and extreme expansion
        for q in [info["normal"][1], info["pathogenic"][0], info["pathogenic"][0] * 10]:
            if q is None:
                continue

            # RepeatNode
            t0 = time.perf_counter_ns()
            node = build_repeat_node(motif, q, h)
            repeat_hash = rope_hash(node)
            t1 = time.perf_counter_ns()
            repeat_ns = t1 - t0

            # Naive: materialize and hash
            t0 = time.perf_counter_ns()
            materialized = motif * q
            naive_hash = h.hash(materialized)
            t1 = time.perf_counter_ns()
            naive_ns = t1 - t0

            # Verify hashes match
            assert repeat_hash == naive_hash, f"Hash mismatch at {name} q={q}"

            # Memory estimate (RepeatNode = 2 Python objects, naive = q*len(motif) bytes)
            repeat_bytes = sys.getsizeof(motif) + 64  # rough: leaf + repeat node overhead
            naive_bytes = len(materialized)

            ratio = naive_bytes / repeat_bytes if repeat_bytes > 0 else 0

            print(f"{name:>8} {motif.decode():>6} {q:>8} {repeat_ns:>12} {naive_ns:>12} {repeat_bytes:>14} {naive_bytes:>14} {ratio:>8.1f}")


def bench_repeat_count_comparison() -> None:
    """Benchmark: comparing two repeat counts via hash (O(1)) vs byte comparison."""
    h = PolynomialHash()
    motif = b"CAG"

    print("\n--- Repeat count comparison ---")
    header = f"{'q1':>8} {'q2':>8} {'hash_ns':>10} {'bytes_ns':>10} {'identical':>10}"
    print(header)
    print("-" * len(header))

    for q1, q2 in [(35, 35), (35, 36), (35, 70), (200, 200), (200, 2000)]:
        node1 = build_repeat_node(motif, q1, h)
        node2 = build_repeat_node(motif, q2, h)

        t0 = time.perf_counter_ns()
        hash_eq = rope_hash(node1) == rope_hash(node2)
        t1 = time.perf_counter_ns()
        hash_ns = t1 - t0

        mat1 = motif * q1
        mat2 = motif * q2
        t0 = time.perf_counter_ns()
        byte_eq = mat1 == mat2
        t1 = time.perf_counter_ns()
        bytes_ns = t1 - t0

        assert hash_eq == byte_eq
        print(f"{q1:>8} {q2:>8} {hash_ns:>10} {bytes_ns:>10} {str(hash_eq):>10}")


if __name__ == "__main__":
    bench_repeat_compression()
    bench_repeat_count_comparison()
