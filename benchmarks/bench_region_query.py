"""E-G1: Genome-Scale substr_hash Query Performance.

Measures region identity query time at various region sizes
on a pre-built genome rope, compared to samtools + hash and
direct read + hash baselines.

Usage:
    python benchmarks/bench_region_query.py --fasta data/chr22.fa --iterations 1000
"""

from __future__ import annotations

import argparse
import random
import time
from pathlib import Path

# TODO: implement when hashrope_bio.genomics is ready


def bench_region_queries(
    fasta_path: Path,
    chunk_size: int,
    region_sizes: list[int],
    iterations: int,
) -> None:
    """Benchmark substr_hash queries at various region sizes.

    For each region size L:
        1. Generate `iterations` random (start, L) pairs
        2. Time substr_hash per query
        3. Time baseline: read bytes + PolynomialHash.hash()
        4. Time baseline: samtools faidx + SHA-256 (if pysam available)

    Reports: region_size, hashrope_ns, read_hash_ns, samtools_ns, speedup
    """
    random.seed(42)

    print(f"Building rope from {fasta_path} (chunk_size={chunk_size})...")
    # TODO: build rope

    header = f"{'region_L':>10} {'hashrope_ns':>14} {'read+hash_ns':>14} {'speedup':>10}"
    print(header)
    print("-" * len(header))

    for L in region_sizes:
        # TODO: implement query benchmark
        print(f"{L:>10} {'TODO':>14} {'TODO':>14} {'TODO':>10}")


def main():
    parser = argparse.ArgumentParser(description="E-G1: Region Query Benchmark")
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--chunk-size", type=int, default=4096)
    parser.add_argument("--iterations", type=int, default=1000)
    parser.add_argument("--region-sizes", type=str, default="100,500,1000,5000,10000,50000,100000")
    args = parser.parse_args()

    region_sizes = [int(x) for x in args.region_sizes.split(",")]
    bench_region_queries(args.fasta, args.chunk_size, region_sizes, args.iterations)


if __name__ == "__main__":
    main()
