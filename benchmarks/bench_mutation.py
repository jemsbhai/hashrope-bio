"""E-G3: Mutation Localization via Binary Search.

Introduces synthetic SNPs into a reference, then uses binary search
via substr_hash to pinpoint mutation positions. Measures comparison
count and wall-clock time vs linear scan.

Usage:
    python benchmarks/bench_mutation.py --fasta data/chr22.fa --mutations 1 --trials 100
"""

from __future__ import annotations

import argparse
import random
import time
from pathlib import Path

# TODO: implement when hashrope_bio.genomics is ready


def bench_mutation_localization(
    fasta_path: Path,
    region_sizes: list[int],
    num_mutations: int,
    trials: int,
) -> None:
    """Benchmark mutation localization at various region sizes.

    For each region size N:
        1. Build reference rope
        2. Introduce `num_mutations` random SNPs → build sample rope
        3. Run binary search — count comparisons, measure total time
        4. Baseline: linear scan byte-by-byte
        5. Repeat `trials` times with different mutation positions

    Reports: N, comparisons, hashrope_us, linear_us, speedup
    """
    random.seed(42)

    header = f"{'N':>12} {'comparisons':>14} {'hashrope_us':>14} {'linear_us':>14} {'speedup':>10}"
    print(header)
    print("-" * len(header))

    for N in region_sizes:
        # TODO: implement
        print(f"{N:>12} {'TODO':>14} {'TODO':>14} {'TODO':>14} {'TODO':>10}")


def main():
    parser = argparse.ArgumentParser(description="E-G3: Mutation Localization Benchmark")
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--region-sizes", type=str, default="1000,10000,100000,1000000,10000000")
    parser.add_argument("--mutations", type=int, default=1)
    parser.add_argument("--trials", type=int, default=100)
    args = parser.parse_args()

    region_sizes = [int(x) for x in args.region_sizes.split(",")]
    bench_mutation_localization(args.fasta, region_sizes, args.mutations, args.trials)


if __name__ == "__main__":
    main()
