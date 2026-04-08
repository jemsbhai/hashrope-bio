"""E-G4: Rope Construction Cost and Amortization.

Measures the one-time cost of building a genome-scale rope
and computes the amortization point vs per-query savings.

Usage:
    python benchmarks/bench_construction.py --fasta data/chr22.fa --chunk-sizes 256,1024,4096,16384
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

# TODO: implement when hashrope_bio.genomics.fasta is ready


def bench_construction(fasta_path: Path, chunk_sizes: list[int]) -> None:
    """Benchmark rope construction from FASTA at various chunk sizes.

    Reports: chunk_size, leaves, nodes, height, construction_time_s, peak_memory_mb
    """
    print(f"{'chunk_size':>12} {'leaves':>8} {'height':>8} {'time_s':>10} {'mem_mb':>10}")
    print("-" * 58)

    for cs in chunk_sizes:
        # TODO: implement
        # 1. Read FASTA into chunks of size cs
        # 2. Time: build rope from chunks
        # 3. Report metrics
        print(f"{cs:>12} {'TODO':>8} {'TODO':>8} {'TODO':>10} {'TODO':>10}")


def main():
    parser = argparse.ArgumentParser(description="E-G4: Rope Construction Benchmark")
    parser.add_argument("--fasta", type=Path, required=True, help="Path to FASTA file")
    parser.add_argument("--chunk-sizes", type=str, default="256,1024,4096,16384",
                        help="Comma-separated chunk sizes")
    args = parser.parse_args()

    chunk_sizes = [int(x) for x in args.chunk_sizes.split(",")]
    bench_construction(args.fasta, chunk_sizes)


if __name__ == "__main__":
    main()
