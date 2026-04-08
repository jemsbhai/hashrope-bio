"""E-D2: Drug Resistance Mutation Panel.

Benchmarks checking N known resistance sites via substr_hash
on HIV-1 reverse transcriptase gene.

Usage:
    python benchmarks/bench_resistance.py
    python benchmarks/bench_resistance.py --fasta data/hiv1_hxb2.fa
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

from hashrope import PolynomialHash, Leaf, rope_substr_hash

from hashrope_bio.cheminformatics.resistance import (
    check_resistance_panel,
    HIV_RT_NNRTI_PANEL,
    HIV_RT_NRTI_PANEL,
)


def bench_synthetic_panel() -> None:
    """Benchmark resistance panel on synthetic gene data."""
    h = PolynomialHash()

    # Synthetic RT gene: 560 amino acids × 3 nt = 1680 bp
    rt_length = 1680
    ref_gene = bytes([ord("A") + (i % 4) for i in range(rt_length)])

    # Introduce mutations at K103 (pos 103) and M184 (pos 184)
    mut_gene = bytearray(ref_gene)
    mut_gene[102 * 3] = ord("T")      # K103N: codon 103 mutated
    mut_gene[102 * 3 + 1] = ord("T")
    mut_gene[183 * 3] = ord("G")      # M184V: codon 184 mutated
    mut_gene = bytes(mut_gene)

    ref_rope = Leaf(ref_gene, h)
    mut_rope = Leaf(mut_gene, h)

    combined_panel = HIV_RT_NNRTI_PANEL + HIV_RT_NRTI_PANEL

    # Warm up
    for _ in range(100):
        check_resistance_panel(ref_rope, mut_rope, combined_panel, h)

    # Timed run
    iterations = 10_000
    t0 = time.perf_counter_ns()
    for _ in range(iterations):
        results = check_resistance_panel(ref_rope, mut_rope, combined_panel, h)
    t1 = time.perf_counter_ns()

    total_ns = t1 - t0
    per_check_ns = total_ns / iterations
    per_site_ns = per_check_ns / len(combined_panel)

    print(f"Panel size:      {len(combined_panel)} sites")
    print(f"Gene length:     {rt_length} bp")
    print(f"Total per check: {per_check_ns:.0f} ns ({per_check_ns / 1000:.1f} µs)")
    print(f"Per site:        {per_site_ns:.0f} ns")
    print(f"Mutations found: {sum(1 for r in results if r.is_mutant)}")
    print()

    # Baseline: extract codons and compare bytes
    t0 = time.perf_counter_ns()
    for _ in range(iterations):
        for site in combined_panel:
            offset = (site.position - 1) * 3
            _ = ref_gene[offset:offset + 3] == mut_gene[offset:offset + 3]
    t1 = time.perf_counter_ns()

    baseline_ns = (t1 - t0) / iterations
    baseline_per_site = baseline_ns / len(combined_panel)

    print(f"Baseline total:  {baseline_ns:.0f} ns ({baseline_ns / 1000:.1f} µs)")
    print(f"Baseline/site:   {baseline_per_site:.0f} ns")
    print(f"Ratio:           {per_check_ns / baseline_ns:.1f}×")


def main():
    parser = argparse.ArgumentParser(description="E-D2: Resistance Panel Benchmark")
    parser.add_argument("--fasta", type=Path, default=None,
                        help="Optional: real HIV-1 HXB2 FASTA. Uses synthetic data if not provided.")
    args = parser.parse_args()

    if args.fasta:
        print(f"TODO: implement real FASTA panel check from {args.fasta}")
    else:
        bench_synthetic_panel()


if __name__ == "__main__":
    main()
