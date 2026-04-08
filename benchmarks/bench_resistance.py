"""E-D2: Drug Resistance Mutation Panel Benchmark.

Benchmarks checking N known resistance sites via substr_hash on
HIV-1 reverse transcriptase gene.

Runs on both:
    1. Synthetic data (controlled baseline, always works)
    2. Real HIV-1 HXB2 reference sequence (if downloaded)

Usage:
    python benchmarks/bench_resistance.py
    python benchmarks/bench_resistance.py --output ../../results/
"""

from __future__ import annotations

import argparse
import sys
import time
from dataclasses import dataclass
from pathlib import Path

from hashrope import PolynomialHash, Leaf, rope_concat, rope_hash, rope_len, rope_substr_hash

from hashrope_bio.genomics.fasta import load_fasta_to_rope, load_fasta_bytes
from hashrope_bio.cheminformatics.resistance import (
    check_resistance_panel,
    ResistanceSite,
    HIV_RT_NNRTI_PANEL,
    HIV_RT_NRTI_PANEL,
)
from hashrope_bio.result_output import save_results


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

WARMUP_ITERS = 500
TIMING_ITERS = 50_000

# HIV-1 HXB2 RT gene coordinates (GenBank K03455.1)
HIV_RT_START = 2549  # 0-based
HIV_RT_END = 4229    # exclusive (1680 bp = 560 codons)

DATA_DIR = Path(__file__).parent.parent / "data"


@dataclass
class PanelBenchResult:
    """Result of one panel benchmark run."""
    data_source: str
    gene_length_bp: int
    panel_size: int
    mutations_introduced: int
    mutations_detected: int
    hashrope_total_ns: float
    hashrope_per_site_ns: float
    baseline_total_ns: float
    baseline_per_site_ns: float
    speedup: float
    all_correct: bool
    details: str = ""


# ---------------------------------------------------------------------------
# Benchmarks
# ---------------------------------------------------------------------------


def bench_panel(
    ref_bytes: bytes,
    panel: list[ResistanceSite],
    mutation_positions: list[int],
    data_source: str,
    h: PolynomialHash,
) -> PanelBenchResult:
    """Benchmark a resistance panel on given gene bytes."""
    gene_len = len(ref_bytes)

    mut_bytes = bytearray(ref_bytes)
    for pos in mutation_positions:
        offset = (pos - 1) * 3
        original = mut_bytes[offset]
        mut_bytes[offset] = ord("G") if original != ord("G") else ord("A")
    mut_bytes = bytes(mut_bytes)

    ref_rope = Leaf(ref_bytes, h)
    mut_rope = Leaf(mut_bytes, h)

    # Hashrope panel check
    for _ in range(WARMUP_ITERS):
        results = check_resistance_panel(ref_rope, mut_rope, panel, h)

    t0 = time.perf_counter_ns()
    for _ in range(TIMING_ITERS):
        results = check_resistance_panel(ref_rope, mut_rope, panel, h)
    t1 = time.perf_counter_ns()
    hashrope_total = (t1 - t0) / TIMING_ITERS
    hashrope_per_site = hashrope_total / len(panel)

    # Baseline: byte slice comparison
    for _ in range(WARMUP_ITERS):
        for site in panel:
            offset = (site.position - 1) * 3
            _ = ref_bytes[offset:offset + 3] == mut_bytes[offset:offset + 3]

    t0 = time.perf_counter_ns()
    for _ in range(TIMING_ITERS):
        for site in panel:
            offset = (site.position - 1) * 3
            _ = ref_bytes[offset:offset + 3] == mut_bytes[offset:offset + 3]
    t1 = time.perf_counter_ns()
    baseline_total = (t1 - t0) / TIMING_ITERS
    baseline_per_site = baseline_total / len(panel)

    detected = {r.site.position for r in results if r.is_mutant}
    expected = set(mutation_positions)
    all_correct = detected == expected
    details = "all correct" if all_correct else f"FN: {expected - detected}, FP: {detected - expected}"

    return PanelBenchResult(
        data_source=data_source, gene_length_bp=gene_len, panel_size=len(panel),
        mutations_introduced=len(mutation_positions), mutations_detected=len(detected),
        hashrope_total_ns=hashrope_total, hashrope_per_site_ns=hashrope_per_site,
        baseline_total_ns=baseline_total, baseline_per_site_ns=baseline_per_site,
        speedup=baseline_total / hashrope_total if hashrope_total > 0 else 0,
        all_correct=all_correct, details=details,
    )


def bench_panel_on_chunked_rope(
    ref_bytes: bytes,
    panel: list[ResistanceSite],
    mutation_positions: list[int],
    chunk_size: int,
    data_source: str,
    h: PolynomialHash,
) -> PanelBenchResult:
    """Same as bench_panel but with a multi-leaf chunked rope."""
    gene_len = len(ref_bytes)

    mut_bytes = bytearray(ref_bytes)
    for pos in mutation_positions:
        offset = (pos - 1) * 3
        original = mut_bytes[offset]
        mut_bytes[offset] = ord("G") if original != ord("G") else ord("A")
    mut_bytes = bytes(mut_bytes)

    ref_rope = None
    for i in range(0, len(ref_bytes), chunk_size):
        leaf = Leaf(ref_bytes[i:i + chunk_size], h)
        ref_rope = rope_concat(ref_rope, leaf, h) if ref_rope else leaf

    mut_rope = None
    for i in range(0, len(mut_bytes), chunk_size):
        leaf = Leaf(mut_bytes[i:i + chunk_size], h)
        mut_rope = rope_concat(mut_rope, leaf, h) if mut_rope else leaf

    for _ in range(WARMUP_ITERS):
        results = check_resistance_panel(ref_rope, mut_rope, panel, h)

    t0 = time.perf_counter_ns()
    for _ in range(TIMING_ITERS):
        results = check_resistance_panel(ref_rope, mut_rope, panel, h)
    t1 = time.perf_counter_ns()
    hashrope_total = (t1 - t0) / TIMING_ITERS
    hashrope_per_site = hashrope_total / len(panel)

    for _ in range(WARMUP_ITERS):
        for site in panel:
            offset = (site.position - 1) * 3
            _ = ref_bytes[offset:offset + 3] == mut_bytes[offset:offset + 3]

    t0 = time.perf_counter_ns()
    for _ in range(TIMING_ITERS):
        for site in panel:
            offset = (site.position - 1) * 3
            _ = ref_bytes[offset:offset + 3] == mut_bytes[offset:offset + 3]
    t1 = time.perf_counter_ns()
    baseline_total = (t1 - t0) / TIMING_ITERS
    baseline_per_site = baseline_total / len(panel)

    detected = {r.site.position for r in results if r.is_mutant}
    expected = set(mutation_positions)
    all_correct = detected == expected

    return PanelBenchResult(
        data_source=f"{data_source}_chunked_{chunk_size}",
        gene_length_bp=gene_len, panel_size=len(panel),
        mutations_introduced=len(mutation_positions), mutations_detected=len(detected),
        hashrope_total_ns=hashrope_total, hashrope_per_site_ns=hashrope_per_site,
        baseline_total_ns=baseline_total, baseline_per_site_ns=baseline_per_site,
        speedup=baseline_total / hashrope_total if hashrope_total > 0 else 0,
        all_correct=all_correct, details=f"chunk_size={chunk_size}",
    )


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


def load_hiv_rt(fasta_path: Path) -> bytes | None:
    """Extract the RT coding sequence from HIV-1 HXB2 FASTA."""
    try:
        seq_bytes, name = load_fasta_bytes(fasta_path)
        print(f"  Loaded: {name} ({len(seq_bytes)} bp)")
        if len(seq_bytes) < HIV_RT_END:
            print(f"  WARNING: sequence too short ({len(seq_bytes)} < {HIV_RT_END})")
            return None
        rt_bytes = seq_bytes[HIV_RT_START:HIV_RT_END]
        print(f"  RT gene: positions {HIV_RT_START + 1}–{HIV_RT_END} ({len(rt_bytes)} bp, {len(rt_bytes) // 3} codons)")
        return rt_bytes
    except Exception as e:
        print(f"  Failed to load: {e}")
        return None


def make_synthetic_gene(length_bp: int = 1680) -> bytes:
    """Deterministic pseudo-random gene sequence."""
    nucs = b"ACGT"
    state = 42
    gene = bytearray(length_bp)
    for i in range(length_bp):
        state = (state * 6364136223846793005 + 1442695040888963407) & 0xFFFFFFFFFFFFFFFF
        gene[i] = nucs[(state >> 33) & 3]
    return bytes(gene)


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------


def print_results(results: list[PanelBenchResult]) -> None:
    print("\n" + "=" * 120)
    print("E-D2: Drug Resistance Mutation Panel Benchmark")
    print("=" * 120)
    header = (f"{'source':>30} {'gene_bp':>8} {'panel':>6} {'muts':>5} {'det':>5} "
              f"{'hr_total_ns':>12} {'hr/site_ns':>12} {'bl_total_ns':>12} {'bl/site_ns':>12} "
              f"{'speedup':>8} {'correct':>8}")
    print(header)
    print("-" * len(header))
    for r in results:
        print(f"{r.data_source:>30} {r.gene_length_bp:>8} {r.panel_size:>6} "
              f"{r.mutations_introduced:>5} {r.mutations_detected:>5} "
              f"{r.hashrope_total_ns:>12.0f} {r.hashrope_per_site_ns:>12.0f} "
              f"{r.baseline_total_ns:>12.0f} {r.baseline_per_site_ns:>12.0f} "
              f"{r.speedup:>8.2f} {'OK' if r.all_correct else 'FAIL':>8}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description="E-D2: Resistance Panel Benchmark")
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--skip-download", action="store_true")
    args = parser.parse_args()

    print("hashrope-bio E-D2: Drug Resistance Mutation Panel")
    print(f"Python {sys.version}")
    print(f"Timing iterations: {TIMING_ITERS}, warmup: {WARMUP_ITERS}")

    h = PolynomialHash()
    combined_panel = HIV_RT_NNRTI_PANEL + HIV_RT_NRTI_PANEL
    mutation_positions = [103, 184]

    print(f"\nPanel: {len(combined_panel)} sites ({len(HIV_RT_NNRTI_PANEL)} NNRTI + {len(HIV_RT_NRTI_PANEL)} NRTI)")
    print(f"Mutations to introduce: positions {mutation_positions}")

    results: list[PanelBenchResult] = []

    # --- 1. Synthetic ---
    print("\n--- Synthetic gene (1680 bp, deterministic pseudo-random) ---")
    synthetic_gene = make_synthetic_gene(1680)

    results.append(bench_panel(synthetic_gene, combined_panel, mutation_positions, "synthetic", h))
    for cs in [64, 256]:
        results.append(bench_panel_on_chunked_rope(
            synthetic_gene, combined_panel, mutation_positions, cs, "synthetic", h))

    # --- 2. Real HIV-1 HXB2 ---
    hxb2_path = DATA_DIR / "hiv1_hxb2.fa"

    if not hxb2_path.exists() and not args.skip_download:
        print(f"\n--- Downloading HIV-1 HXB2 reference ---")
        try:
            import urllib.request
            DATA_DIR.mkdir(exist_ok=True)
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=K03455.1&rettype=fasta&retmode=text"
            req = urllib.request.Request(url, headers={"User-Agent": "hashrope-bio/0.1"})
            with urllib.request.urlopen(req) as resp:
                data = resp.read()
            with open(hxb2_path, "wb") as f:
                f.write(data)
            print(f"  Downloaded: {hxb2_path} ({len(data)} bytes)")
        except Exception as e:
            print(f"  Download failed: {e}")

    if hxb2_path.exists():
        print(f"\n--- Real data: HIV-1 HXB2 (GenBank K03455.1) ---")
        rt_bytes = load_hiv_rt(hxb2_path)

        if rt_bytes is not None:
            results.append(bench_panel(rt_bytes, combined_panel, mutation_positions, "hiv1_hxb2_rt", h))
            for cs in [64, 256]:
                results.append(bench_panel_on_chunked_rope(
                    rt_bytes, combined_panel, mutation_positions, cs, "hiv1_hxb2_rt", h))

            print(f"\n  --- Real codon validation (first 5 panel sites) ---")
            for site in combined_panel[:5]:
                offset = (site.position - 1) * 3
                codon = rt_bytes[offset:offset + 3].decode("ascii")
                print(f"    Position {site.position:>4} ({site.annotation}): codon = {codon}")
    else:
        print(f"\n  Skipping real data — download with: python scripts/download_data.py --quick")

    # --- Report ---
    print_results(results)

    failures = [r for r in results if not r.all_correct]
    if failures:
        print(f"\nFAILED: {len(failures)} correctness failures!", file=sys.stderr)
        for f in failures:
            print(f"  {f.data_source}: {f.details}", file=sys.stderr)
        sys.exit(1)
    else:
        print(f"\nAll {len(results)} benchmarks passed correctness checks.")

    if args.output:
        config = {
            "warmup_iters": WARMUP_ITERS,
            "timing_iters": TIMING_ITERS,
            "panel_sites": len(combined_panel),
            "panel_nnrti": len(HIV_RT_NNRTI_PANEL),
            "panel_nrti": len(HIV_RT_NRTI_PANEL),
            "mutation_positions": mutation_positions,
            "hiv_rt_start_0based": HIV_RT_START,
            "hiv_rt_end_exclusive": HIV_RT_END,
            "hash_base": h.base,
            "hash_prime": str(h.prime),
        }

        save_results("resistance_panel", args.output, results, config=config)


if __name__ == "__main__":
    main()
