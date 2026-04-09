#!/usr/bin/env python3
"""E-CG1: Cross-species gene conservation screening (human vs chimp chr22).

Uses hashrope to screen 643 chr22 genes for human-chimp divergence using
the UCSC LASTZ pairwise alignment. For each gene:
  1. Find alignment blocks overlapping the gene
  2. Build hashrope for human and chimp aligned sequences
  3. Compare hashes: identical = conserved, different = diverged
  4. For diverged genes: binary search to localize first divergence point
  5. Compare speedup vs linear byte-by-byte baseline

Input:
  - data/hg38_panTro6_chr22.axt     — UCSC alignment blocks
  - data/chr22_genes.tsv             — 643 gene annotations from UCSC refGene

Output: results/cross_species_conservation.json (+ timestamped archive)
"""

import json
import sys
import time
import datetime
import platform
from pathlib import Path
from dataclasses import dataclass

try:
    from hashrope import PolynomialHash, Leaf, rope_concat, rope_substr_hash
except ImportError:
    print("ERROR: hashrope package required. Install with: pip install hashrope")
    sys.exit(1)

CODE_DIR = Path(__file__).parent.parent
DATA_DIR = CODE_DIR / "data"
RESULTS_DIR = CODE_DIR / "results"
RESULTS_DIR.mkdir(exist_ok=True)

CHUNK_SIZE = 256  # best query performance from E-G1


@dataclass
class AxtBlock:
    """One pairwise alignment block from the axt file."""
    align_id: int
    t_chr: str      # target (human) chromosome
    t_start: int    # 0-based start
    t_end: int      # 0-based exclusive end
    q_chr: str      # query (chimp) chromosome
    q_start: int
    q_end: int
    strand: str
    score: int
    t_seq: str      # human aligned sequence (with gaps)
    q_seq: str      # chimp aligned sequence (with gaps)


@dataclass
class Gene:
    """Gene annotation."""
    name: str
    strand: str
    tx_start: int   # 0-based
    tx_end: int     # 0-based exclusive
    tx_len: int


def load_axt(path: Path) -> list[AxtBlock]:
    """Load axt alignment blocks."""
    blocks = []
    with open(path) as f:
        state = "header"
        header_parts = None
        t_seq = ""
        for line in f:
            line = line.strip()
            if state == "header":
                if not line:
                    continue
                parts = line.split()
                if len(parts) >= 9:
                    header_parts = parts
                    state = "target_seq"
            elif state == "target_seq":
                t_seq = line.upper()
                state = "query_seq"
            elif state == "query_seq":
                q_seq = line.upper()
                blocks.append(AxtBlock(
                    align_id=int(header_parts[0]),
                    t_chr=header_parts[1],
                    t_start=int(header_parts[2]),
                    t_end=int(header_parts[3]),
                    q_chr=header_parts[4],
                    q_start=int(header_parts[5]),
                    q_end=int(header_parts[6]),
                    strand=header_parts[7],
                    score=int(header_parts[8]),
                    t_seq=t_seq,
                    q_seq=q_seq,
                ))
                state = "blank"
            elif state == "blank":
                state = "header"
    return blocks


def load_genes(path: Path) -> list[Gene]:
    """Load gene annotations TSV."""
    genes = []
    with open(path) as f:
        f.readline()  # skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
            genes.append(Gene(
                name=parts[0],
                strand=parts[2],
                tx_start=int(parts[3]),
                tx_end=int(parts[4]),
                tx_len=int(parts[5]),
            ))
    return genes


def find_overlapping_blocks(gene: Gene, blocks: list[AxtBlock]) -> list[AxtBlock]:
    """Find alignment blocks that overlap a gene (binary search would be better
    but linear is fine for 3,383 blocks × 643 genes — ~2M comparisons)."""
    overlapping = []
    for b in blocks:
        if b.t_start < gene.tx_end and b.t_end > gene.tx_start:
            overlapping.append(b)
    return overlapping


def extract_aligned_region(block: AxtBlock, start: int, end: int) -> tuple[str, str]:
    """Extract aligned human and chimp sequences for a genomic region.

    The axt sequences include gap characters ('-'). Positions in the
    aligned sequence correspond to positions in the target (human) genome,
    but gaps in the target mean the aligned sequence is longer than the
    genomic span.

    We need to map genomic coordinates to alignment positions, accounting
    for gaps in the target.
    """
    # Clamp to block boundaries
    region_start = max(start, block.t_start)
    region_end = min(end, block.t_end)
    if region_start >= region_end:
        return "", ""

    # Walk through the aligned sequence, tracking genomic position
    genome_pos = block.t_start
    align_start = None
    align_end = None

    for i, (t_char, q_char) in enumerate(zip(block.t_seq, block.q_seq)):
        if t_char != '-':
            # This alignment position corresponds to a genomic position
            if genome_pos == region_start and align_start is None:
                align_start = i
            genome_pos += 1
            if genome_pos == region_end:
                align_end = i + 1
                break

    if align_start is None:
        return "", ""
    if align_end is None:
        align_end = len(block.t_seq)

    return block.t_seq[align_start:align_end], block.q_seq[align_start:align_end]


def _build_rope(data: bytes, chunk_size: int, h: PolynomialHash):
    """Build a chunked rope from bytes."""
    rope = None
    for i in range(0, len(data), chunk_size):
        chunk = data[i:i + chunk_size]
        leaf = Leaf(chunk, h)
        rope = rope_concat(rope, leaf, h)
    return rope


def compare_gene(human_seq: str, chimp_seq: str, h: PolynomialHash) -> dict:
    """Compare human and chimp aligned sequences using hashrope.

    Returns dict with comparison results and timing.
    """
    assert len(human_seq) == len(chimp_seq), \
        f"Aligned sequences must be same length: {len(human_seq)} vs {len(chimp_seq)}"

    seq_len = len(human_seq)
    human_bytes = human_seq.encode("ascii")
    chimp_bytes = chimp_seq.encode("ascii")

    # --- Hashrope comparison ---
    t0 = time.perf_counter_ns()
    human_rope = _build_rope(human_bytes, CHUNK_SIZE, h)
    chimp_rope = _build_rope(chimp_bytes, CHUNK_SIZE, h)

    # Full-sequence hash comparison
    human_hash = rope_substr_hash(human_rope, 0, seq_len, h)
    chimp_hash = rope_substr_hash(chimp_rope, 0, seq_len, h)
    is_conserved = (human_hash == chimp_hash)
    hashrope_ns = time.perf_counter_ns() - t0

    # --- Byte baseline ---
    t0 = time.perf_counter_ns()
    byte_match = (human_bytes == chimp_bytes)
    baseline_ns = time.perf_counter_ns() - t0

    # Verify hash correctness
    assert is_conserved == byte_match, \
        f"Hash mismatch! hashrope={is_conserved}, bytes={byte_match}"

    result = {
        "seq_len": seq_len,
        "is_conserved": is_conserved,
        "hashrope_ns": hashrope_ns,
        "baseline_ns": baseline_ns,
    }

    # Binary search to localize first divergence if diverged
    if not is_conserved:
        t0 = time.perf_counter_ns()
        lo, hi = 0, seq_len
        comparisons = 0
        while lo < hi - 1:
            mid = (lo + hi) // 2
            h_hash = rope_substr_hash(human_rope, lo, mid - lo, h)
            c_hash = rope_substr_hash(chimp_rope, lo, mid - lo, h)
            comparisons += 1
            if h_hash != c_hash:
                hi = mid
            else:
                lo = mid
        search_ns = time.perf_counter_ns() - t0

        # Count total mismatches (linear, for ground truth)
        t0_count = time.perf_counter_ns()
        mismatches = sum(1 for a, b in zip(human_bytes, chimp_bytes) if a != b)
        count_ns = time.perf_counter_ns() - t0_count

        result["first_divergence_pos"] = lo
        result["search_comparisons"] = comparisons
        result["search_ns"] = search_ns
        result["total_mismatches"] = mismatches
        result["mismatch_rate_pct"] = round(mismatches / seq_len * 100, 4)
        result["count_ns"] = count_ns

    return result


def main():
    print("=" * 70)
    print("E-CG1: Cross-Species Gene Conservation Screening")
    print("  Human (hg38) vs Chimpanzee (panTro6) — chr22")
    print("=" * 70)

    # Load data
    print("\nLoading alignment blocks...")
    axt_path = DATA_DIR / "hg38_panTro6_chr22.axt"
    blocks = load_axt(axt_path)
    # Sort by target start for efficient overlap search
    blocks.sort(key=lambda b: b.t_start)
    print(f"  Loaded {len(blocks):,} alignment blocks")

    print("Loading gene annotations...")
    gene_path = DATA_DIR / "chr22_genes.tsv"
    genes = load_genes(gene_path)
    print(f"  Loaded {len(genes)} genes")

    h = PolynomialHash()

    # Screen each gene
    print(f"\nScreening {len(genes)} genes for human-chimp conservation...")
    results = []
    conserved_count = 0
    diverged_count = 0
    no_alignment_count = 0
    total_hashrope_ns = 0
    total_baseline_ns = 0
    hash_errors = 0

    for i, gene in enumerate(genes):
        overlapping = find_overlapping_blocks(gene, blocks)
        if not overlapping:
            no_alignment_count += 1
            results.append({
                "gene": gene.name,
                "tx_start": gene.tx_start,
                "tx_end": gene.tx_end,
                "tx_len": gene.tx_len,
                "status": "no_alignment",
            })
            continue

        # Concatenate aligned sequences across overlapping blocks
        all_human = []
        all_chimp = []
        total_aligned = 0
        for block in overlapping:
            h_seq, c_seq = extract_aligned_region(block, gene.tx_start, gene.tx_end)
            if h_seq:
                all_human.append(h_seq)
                all_chimp.append(c_seq)
                total_aligned += len(h_seq)

        if total_aligned == 0:
            no_alignment_count += 1
            results.append({
                "gene": gene.name,
                "tx_start": gene.tx_start,
                "tx_end": gene.tx_end,
                "tx_len": gene.tx_len,
                "status": "no_aligned_bases",
            })
            continue

        human_concat = "".join(all_human)
        chimp_concat = "".join(all_chimp)

        try:
            comp = compare_gene(human_concat, chimp_concat, h)
        except AssertionError as e:
            hash_errors += 1
            results.append({
                "gene": gene.name,
                "status": "error",
                "error": str(e),
            })
            continue

        total_hashrope_ns += comp["hashrope_ns"]
        total_baseline_ns += comp["baseline_ns"]

        if comp["is_conserved"]:
            conserved_count += 1
        else:
            diverged_count += 1

        result = {
            "gene": gene.name,
            "tx_start": gene.tx_start,
            "tx_end": gene.tx_end,
            "tx_len": gene.tx_len,
            "aligned_len": total_aligned,
            "n_blocks": len(overlapping),
            "status": "conserved" if comp["is_conserved"] else "diverged",
            **comp,
        }
        results.append(result)

        if (i + 1) % 100 == 0:
            print(f"  {i + 1}/{len(genes)} genes processed...")

    # Summary
    print(f"\n{'=' * 70}")
    print("Results:")
    print(f"  Genes screened:     {len(genes)}")
    print(f"  With alignment:     {conserved_count + diverged_count}")
    print(f"  No alignment:       {no_alignment_count}")
    print(f"  Hash errors:        {hash_errors}")
    print(f"  ---")
    print(f"  Conserved (identical): {conserved_count}")
    print(f"  Diverged:              {diverged_count}")
    pct_conserved = conserved_count / (conserved_count + diverged_count) * 100 \
        if (conserved_count + diverged_count) > 0 else 0
    print(f"  Conservation rate:     {pct_conserved:.1f}%")
    print(f"  ---")
    print(f"  Total hashrope time: {total_hashrope_ns / 1e6:.1f} ms")
    print(f"  Total baseline time: {total_baseline_ns / 1e6:.1f} ms")

    # Top diverged genes by mismatch count
    diverged_genes = [r for r in results if r.get("status") == "diverged"
                      and "total_mismatches" in r]
    diverged_genes.sort(key=lambda r: r.get("total_mismatches", 0), reverse=True)

    if diverged_genes:
        print(f"\n  Top 10 most diverged genes:")
        print(f"  {'Gene':<15} {'Aligned':>8} {'Mismatches':>11} {'Rate':>7} "
              f"{'Search comps':>13}")
        print(f"  {'-'*15:<15} {'-'*8:>8} {'-'*11:>11} {'-'*7:>7} {'-'*13:>13}")
        for r in diverged_genes[:10]:
            print(f"  {r['gene']:<15} {r['aligned_len']:>8,} "
                  f"{r['total_mismatches']:>11,} "
                  f"{r['mismatch_rate_pct']:>6.2f}% "
                  f"{r.get('search_comparisons', 0):>13}")

    # Conserved genes (potentially interesting — under strong purifying selection)
    conserved_genes = [r for r in results if r.get("status") == "conserved"]
    if conserved_genes:
        conserved_genes.sort(key=lambda r: r.get("aligned_len", 0), reverse=True)
        print(f"\n  Conserved genes (identical human-chimp, top 10 by size):")
        print(f"  {'Gene':<15} {'Aligned bp':>11}")
        for r in conserved_genes[:10]:
            print(f"  {r['gene']:<15} {r['aligned_len']:>11,}")

    # Speedup analysis for diverged genes
    if diverged_genes:
        search_genes = [r for r in diverged_genes if "search_ns" in r and "count_ns" in r]
        if search_genes:
            print(f"\n  Binary search localization (diverged genes):")
            print(f"  {'Gene':<15} {'Aligned':>8} {'Search':>10} {'Linear':>10} "
                  f"{'Speedup':>8} {'Comps':>6}")
            for r in search_genes[:10]:
                speedup = r["count_ns"] / r["search_ns"] if r["search_ns"] > 0 else 0
                print(f"  {r['gene']:<15} {r['aligned_len']:>8,} "
                      f"{r['search_ns']/1000:>9.1f}µs "
                      f"{r['count_ns']/1000:>9.1f}µs "
                      f"{speedup:>7.1f}× "
                      f"{r['search_comparisons']:>6}")

    print(f"{'=' * 70}")

    # Save results
    output = {
        "benchmark": "cross_species_conservation",
        "environment": {
            "timestamp_iso": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "os": platform.system(),
            "python_version": platform.python_version(),
            "language": "python",
        },
        "config": {
            "human_assembly": "hg38 (GRCh38)",
            "chimp_assembly": "panTro6 (Clint_PTRv2)",
            "alignment_source": "UCSC LASTZ net axt",
            "alignment_citation": "Kent WJ et al. PNAS 2003;100(20):11484-11489",
            "chromosome": "chr22",
            "n_alignment_blocks": len(blocks),
            "n_genes": len(genes),
            "chunk_size": CHUNK_SIZE,
        },
        "results": {
            "genes_with_alignment": conserved_count + diverged_count,
            "no_alignment": no_alignment_count,
            "hash_errors": hash_errors,
            "conserved": conserved_count,
            "diverged": diverged_count,
            "conservation_rate_pct": round(pct_conserved, 2),
            "total_hashrope_ns": total_hashrope_ns,
            "total_baseline_ns": total_baseline_ns,
        },
        "per_gene": results,
    }

    out_path = RESULTS_DIR / "cross_species_conservation.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved: {out_path}")

    ts = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    archive = RESULTS_DIR / f"cross_species_conservation_{ts}.json"
    with open(archive, "w") as f:
        json.dump(output, f, indent=2)
    print(f"Saved: {archive}")


if __name__ == "__main__":
    main()
