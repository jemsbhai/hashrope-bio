#!/usr/bin/env python3
"""Download data for E-CG1: Cross-species divergence (human vs chimp chr22).

Downloads the hg38↔panTro6 whole-genome net alignment in axt format (1.6 GB
compressed), stream-filters for chr22 blocks only, and saves the chr22 portion.

The axt alignment is the UCSC LASTZ + chaining + netting pipeline — the
gold standard for pairwise whole-genome alignment in comparative genomics.
Each axt block contains properly aligned human and chimp sequences (including
gap characters), so no separate chimp FASTA is needed.

Citation: Kent WJ, Baertsch R, Hinrichs A, Miller W, Haussler D (2003).
  Evolution's cauldron: duplication, deletion, and rearrangement in the
  mouse and human genomes. PNAS 100(20):11484-11489.

Output files (in data/):
  - hg38_panTro6_chr22.axt     — chr22 alignment blocks (~60-80 MB)
"""

import gzip
import io
import sys
import time
import urllib.request
from pathlib import Path

DATA_DIR = Path(__file__).parent.parent / "data"
DATA_DIR.mkdir(exist_ok=True)

# hg38 as target (reference), panTro6 as query
AXT_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsPanTro6/hg38.panTro6.net.axt.gz"


def stream_filter_axt_chr22():
    """Download whole-genome axt, stream-filter for chr22 blocks.

    AXT format (one block = 4 lines):
      header: align_id tChr tStart tEnd qChr qStart qEnd strand score
      target_sequence (human, may contain dashes for gaps)
      query_sequence (chimp, may contain dashes for gaps)
      blank_line

    We keep blocks where tChr == "chr22".
    """
    out_path = DATA_DIR / "hg38_panTro6_chr22.axt"

    if out_path.exists():
        size_mb = out_path.stat().st_size / 1024 / 1024
        print(f"  Already exists: {out_path} ({size_mb:.1f} MB)")
        return out_path

    print(f"  Downloading + filtering: hg38-panTro6 axt alignment (~1.6 GB)")
    print(f"  Extracting only chr22 blocks (discarding the rest)")
    print(f"  URL: {AXT_URL}")
    print(f"  This will take several minutes...")

    req = urllib.request.Request(AXT_URL)
    t0 = time.time()
    total_blocks = 0
    chr22_blocks = 0

    with urllib.request.urlopen(req, timeout=600) as resp:
        total = int(resp.headers.get("Content-Length", 0))
        total_mb = total / 1024 / 1024 if total else 0
        print(f"  Remote file size: {total_mb:.0f} MB")

        decompressor = gzip.GzipFile(fileobj=resp)
        buf = io.TextIOWrapper(decompressor, encoding="ascii", errors="replace")

        with open(out_path, "w") as out_f:
            state = "header"
            current_header = ""
            current_target = ""
            current_query = ""
            is_chr22 = False
            last_report = time.time()

            for line in buf:
                line = line.rstrip("\n\r")

                if state == "header":
                    if not line or line.startswith("#"):
                        continue
                    parts = line.split()
                    if len(parts) >= 2:
                        t_chr = parts[1]
                        is_chr22 = (t_chr == "chr22")
                        total_blocks += 1
                        if is_chr22:
                            chr22_blocks += 1
                        current_header = line
                        state = "target_seq"

                        now = time.time()
                        if now - last_report > 5:
                            elapsed = now - t0
                            print(f"\r  Scanned {total_blocks:,} blocks, "
                                  f"found {chr22_blocks:,} chr22 blocks "
                                  f"({elapsed:.0f}s elapsed)", end="", flush=True)
                            last_report = now

                elif state == "target_seq":
                    current_target = line
                    state = "query_seq"

                elif state == "query_seq":
                    current_query = line
                    state = "blank"

                    if is_chr22:
                        out_f.write(f"{current_header}\n")
                        out_f.write(f"{current_target}\n")
                        out_f.write(f"{current_query}\n")
                        out_f.write("\n")

                elif state == "blank":
                    state = "header"

    elapsed = time.time() - t0
    size_mb = out_path.stat().st_size / 1024 / 1024
    print(f"\n  Done in {elapsed:.0f}s ({elapsed / 60:.1f} min)")
    print(f"  Total alignment blocks scanned: {total_blocks:,}")
    print(f"  Chr22 blocks extracted: {chr22_blocks:,}")
    print(f"  Saved: {out_path} ({size_mb:.1f} MB)")

    return out_path


def validate_axt(axt_path: Path) -> dict:
    """Validate the extracted axt file and report statistics."""
    print(f"\nValidating {axt_path.name}...")
    n_blocks = 0
    total_aligned_bp = 0
    total_matches = 0
    total_mismatches = 0
    total_human_gaps = 0
    total_chimp_gaps = 0
    min_block = float("inf")
    max_block = 0
    chimp_chr_set = set()

    with open(axt_path) as f:
        state = "header"
        current_target = ""
        for line in f:
            line = line.strip()
            if state == "header":
                if not line:
                    continue
                parts = line.split()
                if len(parts) >= 8:
                    chimp_chr_set.add(parts[4])
                    state = "target_seq"
                    n_blocks += 1
            elif state == "target_seq":
                current_target = line.upper()
                block_len = len(line)
                total_aligned_bp += block_len
                min_block = min(min_block, block_len)
                max_block = max(max_block, block_len)
                state = "query_seq"
            elif state == "query_seq":
                query = line.upper()
                # Count matches, mismatches, gaps
                for h, c in zip(current_target, query):
                    if h == '-':
                        total_human_gaps += 1
                    elif c == '-':
                        total_chimp_gaps += 1
                    elif h == c:
                        total_matches += 1
                    else:
                        total_mismatches += 1
                state = "blank"
            elif state == "blank":
                state = "header"

    identity = total_matches / (total_matches + total_mismatches) * 100 if (total_matches + total_mismatches) > 0 else 0
    divergence = 100 - identity

    print(f"  Alignment blocks: {n_blocks:,}")
    print(f"  Total aligned positions: {total_aligned_bp:,} ({total_aligned_bp / 1e6:.1f} Mbp)")
    print(f"  Block size range: {min_block:,} – {max_block:,} bp")
    print(f"  Chimp chromosomes mapped: {sorted(chimp_chr_set)[:5]}{'...' if len(chimp_chr_set) > 5 else ''}")
    print(f"  ---")
    print(f"  Matches:           {total_matches:,}")
    print(f"  Mismatches (SNVs): {total_mismatches:,}")
    print(f"  Human gaps (ins):  {total_human_gaps:,}")
    print(f"  Chimp gaps (del):  {total_chimp_gaps:,}")
    print(f"  Sequence identity: {identity:.3f}%")
    print(f"  Divergence:        {divergence:.3f}%")

    return {
        "n_blocks": n_blocks,
        "total_aligned_bp": total_aligned_bp,
        "min_block": min_block if min_block != float("inf") else 0,
        "max_block": max_block,
        "matches": total_matches,
        "mismatches": total_mismatches,
        "human_gaps": total_human_gaps,
        "chimp_gaps": total_chimp_gaps,
        "identity_pct": round(identity, 4),
        "divergence_pct": round(divergence, 4),
    }


def main():
    print("=" * 70)
    print("E-CG1 Data Download: Human vs Chimp chr22 (hg38 vs panTro6)")
    print("=" * 70)
    print()
    print("Source: UCSC LASTZ + chaining + netting pairwise alignment")
    print("  Human: hg38 (GRCh38, Dec 2013)")
    print("  Chimp: panTro6 (Clint_PTRv2, Jan 2018)")

    # Download axt alignment filtered for chr22
    print("\nDownloading + filtering alignment for chr22...")
    axt_path = stream_filter_axt_chr22()

    # Validate
    print("\n" + "=" * 70)
    print("Validation & Statistics")
    print("=" * 70)
    stats = validate_axt(axt_path)

    # Summary
    print(f"\n{'=' * 70}")
    print("Summary:")
    print(f"  Alignment: {stats['n_blocks']:,} blocks, "
          f"{stats['total_aligned_bp']:,} aligned positions "
          f"({stats['total_aligned_bp'] / 1e6:.1f} Mbp)")
    print(f"  Identity: {stats['identity_pct']:.3f}% "
          f"(divergence: {stats['divergence_pct']:.3f}%)")
    print(f"  File: {axt_path}")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
