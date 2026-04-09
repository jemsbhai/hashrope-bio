#!/usr/bin/env python3
"""Download UCSC refGene annotations and extract chr22 gene/exon coordinates.

Downloads refGene.txt.gz from UCSC hg38 database, extracts chr22 entries,
and saves a TSV with gene coordinates for E-G1 real-data benchmarking.

Output: data/chr22_genes.tsv
Columns: gene_name, transcript, strand, tx_start, tx_end, tx_len,
         cds_start, cds_end, exon_count, exon_starts, exon_ends, exon_lengths

Coordinates are 0-based (UCSC convention), matching chr22.fa indexing.
"""

import gzip
import io
import sys
import urllib.request
from pathlib import Path

DATA_DIR = Path(__file__).parent.parent / "data"
DATA_DIR.mkdir(exist_ok=True)

REFGENE_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz"
OUTPUT_TSV = DATA_DIR / "chr22_genes.tsv"
REFGENE_GZ = DATA_DIR / "refGene.txt.gz"


def download_refgene():
    """Download refGene.txt.gz if not already present."""
    if REFGENE_GZ.exists():
        print(f"  Already exists: {REFGENE_GZ}")
        return
    print(f"  Downloading from {REFGENE_URL}...")
    urllib.request.urlretrieve(REFGENE_URL, REFGENE_GZ)
    print(f"  Saved: {REFGENE_GZ} ({REFGENE_GZ.stat().st_size / 1024:.0f} KB)")


def extract_chr22_genes():
    """Parse refGene.txt.gz and extract chr22 entries.

    refGene columns (genePred extended format):
      0: bin, 1: name (transcript), 2: chrom, 3: strand,
      4: txStart, 5: txEnd, 6: cdsStart, 7: cdsEnd,
      8: exonCount, 9: exonStarts, 10: exonEnds,
      11: score, 12: name2 (gene symbol)
    """
    genes = []
    total = 0
    with gzip.open(REFGENE_GZ, "rt") as f:
        for line in f:
            total += 1
            fields = line.strip().split("\t")
            if len(fields) < 13:
                continue
            chrom = fields[2]
            if chrom != "chr22":
                continue

            name2 = fields[12]  # gene symbol
            transcript = fields[1]
            strand = fields[3]
            tx_start = int(fields[4])
            tx_end = int(fields[5])
            cds_start = int(fields[6])
            cds_end = int(fields[7])
            exon_count = int(fields[8])

            # Parse exon coordinates (comma-separated, trailing comma)
            exon_starts = [int(x) for x in fields[9].rstrip(",").split(",") if x]
            exon_ends = [int(x) for x in fields[10].rstrip(",").split(",") if x]
            exon_lengths = [e - s for s, e in zip(exon_starts, exon_ends)]

            genes.append({
                "gene_name": name2,
                "transcript": transcript,
                "strand": strand,
                "tx_start": tx_start,
                "tx_end": tx_end,
                "tx_len": tx_end - tx_start,
                "cds_start": cds_start,
                "cds_end": cds_end,
                "exon_count": exon_count,
                "exon_starts": exon_starts,
                "exon_ends": exon_ends,
                "exon_lengths": exon_lengths,
            })

    print(f"  Total refGene entries: {total}")
    print(f"  chr22 transcripts: {len(genes)}")

    # Deduplicate by gene name (keep longest transcript per gene)
    by_gene = {}
    for g in genes:
        name = g["gene_name"]
        if name not in by_gene or g["tx_len"] > by_gene[name]["tx_len"]:
            by_gene[name] = g
    unique_genes = sorted(by_gene.values(), key=lambda g: g["tx_start"])
    print(f"  Unique genes (longest transcript): {len(unique_genes)}")

    return unique_genes


def save_tsv(genes, path):
    """Save gene annotations as TSV."""
    with open(path, "w") as f:
        f.write("gene_name\ttranscript\tstrand\ttx_start\ttx_end\ttx_len\t"
                "cds_start\tcds_end\texon_count\texon_starts\texon_ends\texon_lengths\n")
        for g in genes:
            es = ",".join(str(x) for x in g["exon_starts"])
            ee = ",".join(str(x) for x in g["exon_ends"])
            el = ",".join(str(x) for x in g["exon_lengths"])
            f.write(f"{g['gene_name']}\t{g['transcript']}\t{g['strand']}\t"
                    f"{g['tx_start']}\t{g['tx_end']}\t{g['tx_len']}\t"
                    f"{g['cds_start']}\t{g['cds_end']}\t{g['exon_count']}\t"
                    f"{es}\t{ee}\t{el}\n")


def main():
    print("=" * 60)
    print("UCSC refGene: Download and Extract chr22 Gene Annotations")
    print("=" * 60)

    print("\nStep 1: Download refGene.txt.gz...")
    download_refgene()

    print("\nStep 2: Extract chr22 genes...")
    genes = extract_chr22_genes()

    # Summary statistics
    tx_lens = [g["tx_len"] for g in genes]
    exon_lens = [el for g in genes for el in g["exon_lengths"]]
    total_exons = sum(g["exon_count"] for g in genes)

    print(f"\n  Gene body lengths:")
    print(f"    Min: {min(tx_lens):,} bp")
    print(f"    Max: {max(tx_lens):,} bp")
    print(f"    Median: {sorted(tx_lens)[len(tx_lens)//2]:,} bp")

    print(f"\n  Exon statistics:")
    print(f"    Total exons: {total_exons}")
    print(f"    Exon lengths: {min(exon_lens)}-{max(exon_lens)} bp "
          f"(median {sorted(exon_lens)[len(exon_lens)//2]} bp)")

    # Size distribution
    size_bins = [(0, 1000, "<1 Kbp"), (1000, 10000, "1-10 Kbp"),
                 (10000, 100000, "10-100 Kbp"), (100000, 1000000, "100 Kbp-1 Mbp"),
                 (1000000, 10000000, ">1 Mbp")]
    print(f"\n  Gene body size distribution:")
    for lo, hi, label in size_bins:
        count = sum(1 for l in tx_lens if lo <= l < hi)
        if count > 0:
            print(f"    {label}: {count} genes")

    print(f"\nStep 3: Save TSV...")
    save_tsv(genes, OUTPUT_TSV)
    print(f"  Saved: {OUTPUT_TSV}")

    print(f"\n{'=' * 60}")
    print(f"Done. {len(genes)} chr22 genes with {total_exons} exons.")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
