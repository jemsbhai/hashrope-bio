"""Extract ClinVar pathogenic SNVs on chr22 for E-G3 real-data validation.

Downloads the ClinVar VCF (GRCh38), filters for:
  - Chromosome 22
  - Single nucleotide variants (SNVs) only
  - Pathogenic or Likely_pathogenic significance
  - REF allele matches chr22.fa at that position (coordinate validation)

Outputs a TSV file for use by Rust/Python E-G3 benchmarks.

Usage:
    python scripts/extract_clinvar_chr22.py --fasta data/chr22.fa --output data/clinvar_chr22_pathogenic_snvs.tsv

Data source:
    ClinVar VCF (GRCh38): https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
    Citation: Landrum MJ et al. Nucleic Acids Res. 2018;46(D1):D1062-D1067. PMID: 29165669
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import os
import sys
import urllib.request
from pathlib import Path

CLINVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
CLINVAR_FILENAME = "clinvar.vcf.gz"

# ClinVar VCF uses bare chromosome names ("22"), UCSC FASTA uses "chr22"
TARGET_CHROM_VCF = "22"  # What ClinVar calls it
TARGET_CHROM_FASTA = "chr22"  # What our FASTA calls it

PATHOGENIC_TERMS = {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}


def download_clinvar(data_dir: Path) -> Path:
    """Download ClinVar VCF if not already present."""
    dest = data_dir / CLINVAR_FILENAME
    if dest.exists():
        size_mb = dest.stat().st_size / 1e6
        print(f"  ClinVar VCF already exists: {dest} ({size_mb:.1f} MB)")
        return dest

    print(f"  Downloading ClinVar VCF from {CLINVAR_URL}")
    print(f"  This is ~173 MB — may take a few minutes...")

    req = urllib.request.Request(CLINVAR_URL, headers={"User-Agent": "hashrope-bio/0.1"})
    with urllib.request.urlopen(req) as response, open(dest, "wb") as out:
        total = response.headers.get("Content-Length")
        downloaded = 0
        block_size = 65536
        while True:
            data = response.read(block_size)
            if not data:
                break
            out.write(data)
            downloaded += len(data)
            if total:
                pct = downloaded / int(total) * 100
                print(f"\r  {downloaded / 1e6:.1f} / {int(total) / 1e6:.1f} MB ({pct:.0f}%)", end="", flush=True)
            else:
                print(f"\r  {downloaded / 1e6:.1f} MB", end="", flush=True)

    print(f"\n  Done: {dest.stat().st_size / 1e6:.1f} MB")
    return dest


def load_chr22_sequence(fasta_path: Path) -> bytes:
    """Load chr22 sequence from FASTA (for REF allele validation)."""
    print(f"  Loading reference sequence from {fasta_path} ...")
    seq_parts = []
    found_header = False
    with open(fasta_path, "r") as f:
        for line in f:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if found_header:
                    break
                found_header = True
                continue
            if found_header:
                seq_parts.append(line.upper())

    seq = "".join(seq_parts).encode("ascii")
    print(f"  Loaded {len(seq):,} bp")
    return seq


def parse_info_field(info: str) -> dict[str, str]:
    """Parse VCF INFO field into a dict."""
    result = {}
    for item in info.split(";"):
        if "=" in item:
            key, val = item.split("=", 1)
            result[key] = val
        else:
            result[item] = ""
    return result


def extract_pathogenic_snvs(
    vcf_path: Path,
    chr22_seq: bytes,
    output_path: Path,
) -> int:
    """Parse ClinVar VCF, extract chr22 pathogenic SNVs, validate, and write TSV."""

    total_variants = 0
    chr22_variants = 0
    chr22_snvs = 0
    pathogenic_snvs = 0
    ref_validated = 0
    ref_mismatches = 0
    written = 0

    print(f"\n  Parsing {vcf_path} ...")

    with gzip.open(vcf_path, "rt") as vcf, open(output_path, "w") as out:
        # Write TSV header
        out.write("# ClinVar pathogenic/likely_pathogenic SNVs on chr22 (GRCh38)\n")
        out.write(f"# Source: {CLINVAR_URL}\n")
        out.write(f"# Extracted by: hashrope-bio extract_clinvar_chr22.py\n")
        out.write("#\n")
        out.write("pos_1based\tpos_0based\tref\talt\tclinvar_id\tclnsig\tgene\tdisease\tref_validated\n")

        for line in vcf:
            if line.startswith("#"):
                continue

            total_variants += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue

            chrom = fields[0]
            if chrom != TARGET_CHROM_VCF:
                continue
            chr22_variants += 1

            pos_1based = int(fields[1])  # VCF is 1-based
            vcf_id = fields[2]
            ref_allele = fields[3].upper()
            alt_allele = fields[4].upper()
            info = parse_info_field(fields[7])

            # Filter: SNVs only (single nucleotide, len(REF)==1, len(ALT)==1)
            if len(ref_allele) != 1 or len(alt_allele) != 1:
                continue
            if ref_allele not in "ACGT" or alt_allele not in "ACGT":
                continue
            chr22_snvs += 1

            # Filter: pathogenic significance
            clnsig = info.get("CLNSIG", "")
            if clnsig not in PATHOGENIC_TERMS:
                continue
            pathogenic_snvs += 1

            # Validate REF allele against chr22.fa
            pos_0based = pos_1based - 1
            ref_ok = False
            if 0 <= pos_0based < len(chr22_seq):
                actual_ref = chr(chr22_seq[pos_0based])
                if actual_ref == ref_allele:
                    ref_ok = True
                    ref_validated += 1
                else:
                    ref_mismatches += 1
                    # Still include it but flag as not validated
            else:
                ref_mismatches += 1

            gene = info.get("GENEINFO", "").split(":")[0] if "GENEINFO" in info else ""
            disease = info.get("CLNDN", "").replace("\t", " ")

            out.write(f"{pos_1based}\t{pos_0based}\t{ref_allele}\t{alt_allele}\t{vcf_id}\t{clnsig}\t{gene}\t{disease}\t{ref_ok}\n")
            written += 1

    print(f"\n  --- Extraction Summary ---")
    print(f"  Total ClinVar variants parsed: {total_variants:,}")
    print(f"  Chromosome 22 variants: {chr22_variants:,}")
    print(f"  Chr22 SNVs (single nucleotide): {chr22_snvs:,}")
    print(f"  Pathogenic/Likely_pathogenic SNVs: {pathogenic_snvs:,}")
    print(f"  REF allele validated against chr22.fa: {ref_validated:,}")
    print(f"  REF allele mismatches: {ref_mismatches:,}")
    print(f"  Written to TSV: {written:,}")
    print(f"  Output: {output_path}")

    if ref_mismatches > 0:
        print(f"\n  WARNING: {ref_mismatches} REF allele mismatches detected.")
        print(f"  These may be due to coordinate system differences or N regions in chr22.")

    return written


def main():
    parser = argparse.ArgumentParser(
        description="Extract ClinVar pathogenic SNVs on chr22 for E-G3 validation"
    )
    parser.add_argument("--fasta", type=Path, required=True,
                        help="Path to chr22.fa (GRCh38)")
    parser.add_argument("--output", type=Path, default=None,
                        help="Output TSV path (default: data/clinvar_chr22_pathogenic_snvs.tsv)")
    parser.add_argument("--skip-download", action="store_true",
                        help="Skip download, use existing ClinVar VCF")
    args = parser.parse_args()

    data_dir = args.fasta.parent
    output_path = args.output or (data_dir / "clinvar_chr22_pathogenic_snvs.tsv")

    print("=== ClinVar chr22 Pathogenic SNV Extraction ===")
    print(f"  FASTA: {args.fasta}")
    print(f"  Output: {output_path}")
    print(f"  Source: {CLINVAR_URL}")

    # Step 1: Download ClinVar VCF
    if not args.skip_download:
        vcf_path = download_clinvar(data_dir)
    else:
        vcf_path = data_dir / CLINVAR_FILENAME
        if not vcf_path.exists():
            print(f"  ERROR: {vcf_path} not found. Remove --skip-download to fetch it.")
            sys.exit(1)

    # Step 2: Load chr22 reference
    chr22_seq = load_chr22_sequence(args.fasta)

    # Step 3: Extract and validate
    count = extract_pathogenic_snvs(vcf_path, chr22_seq, output_path)

    if count == 0:
        print("\n  WARNING: No pathogenic SNVs found. Check ClinVar VCF format.")
        sys.exit(1)

    print(f"\n  SUCCESS: {count} pathogenic SNVs ready for E-G3 real-data validation.")


if __name__ == "__main__":
    main()
