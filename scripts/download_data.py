"""Download datasets for hashrope-bio experiments.

Downloads only what's needed for a given experiment phase.
All files go into the data/ directory.

Usage:
    python scripts/download_data.py --phase 1          # Genomics (chr22 only, ~60 MB)
    python scripts/download_data.py --phase 1 --full    # Genomics (full GRCh38, ~3.5 GB)
    python scripts/download_data.py --phase 3          # Drug discovery (PubChem subset)
    python scripts/download_data.py --quick             # Just HIV-1 + chr22 (~60 MB)
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import shutil
import sys
import urllib.request
from pathlib import Path

DATA_DIR = Path(__file__).parent.parent / "data"

# --- Dataset registry ---

DATASETS = {
    # Quick / always needed
    "hiv1_hxb2": {
        "url": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=K03455.1&rettype=fasta&retmode=text",
        "filename": "hiv1_hxb2.fa",
        "phase": 0,  # always available
        "size_mb": 0.01,
        "description": "HIV-1 HXB2 reference (9,719 bp) — drug resistance benchmarks",
        "sha256": None,  # will vary by fetch time due to NCBI formatting
    },

    # Phase 1: Genomics
    "chr22": {
        "url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz",
        "filename": "chr22.fa.gz",
        "uncompressed": "chr22.fa",
        "phase": 1,
        "size_mb": 12,
        "description": "GRCh38 chromosome 22 (~51 MB uncompressed)",
    },
    "grch38_full": {
        "url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
        "filename": "hg38.fa.gz",
        "phase": 1,
        "full_only": True,
        "size_mb": 938,
        "description": "Full GRCh38 reference genome (~3.1 GB uncompressed)",
    },

    # Phase 3: Drug discovery
    "pubchem_1m": {
        "url": "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-SMILES.gz",
        "filename": "CID-SMILES.gz",
        "phase": 3,
        "size_mb": 1500,
        "description": "PubChem canonical SMILES (full, ~110M compounds). Truncate to 1M for experiments.",
    },
}


def download_file(url: str, dest: Path, desc: str = "") -> bool:
    """Download a file with progress reporting. Returns True if downloaded."""
    if dest.exists():
        size_mb = dest.stat().st_size / 1e6
        print(f"  Already exists: {dest} ({size_mb:.1f} MB)")
        return False

    print(f"  Downloading: {desc or url}")
    print(f"  Destination: {dest}")

    try:
        req = urllib.request.Request(url, headers={"User-Agent": "hashrope-bio/0.1"})
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
        return True
    except Exception as e:
        print(f"\n  FAILED: {e}")
        if dest.exists():
            dest.unlink()
        return False


def decompress_gz(gz_path: Path, out_path: Path) -> None:
    """Decompress a .gz file if the output doesn't exist."""
    if out_path.exists():
        print(f"  Already decompressed: {out_path}")
        return

    print(f"  Decompressing: {gz_path} -> {out_path}")
    with gzip.open(gz_path, "rb") as f_in, open(out_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    print(f"  Done: {out_path.stat().st_size / 1e6:.1f} MB")


def checksum_file(path: Path) -> str:
    """Compute SHA-256 of a file."""
    sha = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            data = f.read(65536)
            if not data:
                break
            sha.update(data)
    return sha.hexdigest()


def download_phase(phase: int, full: bool = False) -> None:
    """Download all datasets for a given phase."""
    DATA_DIR.mkdir(exist_ok=True)

    for name, info in DATASETS.items():
        if info["phase"] != phase:
            continue
        if info.get("full_only") and not full:
            continue

        print(f"\n[{name}] {info['description']}")
        dest = DATA_DIR / info["filename"]
        downloaded = download_file(info["url"], dest, info["description"])

        # Decompress if needed
        if "uncompressed" in info:
            decompress_gz(dest, DATA_DIR / info["uncompressed"])

        # Log checksum
        final = DATA_DIR / info.get("uncompressed", info["filename"])
        if final.exists():
            sha = checksum_file(final)
            print(f"  SHA-256: {sha}")


def download_quick() -> None:
    """Download minimal datasets for immediate benchmarking."""
    DATA_DIR.mkdir(exist_ok=True)

    # HIV-1 HXB2 — tiny, always needed
    info = DATASETS["hiv1_hxb2"]
    print(f"\n[hiv1_hxb2] {info['description']}")
    dest = DATA_DIR / info["filename"]
    download_file(info["url"], dest, info["description"])
    if dest.exists():
        print(f"  SHA-256: {checksum_file(dest)}")

    # chr22
    info = DATASETS["chr22"]
    print(f"\n[chr22] {info['description']}")
    dest = DATA_DIR / info["filename"]
    download_file(info["url"], dest, info["description"])
    if "uncompressed" in info:
        decompress_gz(dest, DATA_DIR / info["uncompressed"])
    final = DATA_DIR / info.get("uncompressed", info["filename"])
    if final.exists():
        print(f"  SHA-256: {checksum_file(final)}")


def main():
    parser = argparse.ArgumentParser(description="Download hashrope-bio datasets")
    parser.add_argument("--phase", type=int, choices=[0, 1, 2, 3, 4])
    parser.add_argument("--full", action="store_true", help="Include large datasets")
    parser.add_argument("--quick", action="store_true", help="HIV-1 + chr22 only (~60 MB)")
    parser.add_argument("--list", action="store_true", help="List datasets without downloading")
    args = parser.parse_args()

    if args.list:
        for name, info in DATASETS.items():
            full_tag = " [--full only]" if info.get("full_only") else ""
            print(f"  Phase {info['phase']}: {name} — {info['description']} ({info['size_mb']} MB){full_tag}")
        return

    if args.quick:
        download_quick()
        return

    if args.phase is not None:
        download_phase(args.phase, args.full)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
