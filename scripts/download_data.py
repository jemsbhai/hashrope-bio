"""Download datasets for hashrope-bio experiments.

Downloads only what's needed for a given experiment phase.
All files go into the data/ directory.

Usage:
    python scripts/download_data.py --phase 1          # Genomics (chr22 only, ~60 MB)
    python scripts/download_data.py --phase 1 --full    # Genomics (full GRCh38, ~3.5 GB)
    python scripts/download_data.py --phase 2          # Proteomics (test trajectories)
    python scripts/download_data.py --phase 3          # Drug discovery (PubChem subset)
"""

from __future__ import annotations

import argparse
import hashlib
import urllib.request
from pathlib import Path

DATA_DIR = Path(__file__).parent.parent / "data"

# --- Dataset registry ---

DATASETS = {
    # Phase 1: Genomics
    "chr22": {
        "url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz",
        "filename": "chr22.fa.gz",
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
    "hiv1_hxb2": {
        "url": "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=K03455.1&db=nuccore&report=fasta&rettype=fasta",
        "filename": "hiv1_hxb2.fa",
        "phase": 1,
        "size_mb": 0.01,
        "description": "HIV-1 HXB2 reference (9.7 Kbp)",
    },

    # Phase 2: Proteomics
    "alanine_dipeptide": {
        "url": "https://figshare.com/ndownloader/files/TODO",  # TODO: find stable URL
        "filename": "alanine_dipeptide.xtc",
        "phase": 2,
        "size_mb": 50,
        "description": "Alanine dipeptide MD trajectory (standard benchmark)",
    },

    # Phase 3: Drug discovery
    # PubChem SMILES are large — download a subset
    "pubchem_1m": {
        "url": "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-SMILES.gz",
        "filename": "CID-SMILES.gz",
        "phase": 3,
        "size_mb": 1500,
        "description": "PubChem canonical SMILES (full, ~110M compounds). Truncate to 1M for experiments.",
    },
}


def download_file(url: str, dest: Path) -> None:
    """Download a file with progress reporting."""
    if dest.exists():
        print(f"  Already exists: {dest}")
        return

    print(f"  Downloading: {url}")
    print(f"  Destination: {dest}")

    # TODO: add proper progress bar (tqdm or urllib reporthook)
    urllib.request.urlretrieve(url, dest)
    print(f"  Done: {dest.stat().st_size / 1e6:.1f} MB")


def main():
    parser = argparse.ArgumentParser(description="Download hashrope-bio datasets")
    parser.add_argument("--phase", type=int, required=True, choices=[1, 2, 3, 4])
    parser.add_argument("--full", action="store_true", help="Include large datasets")
    parser.add_argument("--list", action="store_true", help="List datasets without downloading")
    args = parser.parse_args()

    DATA_DIR.mkdir(exist_ok=True)

    for name, info in DATASETS.items():
        if info["phase"] != args.phase:
            continue
        if info.get("full_only") and not args.full:
            continue

        if args.list:
            print(f"  {name}: {info['description']} ({info['size_mb']} MB)")
        else:
            print(f"\n[{name}] {info['description']}")
            dest = DATA_DIR / info["filename"]
            download_file(info["url"], dest)


if __name__ == "__main__":
    main()
