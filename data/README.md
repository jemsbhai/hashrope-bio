# Data directory

This directory holds downloaded datasets for experiments. All files are
gitignored — use the download script to fetch them:

```bash
python scripts/download_data.py --phase 1          # chr22 + HIV-1 (~12 MB)
python scripts/download_data.py --phase 1 --full    # full GRCh38 (~1 GB compressed)
python scripts/download_data.py --phase 2          # MD trajectories
python scripts/download_data.py --phase 3          # PubChem SMILES
```

See `scripts/download_data.py` for exact URLs and checksums.
