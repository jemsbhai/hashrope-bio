#!/usr/bin/env python3
"""Download real HIV-1 RT sequences from NCBI GenBank and annotate
with Stanford HIVDB Sierra for E-D2 real-data validation.

Pipeline:
  1. Search NCBI Entrez for HIV-1 RT nucleotide sequences
  2. Download as FASTA
  3. Submit to Stanford HIVDB Sierra GraphQL API
  4. Extract mutations at our 23 resistance panel positions
  5. Save: FASTA + TSV ground truth

Requirements: requests (pip install requests)
No API key needed for NCBI (with rate limiting) or HIVDB Sierra.

Output files (in data/):
  - hivdb_rt_sequences.fa       — raw FASTA sequences
  - hivdb_rt_sierra_results.json — full Sierra JSON response
  - hivdb_rt_ground_truth.tsv   — per-sequence mutations at panel positions
"""

import json
import sys
import time
import urllib.request
import urllib.parse
from pathlib import Path

try:
    import requests
except ImportError:
    print("ERROR: 'requests' package required. Install with: pip install requests")
    sys.exit(1)

# --- Configuration ---
DATA_DIR = Path(__file__).parent.parent / "data"
DATA_DIR.mkdir(exist_ok=True)

# NCBI Entrez settings
ENTREZ_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
# Search for HIV-1 subtype B RT sequences from drug-resistance studies
# These are clinical genotyping sequences, typically 500-1700 bp covering RT
ENTREZ_QUERY = (
    '"Human immunodeficiency virus 1"[Organism] AND pol[Gene Name] '
    'AND 600:1800[Sequence Length] AND "drug resistant"[All Fields]'
)
MAX_SEQUENCES = 100  # download up to this many
ENTREZ_RETMAX = 200  # search result cap

# Stanford HIVDB Sierra GraphQL
SIERRA_URL = "https://hivdb.stanford.edu/graphql"

# Our 23 resistance panel positions (amino acid, 1-based)
PANEL_POSITIONS = {
    # NNRTI (12 positions)
    100, 101, 103, 106, 108, 138, 181, 188, 190, 225, 227, 230,
    # NRTI (11 positions)
    41, 65, 67, 70, 74, 115, 151, 184, 210, 215, 219,
}


def entrez_search(query: str, retmax: int = 200) -> list[str]:
    """Search NCBI and return list of GI/accession IDs."""
    params = {
        "db": "nucleotide",
        "term": query,
        "retmax": retmax,
        "retmode": "json",
        "usehistory": "y",
    }
    url = f"{ENTREZ_BASE}/esearch.fcgi?{urllib.parse.urlencode(params)}"
    print(f"Searching NCBI: {query[:80]}...")
    with urllib.request.urlopen(url) as resp:
        data = json.loads(resp.read())
    result = data["esearchresult"]
    ids = result["idlist"]
    total = int(result["count"])
    print(f"  Found {total} total, retrieved {len(ids)} IDs")
    return ids


def entrez_fetch_fasta(ids: list[str], batch_size: int = 50) -> str:
    """Download sequences as FASTA from NCBI in batches."""
    all_fasta = []
    for i in range(0, len(ids), batch_size):
        batch = ids[i:i + batch_size]
        params = {
            "db": "nucleotide",
            "id": ",".join(batch),
            "rettype": "fasta",
            "retmode": "text",
        }
        url = f"{ENTREZ_BASE}/efetch.fcgi?{urllib.parse.urlencode(params)}"
        print(f"  Fetching batch {i // batch_size + 1} ({len(batch)} sequences)...")
        with urllib.request.urlopen(url) as resp:
            all_fasta.append(resp.read().decode("utf-8"))
        time.sleep(0.5)  # NCBI rate limit
    return "\n".join(all_fasta)


def parse_fasta(text: str) -> list[tuple[str, str]]:
    """Parse FASTA text into list of (header, sequence) tuples."""
    entries = []
    header = None
    seq_lines = []
    for line in text.strip().split("\n"):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                entries.append((header, "".join(seq_lines)))
            header = line[1:].strip()
            seq_lines = []
        else:
            seq_lines.append(line)
    if header is not None:
        entries.append((header, "".join(seq_lines)))
    return entries


def sierra_analyze(sequences: list[tuple[str, str]], batch_size: int = 20) -> list[dict]:
    """Submit sequences to Stanford HIVDB Sierra GraphQL API.

    Returns list of SequenceAnalysis objects.
    """
    # GraphQL query — request aligned gene sequences and mutations
    query = """
    query SequenceAnalysis($sequences: [UnalignedSequenceInput]!) {
      viewer {
        sequenceAnalysis(sequences: $sequences) {
          inputSequence {
            header
          }
          bestMatchingSubtype {
            displayWithoutDistance
          }
          validationResults {
            level
            message
          }
          alignedGeneSequences {
            gene {
              name
            }
            firstAA
            lastAA
            firstNA
            lastNA
            alignedNAs
            mutations {
              position
              AAs
              isInsertion
              isDeletion
              primaryType
              text
            }
          }
        }
      }
    }
    """

    all_results = []
    for i in range(0, len(sequences), batch_size):
        batch = sequences[i:i + batch_size]
        variables = {
            "sequences": [
                {"header": hdr, "sequence": seq}
                for hdr, seq in batch
            ]
        }
        print(f"  Submitting batch {i // batch_size + 1} to Sierra "
              f"({len(batch)} sequences)...")

        resp = requests.post(
            SIERRA_URL,
            json={"query": query, "variables": variables},
            headers={"Content-Type": "application/json"},
            timeout=120,
        )
        resp.raise_for_status()
        data = resp.json()

        if "errors" in data:
            print(f"  WARNING: Sierra returned errors: {data['errors']}")
            continue

        results = data["data"]["viewer"]["sequenceAnalysis"]
        all_results.extend(results)
        print(f"    Got {len(results)} results")
        time.sleep(1)  # be polite

    return all_results


def extract_panel_mutations(sierra_results: list[dict]) -> list[dict]:
    """Extract mutations at our 23 panel positions from Sierra results.

    Returns list of dicts with:
      - header: sequence header
      - subtype: best matching subtype
      - rt_first_aa, rt_last_aa: coverage range
      - aligned_rt_nas: aligned RT nucleotide sequence
      - panel_mutations: dict of position -> mutation text
      - all_rt_mutations: list of all RT mutation texts
      - warnings: list of validation warnings
    """
    records = []
    for result in sierra_results:
        header = result["inputSequence"]["header"]
        subtype = result["bestMatchingSubtype"]["displayWithoutDistance"] if result.get("bestMatchingSubtype") else "?"

        warnings = [
            v["message"] for v in result.get("validationResults", [])
            if v["level"] in ("WARNING", "CRITICAL")
        ]

        # Find RT gene alignment
        rt_data = None
        for gene_seq in result.get("alignedGeneSequences", []):
            if gene_seq["gene"]["name"] == "RT":
                rt_data = gene_seq
                break

        if rt_data is None:
            print(f"  SKIP: {header} — no RT alignment")
            continue

        # Extract mutations at panel positions
        panel_muts = {}
        all_muts = []
        for mut in rt_data.get("mutations", []):
            pos = mut["position"]
            text = mut["text"]
            all_muts.append(text)
            if pos in PANEL_POSITIONS:
                panel_muts[pos] = text

        records.append({
            "header": header,
            "subtype": subtype,
            "rt_first_aa": rt_data["firstAA"],
            "rt_last_aa": rt_data["lastAA"],
            "aligned_rt_nas": rt_data["alignedNAs"],
            "panel_mutations": panel_muts,
            "all_rt_mutations": all_muts,
            "warnings": warnings,
        })

    return records


def save_ground_truth_tsv(records: list[dict], path: Path):
    """Save ground truth TSV for Rust benchmark consumption."""
    with open(path, "w") as f:
        # Header
        f.write("header\tsubtype\trt_first_aa\trt_last_aa\t"
                "panel_mutation_positions\tpanel_mutation_texts\t"
                "all_rt_mutations\twarnings\n")
        for rec in records:
            panel_pos = ",".join(str(p) for p in sorted(rec["panel_mutations"].keys()))
            panel_txt = ",".join(rec["panel_mutations"][p] for p in sorted(rec["panel_mutations"].keys()))
            all_muts = ",".join(rec["all_rt_mutations"])
            warns = ";".join(rec["warnings"])
            f.write(f"{rec['header']}\t{rec['subtype']}\t"
                    f"{rec['rt_first_aa']}\t{rec['rt_last_aa']}\t"
                    f"{panel_pos}\t{panel_txt}\t{all_muts}\t{warns}\n")


def save_aligned_fasta(records: list[dict], path: Path):
    """Save Sierra-aligned RT nucleotide sequences as FASTA.

    These are aligned to HXB2 RT reading frame — codon positions
    match our panel positions directly.
    """
    with open(path, "w") as f:
        for rec in records:
            # Use a clean header (first token)
            hdr = rec["header"].split()[0]
            seq = rec["aligned_rt_nas"]
            f.write(f">{hdr} subtype={rec['subtype']} "
                    f"RT:{rec['rt_first_aa']}-{rec['rt_last_aa']} "
                    f"panel_muts={len(rec['panel_mutations'])}\n")
            # Wrap at 80 chars
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")


def main():
    print("=" * 60)
    print("E-D2 Real Data: Download HIV-1 RT + Sierra Annotation")
    print("=" * 60)

    # Step 1: Search NCBI
    ids = entrez_search(ENTREZ_QUERY, retmax=ENTREZ_RETMAX)
    if not ids:
        print("ERROR: No sequences found. Try different search terms.")
        sys.exit(1)

    # Limit to MAX_SEQUENCES
    ids = ids[:MAX_SEQUENCES]
    print(f"Using {len(ids)} sequence IDs")

    # Step 2: Download FASTA
    print("\nDownloading FASTA from NCBI...")
    fasta_text = entrez_fetch_fasta(ids)
    sequences = parse_fasta(fasta_text)
    print(f"  Parsed {len(sequences)} sequences")

    # Save raw FASTA
    raw_fasta_path = DATA_DIR / "hivdb_rt_sequences.fa"
    with open(raw_fasta_path, "w") as f:
        f.write(fasta_text)
    print(f"  Saved raw FASTA: {raw_fasta_path}")

    # Step 3: Submit to Sierra
    print("\nSubmitting to Stanford HIVDB Sierra...")
    sierra_results = sierra_analyze(sequences)
    print(f"  Got {len(sierra_results)} Sierra results")

    # Save full Sierra JSON
    sierra_json_path = DATA_DIR / "hivdb_rt_sierra_results.json"
    with open(sierra_json_path, "w") as f:
        json.dump(sierra_results, f, indent=2)
    print(f"  Saved Sierra JSON: {sierra_json_path}")

    # Step 4: Extract panel mutations
    print("\nExtracting panel mutations...")
    records = extract_panel_mutations(sierra_results)
    print(f"  {len(records)} sequences with RT alignment")

    # Count sequences with panel mutations
    with_muts = sum(1 for r in records if r["panel_mutations"])
    total_panel_muts = sum(len(r["panel_mutations"]) for r in records)
    print(f"  {with_muts} sequences have >= 1 panel mutation")
    print(f"  {total_panel_muts} total panel mutations across all sequences")

    # Mutation frequency at each panel position
    pos_counts = {}
    for rec in records:
        for pos in rec["panel_mutations"]:
            pos_counts[pos] = pos_counts.get(pos, 0) + 1
    if pos_counts:
        print("\n  Panel position mutation frequency:")
        for pos in sorted(pos_counts.keys()):
            print(f"    RT:{pos} — {pos_counts[pos]} sequences")

    # Step 5: Save outputs
    tsv_path = DATA_DIR / "hivdb_rt_ground_truth.tsv"
    save_ground_truth_tsv(records, tsv_path)
    print(f"\n  Saved ground truth TSV: {tsv_path}")

    aligned_fasta_path = DATA_DIR / "hivdb_rt_aligned.fa"
    save_aligned_fasta(records, aligned_fasta_path)
    print(f"  Saved aligned FASTA: {aligned_fasta_path}")

    # Summary
    print("\n" + "=" * 60)
    print("Summary:")
    print(f"  Sequences downloaded: {len(sequences)}")
    print(f"  With RT alignment: {len(records)}")
    print(f"  With panel mutations: {with_muts}")
    print(f"  Total panel mutations: {total_panel_muts}")
    print(f"  Output files:")
    print(f"    {raw_fasta_path}")
    print(f"    {sierra_json_path}")
    print(f"    {aligned_fasta_path}")
    print(f"    {tsv_path}")
    print("=" * 60)


if __name__ == "__main__":
    main()
