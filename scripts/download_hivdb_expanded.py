#!/usr/bin/env python3
"""Expanded E-D2: Download 2,000+ HIV-1 RT sequences from NCBI GenBank,
annotate with Stanford HIVDB Sierra (including drug resistance scoring).

This is the scaled-up version of download_hivdb_sequences.py.
Addresses three gaps:
  Gap 1: Scale from 100 → 2,000+ sequences
  Gap 2: Add per-drug resistance scoring (drugResistance from Sierra)
  Gap 3: (RHIVDB treatment history — deferred to separate script)

Pipeline:
  1. Search NCBI Entrez for HIV-1 RT nucleotide sequences
  2. Download as FASTA in batches
  3. Submit to Stanford HIVDB Sierra GraphQL API (batches of 20)
     — includes drugResistance scoring (per-drug level + score)
  4. Extract mutations at 23 resistance panel positions
  5. Save: FASTA + TSV ground truth + drug scoring TSV + full Sierra JSON

Requirements: requests (pip install requests)
No API key needed for NCBI (with rate limiting) or HIVDB Sierra.

Output files (in data/):
  - hivdb_rt_expanded_sequences.fa         — raw FASTA sequences
  - hivdb_rt_expanded_aligned.fa           — Sierra-aligned RT nucleotides
  - hivdb_rt_expanded_sierra_results.json  — full Sierra JSON response
  - hivdb_rt_expanded_ground_truth.tsv     — per-sequence mutations at panel positions
  - hivdb_rt_expanded_drug_scores.tsv      — per-sequence per-drug resistance levels
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
ENTREZ_QUERY = (
    '"Human immunodeficiency virus 1"[Organism] AND pol[Gene Name] '
    'AND 600:1800[Sequence Length] AND "drug resistant"[All Fields]'
)
MAX_SEQUENCES = 2000
ENTREZ_RETMAX = 2500  # overshoot to account for duplicates/failures

# Stanford HIVDB Sierra GraphQL
SIERRA_URL = "https://hivdb.stanford.edu/graphql"
SIERRA_BATCH_SIZE = 20
SIERRA_DELAY = 1.5  # seconds between batches (polite for 100+ batches)

# Our 23 resistance panel positions (amino acid, 1-based)
PANEL_POSITIONS = {
    # NNRTI (12 positions)
    100, 101, 103, 106, 108, 138, 181, 188, 190, 225, 227, 230,
    # NRTI (11 positions)
    41, 65, 67, 70, 74, 115, 151, 184, 210, 215, 219,
}

# Sierra GraphQL query — expanded with drugResistance scoring
SIERRA_QUERY = """
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
      drugResistance {
        gene {
          name
        }
        drugScores {
          drug {
            displayAbbr
            fullName
          }
          score
          level
          text
          partialScores {
            mutations {
              text
            }
            score
          }
        }
      }
    }
  }
}
"""


def entrez_search(query: str, retmax: int = 2500) -> list[str]:
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


def entrez_fetch_fasta(ids: list[str], batch_size: int = 100) -> str:
    """Download sequences as FASTA from NCBI in batches."""
    all_fasta = []
    total_batches = (len(ids) + batch_size - 1) // batch_size
    for i in range(0, len(ids), batch_size):
        batch = ids[i:i + batch_size]
        batch_num = i // batch_size + 1
        params = {
            "db": "nucleotide",
            "id": ",".join(batch),
            "rettype": "fasta",
            "retmode": "text",
        }
        url = f"{ENTREZ_BASE}/efetch.fcgi?{urllib.parse.urlencode(params)}"
        print(f"  Fetching FASTA batch {batch_num}/{total_batches} "
              f"({len(batch)} sequences)...")
        retries = 3
        for attempt in range(retries):
            try:
                with urllib.request.urlopen(url, timeout=60) as resp:
                    all_fasta.append(resp.read().decode("utf-8"))
                break
            except Exception as e:
                if attempt < retries - 1:
                    print(f"    Retry {attempt + 1}/{retries}: {e}")
                    time.sleep(2 ** attempt)
                else:
                    print(f"    FAILED after {retries} attempts: {e}")
                    raise
        time.sleep(0.4)  # NCBI rate limit: 3 req/sec without API key
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


def sierra_analyze(sequences: list[tuple[str, str]]) -> list[dict]:
    """Submit sequences to Stanford HIVDB Sierra GraphQL API.

    Returns list of SequenceAnalysis objects (includes drugResistance).
    """
    all_results = []
    total_batches = (len(sequences) + SIERRA_BATCH_SIZE - 1) // SIERRA_BATCH_SIZE
    failed_batches = 0

    for i in range(0, len(sequences), SIERRA_BATCH_SIZE):
        batch = sequences[i:i + SIERRA_BATCH_SIZE]
        batch_num = i // SIERRA_BATCH_SIZE + 1
        variables = {
            "sequences": [
                {"header": hdr, "sequence": seq}
                for hdr, seq in batch
            ]
        }
        print(f"  Sierra batch {batch_num}/{total_batches} "
              f"({len(batch)} seqs, {len(all_results)} done)...", end="", flush=True)

        retries = 3
        success = False
        for attempt in range(retries):
            try:
                resp = requests.post(
                    SIERRA_URL,
                    json={"query": SIERRA_QUERY, "variables": variables},
                    headers={"Content-Type": "application/json"},
                    timeout=180,
                )
                resp.raise_for_status()
                data = resp.json()

                if "errors" in data:
                    print(f" ERRORS: {data['errors'][:200]}")
                    failed_batches += 1
                    break

                results = data["data"]["viewer"]["sequenceAnalysis"]
                all_results.extend(results)
                print(f" OK ({len(results)})")
                success = True
                break
            except requests.exceptions.Timeout:
                print(f" TIMEOUT (attempt {attempt + 1})")
                if attempt < retries - 1:
                    time.sleep(5)
            except Exception as e:
                print(f" ERROR: {e} (attempt {attempt + 1})")
                if attempt < retries - 1:
                    time.sleep(3)

        if not success and attempt == retries - 1:
            print(f"  SKIPPING batch {batch_num} after {retries} failures")
            failed_batches += 1

        time.sleep(SIERRA_DELAY)

    print(f"\n  Sierra complete: {len(all_results)} results, "
          f"{failed_batches} failed batches")
    return all_results


def extract_panel_mutations(sierra_results: list[dict]) -> list[dict]:
    """Extract mutations at 23 panel positions + drug resistance scoring.

    Returns list of dicts with:
      - header, subtype, rt_first_aa, rt_last_aa, aligned_rt_nas
      - panel_mutations: {position: text}
      - all_rt_mutations: [text]
      - drug_scores: [{drug, score, level, text}] for RT drugs
      - warnings: [str]
    """
    records = []
    for result in sierra_results:
        header = result["inputSequence"]["header"]
        subtype_obj = result.get("bestMatchingSubtype")
        subtype = subtype_obj["displayWithoutDistance"] if subtype_obj else "?"

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
            continue  # skip sequences without RT alignment

        # Extract mutations at panel positions
        panel_muts = {}
        all_muts = []
        for mut in rt_data.get("mutations", []):
            pos = mut["position"]
            text = mut["text"]
            all_muts.append(text)
            if pos in PANEL_POSITIONS:
                panel_muts[pos] = text

        # Extract drug resistance scores for RT
        drug_scores = []
        for dr in result.get("drugResistance", []):
            if dr["gene"]["name"] != "RT":
                continue
            for ds in dr.get("drugScores", []):
                drug_scores.append({
                    "drug_abbr": ds["drug"]["displayAbbr"],
                    "drug_name": ds["drug"]["fullName"],
                    "score": ds["score"],
                    "level": ds["level"],
                    "text": ds["text"],
                    "contributing_mutations": [
                        {
                            "mutations": [m["text"] for m in ps.get("mutations", [])],
                            "score": ps["score"],
                        }
                        for ps in ds.get("partialScores", [])
                        if ps["score"] != 0
                    ],
                })

        records.append({
            "header": header,
            "subtype": subtype,
            "rt_first_aa": rt_data["firstAA"],
            "rt_last_aa": rt_data["lastAA"],
            "aligned_rt_nas": rt_data["alignedNAs"],
            "panel_mutations": panel_muts,
            "all_rt_mutations": all_muts,
            "drug_scores": drug_scores,
            "warnings": warnings,
        })

    return records


def save_ground_truth_tsv(records: list[dict], path: Path):
    """Save ground truth TSV (same format as original for benchmark compatibility)."""
    with open(path, "w") as f:
        f.write("header\tsubtype\trt_first_aa\trt_last_aa\t"
                "panel_mutation_positions\tpanel_mutation_texts\t"
                "all_rt_mutations\twarnings\n")
        for rec in records:
            panel_pos = ",".join(str(p) for p in sorted(rec["panel_mutations"].keys()))
            panel_txt = ",".join(rec["panel_mutations"][p]
                                 for p in sorted(rec["panel_mutations"].keys()))
            all_muts = ",".join(rec["all_rt_mutations"])
            warns = ";".join(rec["warnings"])
            f.write(f"{rec['header']}\t{rec['subtype']}\t"
                    f"{rec['rt_first_aa']}\t{rec['rt_last_aa']}\t"
                    f"{panel_pos}\t{panel_txt}\t{all_muts}\t{warns}\n")


def save_drug_scores_tsv(records: list[dict], path: Path):
    """Save per-sequence per-drug resistance scoring TSV.

    One row per (sequence, drug) pair.
    """
    with open(path, "w") as f:
        f.write("header\tsubtype\tdrug_abbr\tdrug_name\tscore\tlevel\t"
                "text\tcontributing_mutations\n")
        for rec in records:
            accession = rec["header"].split()[0]
            for ds in rec["drug_scores"]:
                # Format contributing mutations as semicolon-separated
                contrib = "; ".join(
                    f"{'+'.join(ps['mutations'])}={ps['score']}"
                    for ps in ds["contributing_mutations"]
                )
                f.write(f"{accession}\t{rec['subtype']}\t"
                        f"{ds['drug_abbr']}\t{ds['drug_name']}\t"
                        f"{ds['score']}\t{ds['level']}\t"
                        f"{ds['text']}\t{contrib}\n")


def save_aligned_fasta(records: list[dict], path: Path):
    """Save Sierra-aligned RT nucleotide sequences as FASTA."""
    with open(path, "w") as f:
        for rec in records:
            hdr = rec["header"].split()[0]
            seq = rec["aligned_rt_nas"]
            f.write(f">{hdr} subtype={rec['subtype']} "
                    f"RT:{rec['rt_first_aa']}-{rec['rt_last_aa']} "
                    f"panel_muts={len(rec['panel_mutations'])}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")


def summarize_drug_resistance(records: list[dict]):
    """Print summary of drug resistance levels across all sequences."""
    level_counts = {}  # drug_abbr → {level → count}
    for rec in records:
        for ds in rec["drug_scores"]:
            abbr = ds["drug_abbr"]
            level = ds["level"]
            if abbr not in level_counts:
                level_counts[abbr] = {}
            level_counts[abbr][level] = level_counts[abbr].get(level, 0) + 1

    if not level_counts:
        print("  (no drug resistance data)")
        return

    # Sort by total non-susceptible
    drugs_sorted = sorted(
        level_counts.keys(),
        key=lambda d: sum(v for k, v in level_counts[d].items()
                          if k != 1),  # level 1 = Susceptible
        reverse=True,
    )

    print(f"\n  Drug resistance summary ({len(records)} sequences):")
    print(f"  {'Drug':<8} {'Suscept':>8} {'Low':>6} {'Interm':>8} {'High':>6}")
    print(f"  {'----':<8} {'-------':>8} {'---':>6} {'------':>8} {'----':>6}")
    for drug in drugs_sorted:
        levels = level_counts[drug]
        s = levels.get(1, 0)
        l = levels.get(2, 0)
        m = levels.get(3, 0)
        h = levels.get(4, 0)
        # Sierra levels: 1=Susceptible, 2=Potential Low-Level, 3=Low-Level,
        # 4=Intermediate, 5=High-Level
        # Actually let's just use the text labels
        pass

    # Redo with text labels from the data
    level_text_counts = {}
    for rec in records:
        for ds in rec["drug_scores"]:
            abbr = ds["drug_abbr"]
            text = ds["text"]
            if abbr not in level_text_counts:
                level_text_counts[abbr] = {}
            level_text_counts[abbr][text] = level_text_counts[abbr].get(text, 0) + 1

    print(f"\n  Drug resistance levels ({len(records)} sequences):")
    for drug in sorted(level_text_counts.keys()):
        levels = level_text_counts[drug]
        parts = ", ".join(f"{lbl}: {cnt}" for lbl, cnt in
                          sorted(levels.items(), key=lambda x: -x[1]))
        print(f"    {drug:<8} {parts}")


def main():
    print("=" * 70)
    print("E-D2 Expanded: 2,000+ HIV-1 RT Sequences + Drug Resistance Scoring")
    print("=" * 70)

    # Step 1: Search NCBI
    ids = entrez_search(ENTREZ_QUERY, retmax=ENTREZ_RETMAX)
    if not ids:
        print("ERROR: No sequences found.")
        sys.exit(1)

    # Limit to MAX_SEQUENCES
    ids = ids[:MAX_SEQUENCES]
    print(f"\nUsing {len(ids)} sequence IDs (of {ENTREZ_RETMAX} searched)")

    # Step 2: Download FASTA
    print("\nDownloading FASTA from NCBI...")
    fasta_text = entrez_fetch_fasta(ids, batch_size=100)
    sequences = parse_fasta(fasta_text)
    print(f"  Parsed {len(sequences)} sequences")

    # Deduplicate by accession (first token of header)
    seen = set()
    unique_sequences = []
    for hdr, seq in sequences:
        acc = hdr.split()[0]
        if acc not in seen:
            seen.add(acc)
            unique_sequences.append((hdr, seq))
    if len(unique_sequences) < len(sequences):
        print(f"  Deduplicated: {len(sequences)} → {len(unique_sequences)}")
    sequences = unique_sequences

    # Save raw FASTA
    raw_fasta_path = DATA_DIR / "hivdb_rt_expanded_sequences.fa"
    with open(raw_fasta_path, "w") as f:
        f.write(fasta_text)
    print(f"  Saved raw FASTA: {raw_fasta_path}")

    # Step 3: Submit to Sierra (with drug resistance scoring)
    print(f"\nSubmitting {len(sequences)} sequences to Stanford HIVDB Sierra...")
    print(f"  Batch size: {SIERRA_BATCH_SIZE}, delay: {SIERRA_DELAY}s")
    est_time = len(sequences) / SIERRA_BATCH_SIZE * SIERRA_DELAY
    print(f"  Estimated time: ~{est_time / 60:.1f} minutes (network dependent)")

    t0 = time.time()
    sierra_results = sierra_analyze(sequences)
    elapsed = time.time() - t0
    print(f"  Sierra done in {elapsed:.1f}s ({elapsed / 60:.1f} min)")

    # Save full Sierra JSON
    sierra_json_path = DATA_DIR / "hivdb_rt_expanded_sierra_results.json"
    with open(sierra_json_path, "w") as f:
        json.dump(sierra_results, f, indent=2)
    print(f"  Saved Sierra JSON: {sierra_json_path}")

    # Step 4: Extract panel mutations + drug scores
    print("\nExtracting panel mutations and drug resistance scores...")
    records = extract_panel_mutations(sierra_results)
    print(f"  {len(records)} sequences with RT alignment")

    with_muts = sum(1 for r in records if r["panel_mutations"])
    total_panel_muts = sum(len(r["panel_mutations"]) for r in records)
    with_drug_scores = sum(1 for r in records if r["drug_scores"])
    print(f"  {with_muts} sequences have >= 1 panel mutation")
    print(f"  {total_panel_muts} total panel mutations across all sequences")
    print(f"  {with_drug_scores} sequences have drug resistance scores")

    # Panel position frequency
    pos_counts = {}
    for rec in records:
        for pos in rec["panel_mutations"]:
            pos_counts[pos] = pos_counts.get(pos, 0) + 1
    if pos_counts:
        print("\n  Panel position mutation frequency:")
        for pos in sorted(pos_counts.keys()):
            pct = pos_counts[pos] / len(records) * 100
            print(f"    RT:{pos:>3d} — {pos_counts[pos]:>5d} sequences ({pct:.1f}%)")

    # Drug resistance summary
    summarize_drug_resistance(records)

    # Step 5: Save outputs
    tsv_path = DATA_DIR / "hivdb_rt_expanded_ground_truth.tsv"
    save_ground_truth_tsv(records, tsv_path)
    print(f"\n  Saved ground truth TSV: {tsv_path}")

    aligned_fasta_path = DATA_DIR / "hivdb_rt_expanded_aligned.fa"
    save_aligned_fasta(records, aligned_fasta_path)
    print(f"  Saved aligned FASTA: {aligned_fasta_path}")

    drug_scores_path = DATA_DIR / "hivdb_rt_expanded_drug_scores.tsv"
    save_drug_scores_tsv(records, drug_scores_path)
    print(f"  Saved drug scores TSV: {drug_scores_path}")

    # Summary
    print(f"\n{'=' * 70}")
    print("Summary:")
    print(f"  Sequences downloaded:  {len(sequences)}")
    print(f"  With RT alignment:     {len(records)}")
    print(f"  With panel mutations:  {with_muts}")
    print(f"  Total panel mutations: {total_panel_muts}")
    print(f"  With drug scores:      {with_drug_scores}")
    print(f"  Sierra elapsed:        {elapsed:.1f}s")
    print(f"\n  Output files:")
    print(f"    {raw_fasta_path}")
    print(f"    {sierra_json_path}")
    print(f"    {aligned_fasta_path}")
    print(f"    {tsv_path}")
    print(f"    {drug_scores_path}")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
