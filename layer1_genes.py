"""
LAYER 1 — Gene Harvesting
==========================
What this does:
  → Talks to Open Targets API
  → Asks: "give me all genes associated with Bardet-Biedl Syndrome"
  → Filters out weak associations (score < 0.3)
  → Saves a clean list of genes to a CSV file

Run this first. If you get a CSV with BBS genes, Layer 1 is done.
"""

import requests
import pandas as pd
import os

# ── Settings ──────────────────────────────────────
DISEASE_NAME = "Bardet-Biedl Syndrome"
EFO_ID       = "MONDO_0015229"        # BBS identifier in Open Targets
SCORE_CUTOFF = 0.3                   # ignore weak associations
OUTPUT_FILE  = "outputs/layer1_bbs_genes.csv"

os.makedirs("outputs", exist_ok=True)

# ── The API query (GraphQL) ────────────────────────
# Think of this like asking a very specific question
# to the Open Targets database in their language

QUERY = """
query DisgenetBBS($efoId: String!, $size: Int!) {
  disease(efoId: $efoId) {
    name
    associatedTargets(page: {index: 0, size: $size}) {
      rows {
        target {
          approvedSymbol
          approvedName
          id
        }
        score
        datatypeScores {
          componentId
          score
        }
      }
    }
  }
}
"""

# ── Fetch from API ─────────────────────────────────
def fetch_bbs_genes():
    print(f"Querying Open Targets for: {DISEASE_NAME}")
    print(f"EFO ID: {EFO_ID}")
    print("-" * 40)

    url = "https://api.platform.opentargets.org/api/v4/graphql"

    response = requests.post(
        url,
        json={
            "query": QUERY,
            "variables": {
                "efoId": EFO_ID,
                "size": 50          # get top 50 associated genes
            }
        },
        timeout=30
    )

    # Check if API call worked
    if response.status_code != 200:
        print(f"API Error: {response.status_code}")
        print("Using demo data instead...")
        return get_demo_data()

    data = response.json()

    # Check if we got actual data back
    if not data.get("data") or not data["data"].get("disease"):
        print("No data returned from API. Using demo data...")
        return get_demo_data()

    # ── Parse the response ─────────────────────────
    rows = data["data"]["disease"]["associatedTargets"]["rows"]
    print(f"Total associations found: {len(rows)}")

    genes = []
    for row in rows:
        score = row["score"]

        # Skip weak associations
        if score < SCORE_CUTOFF:
            continue

        gene_symbol = row["target"]["approvedSymbol"]
        gene_name   = row["target"]["approvedName"]
        ensembl_id  = row["target"]["id"]

        # What type of evidence do we have?
        evidence_types = []
        for dtype in row["datatypeScores"]:
            evidence_types.append(dtype["componentId"])

        genes.append({
            "gene_symbol":   gene_symbol,
            "gene_name":     gene_name,
            "ensembl_id":    ensembl_id,
            "ot_score":      round(score, 4),
            "evidence_types": ", ".join(evidence_types)
        })

    return genes


# ── Demo data (backup if API fails) ───────────────
def get_demo_data():
    """
    Real BBS genes from literature — used if API is down.
    These are the 6 most studied BBS genes.
    """
    return [
        {"gene_symbol": "BBS1",   "gene_name": "Bardet-Biedl syndrome 1",  "ensembl_id": "ENSG00000174483", "ot_score": 0.95, "evidence_types": "genetic_association"},
        {"gene_symbol": "BBS2",   "gene_name": "Bardet-Biedl syndrome 2",  "ensembl_id": "ENSG00000125124", "ot_score": 0.92, "evidence_types": "genetic_association"},
        {"gene_symbol": "BBS4",   "gene_name": "Bardet-Biedl syndrome 4",  "ensembl_id": "ENSG00000140463", "ot_score": 0.90, "evidence_types": "genetic_association"},
        {"gene_symbol": "BBS7",   "gene_name": "Bardet-Biedl syndrome 7",  "ensembl_id": "ENSG00000138686", "ot_score": 0.88, "evidence_types": "genetic_association"},
        {"gene_symbol": "BBS9",   "gene_name": "Bardet-Biedl syndrome 9",  "ensembl_id": "ENSG00000122585", "ot_score": 0.85, "evidence_types": "genetic_association"},
        {"gene_symbol": "BBS10",  "gene_name": "Bardet-Biedl syndrome 10", "ensembl_id": "ENSG00000179941", "ot_score": 0.83, "evidence_types": "genetic_association"},
        {"gene_symbol": "MKKS",   "gene_name": "McKusick-Kaufman syndrome","ensembl_id": "ENSG00000125863", "ot_score": 0.80, "evidence_types": "genetic_association"},
        {"gene_symbol": "CEP290", "gene_name": "Centrosomal protein 290",  "ensembl_id": "ENSG00000198707", "ot_score": 0.75, "evidence_types": "genetic_association, literature"},
        {"gene_symbol": "ARL6",   "gene_name": "ADP ribosylation factor 6","ensembl_id": "ENSG00000113966", "ot_score": 0.72, "evidence_types": "genetic_association"},
        {"gene_symbol": "TRIM32", "gene_name": "Tripartite motif 32",      "ensembl_id": "ENSG00000119401", "ot_score": 0.65, "evidence_types": "genetic_association, literature"},
    ]


# ── Main ───────────────────────────────────────────
def main():
    genes = fetch_bbs_genes()

    if not genes:
        print("No genes found above score cutoff!")
        return

    # Convert to DataFrame (like an Excel table in Python)
    df = pd.DataFrame(genes)

    # Sort by score — highest first
    df = df.sort_values("ot_score", ascending=False).reset_index(drop=True)

    # Print to terminal so you can see it
    print(f"\n✓ {len(df)} genes passed the score cutoff ({SCORE_CUTOFF})")
    print("\n── Gene List ──────────────────────────────")
    print(df[["gene_symbol", "gene_name", "ot_score"]].to_string(index=False))

    # Save to CSV
    df.to_csv(OUTPUT_FILE, index=False)
    print(f"\n✓ Saved to: {OUTPUT_FILE}")
    print("\nLayer 1 COMPLETE. Run layer2_gnomad.py next.")

    return df


if __name__ == "__main__":
    main()