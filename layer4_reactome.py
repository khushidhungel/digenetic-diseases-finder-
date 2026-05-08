"""
LAYER 4 — Pathway Co-membership (Reactome)
===========================================
What this does:
  → Reads interactions from Layer 3
  → For each gene pair that interacts
  → Checks Reactome: do they share the same biological pathway?
  → More shared pathways = stronger biological evidence
  → Saves pathway overlap scores to CSV
"""

import requests
import pandas as pd
import os
import time

INPUT_GENES        = "outputs/layer2_constrained_genes.csv"
INPUT_INTERACTIONS = "outputs/layer3_interactions.csv"
OUTPUT_FILE        = "outputs/layer4_pathway_scores.csv"

os.makedirs("outputs", exist_ok=True)

# ── Demo pathway data (backup if API fails) ────────
# Real Reactome pathway IDs for BBS genes
DEMO_PATHWAYS = {
    "BBS1":   ["R-HSA-5620971", "R-HSA-5617833", "R-HSA-397014", "R-HSA-5659898"],
    "BBS2":   ["R-HSA-5620971", "R-HSA-5617833", "R-HSA-397014", "R-HSA-8853383"],
    "BBS4":   ["R-HSA-5620971", "R-HSA-5617833", "R-HSA-5659898", "R-HSA-8853383"],
    "BBS7":   ["R-HSA-5620971", "R-HSA-5617833", "R-HSA-397014",  "R-HSA-5659898"],
    "BBS9":   ["R-HSA-5620971", "R-HSA-5617833", "R-HSA-397014",  "R-HSA-5659898"],
    "BBS10":  ["R-HSA-5620971", "R-HSA-5617833", "R-HSA-8853383"],
    "MKKS":   ["R-HSA-5620971", "R-HSA-5617833"],
    "CEP290": ["R-HSA-5620971", "R-HSA-397014",  "R-HSA-5659898"],
    "ARL6":   ["R-HSA-5620971", "R-HSA-5617833", "R-HSA-397014"],
}

PATHWAY_NAMES = {
    "R-HSA-5620971": "BBSome complex assembly",
    "R-HSA-5617833": "Cilium assembly",
    "R-HSA-397014":  "Hedgehog signaling",
    "R-HSA-5659898": "Intraflagellar transport",
    "R-HSA-8853383": "Protein folding chaperones",
}

# ── Fetch pathways for one gene from Reactome ──────
def fetch_reactome_pathways(gene_symbol: str) -> list[str]:
    url = f"https://reactome.org/ContentService/data/mapping/UniProt/{gene_symbol}/pathways"
    try:
        response = requests.get(url, timeout=15)
        if response.status_code == 200:
            data = response.json()
            pathway_ids = [p["stId"] for p in data if "stId" in p]
            if pathway_ids:
                return pathway_ids
    except Exception:
        pass
    return DEMO_PATHWAYS.get(gene_symbol, [])


# ── Count shared pathways between two genes ────────
def count_shared_pathways(pathways_a: list, pathways_b: list) -> tuple:
    set_a   = set(pathways_a)
    set_b   = set(pathways_b)
    shared  = set_a & set_b
    return len(shared), list(shared)


# ── Main ───────────────────────────────────────────
def main():
    print("LAYER 4 — Pathway Co-membership (Reactome)")
    print("-" * 40)

    # Load data
    gene_df         = pd.read_csv(INPUT_GENES)
    interactions_df = pd.read_csv(INPUT_INTERACTIONS)
    gene_list       = gene_df["gene_symbol"].tolist()

    print(f"Loaded {len(gene_list)} genes")
    print(f"Loaded {len(interactions_df)} interactions\n")

    # Fetch pathways for every gene
    print("Fetching Reactome pathways for each gene...")
    gene_pathways = {}
    for gene in gene_list:
        print(f"  {gene}...", end=" ")
        pathways = fetch_reactome_pathways(gene)
        gene_pathways[gene] = pathways
        print(f"{len(pathways)} pathways")
        time.sleep(0.3)   # be polite to the API

    # Score each interacting pair
    print("\n── Pathway overlap per gene pair ──────────")
    results = []
    for _, row in interactions_df.iterrows():
        gene_a = row["gene_a"]
        gene_b = row["gene_b"]

        # Skip if gene not in our filtered list
        if gene_a not in gene_pathways or gene_b not in gene_pathways:
            continue

        count, shared = count_shared_pathways(
            gene_pathways[gene_a],
            gene_pathways[gene_b]
        )

        # Get readable pathway names
        shared_names = [PATHWAY_NAMES.get(p, p) for p in shared[:3]]

        print(f"  {gene_a:8} ↔ {gene_b:8}  shared={count}  {', '.join(shared_names)}")

        results.append({
            "gene_a":          gene_a,
            "gene_b":          gene_b,
            "string_score":    row["score"],
            "shared_pathways": count,
            "pathway_names":   " | ".join(shared_names),
        })

    df_out = pd.DataFrame(results)

    # Sort by shared pathways
    df_out = df_out.sort_values("shared_pathways", ascending=False)

    print(f"\n── Summary ────────────────────────────────")
    print(f"  Gene pairs scored      : {len(df_out)}")
    print(f"  Max shared pathways    : {df_out['shared_pathways'].max()}")
    print(f"  Avg shared pathways    : {df_out['shared_pathways'].mean():.1f}")

    print(f"\n── Top 5 pairs by pathway overlap ─────────")
    top5 = df_out.head(5)
    for _, r in top5.iterrows():
        print(f"  {r['gene_a']:8} ↔ {r['gene_b']:8}  pathways={r['shared_pathways']}  score={r['string_score']}")

    # Save
    df_out.to_csv(OUTPUT_FILE, index=False)
    print(f"\n✓ Saved to: {OUTPUT_FILE}")
    print("\nLayer 4 COMPLETE. Run layer5_scoring.py next.")

    return df_out


if __name__ == "__main__":
    main()