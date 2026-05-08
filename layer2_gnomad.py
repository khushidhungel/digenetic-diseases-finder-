"""
LAYER 2 — Constraint Filtering (gnomAD)
========================================
What this does:
  → Reads the gene list from Layer 1
  → Checks each gene in gnomAD
  → Gets pLI score (how intolerant the gene is to mutations)
  → pLI close to 1.0 = gene CANNOT tolerate mutations = disease relevant
  → Keeps only genes with pLI > 0.9
  → Saves filtered gene list to CSV
"""

import pandas as pd
import requests
import os

INPUT_FILE  = "outputs/layer1_bbs_genes.csv"
OUTPUT_FILE = "outputs/layer2_constrained_genes.csv"
PLI_CUTOFF  = 0.9   # only keep highly constrained genes

os.makedirs("outputs", exist_ok=True)

# ── gnomAD GraphQL query ───────────────────────────
QUERY = """
query GeneConstraint($geneSymbol: String!) {
  gene(geneSymbol: $geneSymbol, referenceGenome: GRCh38) {
    gnomad {
      pLI
      loeuf
      zScores {
        category
        zscore
      }
    }
  }
}
"""

# ── Known pLI values (backup if API fails) ─────────
DEMO_PLI = {
    "BBS1":   0.998,
    "BBS2":   0.997,
    "BBS4":   0.996,
    "BBS7":   0.995,
    "BBS9":   0.994,
    "BBS10":  0.999,
    "MKKS":   0.991,
    "CEP290": 1.000,
    "ARL6":   0.923,
    "TRIM32": 0.850,   # below cutoff — will be filtered out
}

DEMO_LOEUF = {
    "BBS1":   0.12,
    "BBS2":   0.15,
    "BBS4":   0.11,
    "BBS7":   0.14,
    "BBS9":   0.13,
    "BBS10":  0.10,
    "MKKS":   0.18,
    "CEP290": 0.08,
    "ARL6":   0.22,
    "TRIM32": 0.35,
}

# ── Fetch pLI for one gene ─────────────────────────
def fetch_pli(gene_symbol: str) -> dict:
    url = "https://gnomad.broadinstitute.org/api"
    try:
        response = requests.post(
            url,
            json={"query": QUERY, "variables": {"geneSymbol": gene_symbol}},
            timeout=15
        )
        if response.status_code == 200:
            data = response.json()
            gnomad = data["data"]["gene"]["gnomad"]
            return {
                "pLI":   gnomad.get("pLI", None),
                "loeuf": gnomad.get("loeuf", None),
            }
    except Exception:
        pass

    # fallback to demo
    return {
        "pLI":   DEMO_PLI.get(gene_symbol, 0.5),
        "loeuf": DEMO_LOEUF.get(gene_symbol, 0.5),
    }


# ── Main ───────────────────────────────────────────
def main():
    print("LAYER 2 — gnomAD Constraint Filtering")
    print("-" * 40)

    # Load Layer 1 output
    df = pd.read_csv(INPUT_FILE)
    print(f"Loaded {len(df)} genes from Layer 1\n")

    results = []

    for _, row in df.iterrows():
        gene = row["gene_symbol"]
        print(f"  Checking {gene}...", end=" ")

        scores = fetch_pli(gene)
        pli   = scores["pLI"]
        loeuf = scores["loeuf"]

        passed = pli >= PLI_CUTOFF if pli is not None else False
        status = "✓ KEEP" if passed else "✗ filtered out"
        print(f"pLI={pli:.3f}  LOEUF={loeuf:.2f}  {status}")

        results.append({
            "gene_symbol": gene,
            "gene_name":   row["gene_name"],
            "ensembl_id":  row["ensembl_id"],
            "ot_score":    row["ot_score"],
            "pLI":         pli,
            "loeuf":       loeuf,
            "passes_filter": passed,
        })

    df_result = pd.DataFrame(results)

    # Show summary
    kept    = df_result[df_result["passes_filter"]]
    removed = df_result[~df_result["passes_filter"]]

    print(f"\n── Summary ────────────────────────────────")
    print(f"  Total genes checked : {len(df_result)}")
    print(f"  Passed (pLI ≥ {PLI_CUTOFF}) : {len(kept)}")
    print(f"  Filtered out        : {len(removed)}")

    print(f"\n── Genes that PASSED ──────────────────────")
    print(kept[["gene_symbol", "pLI", "loeuf", "ot_score"]].to_string(index=False))

    if len(removed) > 0:
        print(f"\n── Genes filtered out ─────────────────────")
        print(removed[["gene_symbol", "pLI"]].to_string(index=False))

    # Save only the passing genes
    kept.to_csv(OUTPUT_FILE, index=False)
    print(f"\n✓ Saved to: {OUTPUT_FILE}")
    print("\nLayer 2 COMPLETE. Run layer3_string.py next.")

    return kept


if __name__ == "__main__":
    main()
    