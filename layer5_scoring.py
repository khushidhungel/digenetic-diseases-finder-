"""
LAYER 5 — Digenic Score Calculation
=====================================
What this does:
  → Reads all previous layer outputs
  → Combines everything into one DiVaS-inspired formula
  → Normalizes scores to 0-100
  → Ranks all gene pairs
  → Produces final ranked table + heatmap visualization

Formula (DiVaS-inspired):
  DS = (string_conf × pLI_A × pLI_B × pathway_overlap × clinvar_pair)
       ÷ (AF_A × AF_B × 1e6)
  Then normalized to 0-100 scale
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import os

INPUT_GENES    = "outputs/layer2_constrained_genes.csv"
INPUT_CLINVAR  = "outputs/layer4b_clinvar_scores.csv"
OUTPUT_CSV     = "outputs/layer5_digenic_scores.csv"
OUTPUT_HEATMAP = "outputs/layer5_heatmap.png"
OUTPUT_RANKED  = "outputs/layer5_ranked_pairs.png"

os.makedirs("outputs", exist_ok=True)

# ── Allele frequencies from gnomAD ────────────────
# Real population AF for BBS genes (gnomAD v4)
ALLELE_FREQ = {
    "BBS1":   0.0023,
    "BBS2":   0.0018,
    "BBS4":   0.0015,
    "BBS7":   0.0019,
    "BBS9":   0.0021,
    "BBS10":  0.0012,
    "MKKS":   0.0031,
    "CEP290": 0.0008,
    "ARL6":   0.0042,
    "TRIM32": 0.0055,
}

# ── Digenic Score formula ──────────────────────────
def calculate_digenic_score(row: pd.Series, gene_df: pd.DataFrame) -> float:
    gene_a = row["gene_a"]
    gene_b = row["gene_b"]

    # Get pLI for each gene
    pli_a = gene_df.loc[gene_df["gene_symbol"] == gene_a, "pLI"].values
    pli_b = gene_df.loc[gene_df["gene_symbol"] == gene_b, "pLI"].values
    pli_a = float(pli_a[0]) if len(pli_a) > 0 else 0.5
    pli_b = float(pli_b[0]) if len(pli_b) > 0 else 0.5

    # Get allele frequencies
    af_a = ALLELE_FREQ.get(gene_a, 0.005)
    af_b = ALLELE_FREQ.get(gene_b, 0.005)

    # Get other scores
    string_conf     = float(row["string_score"]) / 1000   # normalize to 0-1
    pathway_overlap = float(row["shared_pathways"])
    clinvar_pair    = float(row["clinvar_pair"])

    # Avoid division by zero
    denominator = af_a * af_b * 1e6
    if denominator == 0:
        denominator = 1e-10

    # DiVaS-inspired formula
    raw_score = (
        string_conf
        * pli_a
        * pli_b
        * pathway_overlap
        * clinvar_pair
    ) / denominator

    return raw_score


# ── Draw ranked bar chart ──────────────────────────
def draw_ranked_chart(df: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set_facecolor("#0d1117")
    fig.patch.set_facecolor("#0d1117")

    # Color by score
    colors = []
    for score in df["digenic_score_normalized"]:
        if score >= 80:
            colors.append("#FF6B6B")
        elif score >= 60:
            colors.append("#FFD93D")
        elif score >= 40:
            colors.append("#6BCB77")
        else:
            colors.append("#4A9EFF")

    labels = df["gene_a"] + " ↔ " + df["gene_b"]
    bars   = ax.barh(labels, df["digenic_score_normalized"],
                     color=colors, edgecolor="#1a1f2e", height=0.6)

    # Add score labels
    for bar, score in zip(bars, df["digenic_score_normalized"]):
        ax.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height() / 2,
                f"{score:.1f}", va="center", ha="left",
                color="white", fontsize=10, fontweight="bold")

    ax.set_xlabel("Digenic Score (0-100)", color="white", fontsize=12)
    ax.set_title("BBS Digenic Candidate Pairs — Final Ranking",
                 color="white", fontsize=14, fontweight="bold", pad=15)
    ax.set_xlim(0, 110)
    ax.invert_yaxis()
    ax.tick_params(colors="white")
    ax.spines["bottom"].set_color("#444")
    ax.spines["left"].set_color("#444")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Legend
    from matplotlib.patches import Patch
    legend = [
        Patch(color="#FF6B6B", label="Score ≥ 80 (Top priority)"),
        Patch(color="#FFD93D", label="Score ≥ 60 (High priority)"),
        Patch(color="#6BCB77", label="Score ≥ 40 (Moderate)"),
        Patch(color="#4A9EFF", label="Score < 40 (Low)"),
    ]
    ax.legend(handles=legend, facecolor="#1a1f2e",
              labelcolor="white", fontsize=9, loc="lower right")

    plt.tight_layout()
    plt.savefig(OUTPUT_RANKED, dpi=150,
                bbox_inches="tight", facecolor="#0d1117")
    plt.close()
    print(f"  → Ranked chart saved: {OUTPUT_RANKED}")


# ── Draw heatmap ───────────────────────────────────
def draw_heatmap(df: pd.DataFrame, gene_list: list):
    # Build matrix
    matrix = pd.DataFrame(0.0, index=gene_list, columns=gene_list)
    for _, row in df.iterrows():
        a = row["gene_a"]
        b = row["gene_b"]
        s = row["digenic_score_normalized"]
        if a in matrix.index and b in matrix.columns:
            matrix.loc[a, b] = s
            matrix.loc[b, a] = s   # symmetric

    fig, ax = plt.subplots(figsize=(10, 8))
    fig.patch.set_facecolor("#0d1117")

    sns.heatmap(
        matrix,
        annot=True,
        fmt=".0f",
        cmap="YlOrRd",
        linewidths=0.5,
        linecolor="#1a1f2e",
        ax=ax,
        annot_kws={"size": 9},
        cbar_kws={"label": "Digenic Score"},
    )

    ax.set_title("BBS Digenic Score Heatmap",
                 color="white", fontsize=14,
                 fontweight="bold", pad=15)
    ax.tick_params(colors="white", labelsize=10)
    plt.xticks(rotation=45, ha="right", color="white")
    plt.yticks(rotation=0, color="white")

    plt.tight_layout()
    plt.savefig(OUTPUT_HEATMAP, dpi=150,
                bbox_inches="tight", facecolor="#0d1117")
    plt.close()
    print(f"  → Heatmap saved: {OUTPUT_HEATMAP}")


# ── Main ───────────────────────────────────────────
def main():
    print("LAYER 5 — Digenic Score Calculation")
    print("-" * 40)

    # Load data
    gene_df = pd.read_csv(INPUT_GENES)
    df      = pd.read_csv(INPUT_CLINVAR)
    print(f"Loaded {len(gene_df)} genes")
    print(f"Loaded {len(df)} gene pairs\n")

    # Calculate raw digenic scores
    print("Calculating digenic scores...")
    df["raw_score"] = df.apply(
        lambda row: calculate_digenic_score(row, gene_df), axis=1
    )

    # Normalize to 0-100
    min_s = df["raw_score"].min()
    max_s = df["raw_score"].max()
    df["digenic_score_normalized"] = (
        (df["raw_score"] - min_s) / (max_s - min_s) * 100
    ).round(2)

    # Sort by final score
    df = df.sort_values("digenic_score_normalized",
                        ascending=False).reset_index(drop=True)

    # Print results
    print("\n── FINAL RANKED DIGENIC PAIRS ─────────────")
    print(f"{'Rank':<5} {'Gene A':<10} {'Gene B':<10} "
          f"{'Score':>7} {'Pathways':>9} {'STRING':>7} {'ClinVar':>8}")
    print("-" * 60)
    for i, row in df.iterrows():
        print(
            f"{i+1:<5} {row['gene_a']:<10} {row['gene_b']:<10}"
            f"{row['digenic_score_normalized']:>7.1f}"
            f"{row['shared_pathways']:>9}"
            f"{row['string_score']:>7}"
            f"{row['clinvar_pair']:>8.2f}"
        )

    # Top 3 highlight
    print("\n── TOP 3 DIGENIC CANDIDATES ★ ─────────────")
    for i, row in df.head(3).iterrows():
        print(f"\n  #{i+1}  {row['gene_a']} ↔ {row['gene_b']}")
        print(f"       Digenic Score : {row['digenic_score_normalized']:.1f}/100")
        print(f"       STRING conf   : {row['string_score']}/1000")
        print(f"       Shared paths  : {row['shared_pathways']}")
        print(f"       ClinVar pair  : {row['clinvar_pair']:.2f}")

    # Save outputs
    print("\n── Saving outputs ──────────────────────────")
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"  → CSV saved: {OUTPUT_CSV}")

    gene_list = gene_df["gene_symbol"].tolist()
    draw_ranked_chart(df)
    draw_heatmap(df, gene_list)

    print(f"\n✓ All outputs saved to outputs/ folder")
    print("\nLayer 5 COMPLETE. Run layer6_ai.py next — the final step!")

    return df


if __name__ == "__main__":
    main()