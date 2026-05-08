"""
LAYER 3 — PPI Network Construction (STRING DB)
===============================================
What this does:
  → Reads constrained genes from Layer 2
  → Queries STRING DB for protein-protein interactions
  → Filters by confidence score > 700 (high confidence only)
  → Builds a NetworkX graph
  → Saves network image + interaction table to outputs/
"""

import requests
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os

INPUT_FILE        = "outputs/layer2_constrained_genes.csv"
OUTPUT_CSV        = "outputs/layer3_interactions.csv"
OUTPUT_IMAGE      = "outputs/layer3_network.png"
STRING_SCORE_MIN  = 700   # out of 1000 — high confidence only
SPECIES_ID        = 9606  # 9606 = Human

os.makedirs("outputs", exist_ok=True)

# ── Demo interactions (backup if API fails) ────────
DEMO_INTERACTIONS = [
    {"gene_a": "BBS1",   "gene_b": "BBS2",   "score": 980},
    {"gene_a": "BBS1",   "gene_b": "BBS4",   "score": 965},
    {"gene_a": "BBS1",   "gene_b": "BBS7",   "score": 972},
    {"gene_a": "BBS1",   "gene_b": "BBS9",   "score": 958},
    {"gene_a": "BBS2",   "gene_b": "BBS7",   "score": 975},
    {"gene_a": "BBS2",   "gene_b": "BBS9",   "score": 961},
    {"gene_a": "BBS4",   "gene_b": "BBS8",   "score": 944},
    {"gene_a": "BBS7",   "gene_b": "BBS9",   "score": 968},
    {"gene_a": "BBS7",   "gene_b": "BBS10",  "score": 955},
    {"gene_a": "BBS9",   "gene_b": "BBS10",  "score": 949},
    {"gene_a": "BBS10",  "gene_b": "MKKS",   "score": 938},
    {"gene_a": "MKKS",   "gene_b": "BBS2",   "score": 921},
    {"gene_a": "CEP290", "gene_b": "BBS9",   "score": 912},
    {"gene_a": "CEP290", "gene_b": "ARL6",   "score": 903},
    {"gene_a": "ARL6",   "gene_b": "BBS1",   "score": 895},
    {"gene_a": "ARL6",   "gene_b": "BBS7",   "score": 878},
]

# ── Fetch interactions from STRING ─────────────────
def fetch_string_interactions(gene_list: list[str]) -> list[dict]:
    print("Querying STRING DB for interactions...")
    genes_str = "%0d".join(gene_list)
    url = (
        f"https://string-db.org/api/json/network"
        f"?identifiers={genes_str}"
        f"&species={SPECIES_ID}"
        f"&required_score={STRING_SCORE_MIN}"
        f"&caller_identity=bbs_digenic_hackathon"
    )
    try:
        response = requests.get(url, timeout=20)
        if response.status_code == 200:
            data = response.json()
            if data:
                interactions = []
                for item in data:
                    interactions.append({
                        "gene_a": item["preferredName_A"],
                        "gene_b": item["preferredName_B"],
                        "score":  int(item["score"] * 1000),
                    })
                print(f"  → {len(interactions)} interactions found from STRING API")
                return interactions
    except Exception as e:
        print(f"  STRING API error: {e}")

    print("  Using demo interaction data...")
    return DEMO_INTERACTIONS


# ── Build NetworkX graph ───────────────────────────
def build_graph(interactions: list[dict], gene_df: pd.DataFrame) -> nx.Graph:
    G = nx.Graph()

    # Add nodes with pLI as attribute
    pli_map = dict(zip(gene_df["gene_symbol"], gene_df["pLI"]))
    for gene in gene_df["gene_symbol"]:
        G.add_node(gene, pLI=pli_map.get(gene, 0.5))

    # Add edges with score as weight
    for inter in interactions:
        a = inter["gene_a"]
        b = inter["gene_b"]
        s = inter["score"]
        if a in G.nodes and b in G.nodes:
            G.add_edge(a, b, weight=s)

    return G


# ── Draw and save network image ────────────────────
def draw_network(G: nx.Graph, gene_df: pd.DataFrame):
    print("Drawing network visualization...")

    fig, ax = plt.subplots(figsize=(12, 9))
    ax.set_facecolor("#0d1117")
    fig.patch.set_facecolor("#0d1117")

    # Layout
    pos = nx.spring_layout(G, seed=42, k=2.5)

    # Node color by pLI
    pli_map   = nx.get_node_attributes(G, "pLI")
    node_colors = []
    for node in G.nodes:
        pli = pli_map.get(node, 0.5)
        if pli >= 0.99:
            node_colors.append("#FF6B6B")   # red = extremely constrained
        elif pli >= 0.95:
            node_colors.append("#FFD93D")   # yellow = highly constrained
        else:
            node_colors.append("#6BCB77")   # green = moderately constrained

    # Edge width by score
    edge_weights = [G[u][v]["weight"] / 300 for u, v in G.edges()]

    # Draw
    nx.draw_networkx_edges(G, pos, width=edge_weights,
                           edge_color="#4A9EFF", alpha=0.6, ax=ax)
    nx.draw_networkx_nodes(G, pos, node_color=node_colors,
                           node_size=1800, alpha=0.95, ax=ax)
    nx.draw_networkx_labels(G, pos, font_size=11,
                            font_color="white", font_weight="bold", ax=ax)

    # Legend
    legend = [
        mpatches.Patch(color="#FF6B6B", label="pLI ≥ 0.99 (extremely constrained)"),
        mpatches.Patch(color="#FFD93D", label="pLI ≥ 0.95 (highly constrained)"),
        mpatches.Patch(color="#6BCB77", label="pLI < 0.95 (moderately constrained)"),
    ]
    ax.legend(handles=legend, loc="upper left",
              facecolor="#1a1f2e", labelcolor="white", fontsize=9)

    ax.set_title("BBS Gene Interaction Network (STRING DB)",
                 color="white", fontsize=14, fontweight="bold", pad=15)
    ax.axis("off")

    plt.tight_layout()
    plt.savefig(OUTPUT_IMAGE, dpi=150,
                bbox_inches="tight", facecolor="#0d1117")
    plt.close()
    print(f"  → Network image saved: {OUTPUT_IMAGE}")


# ── Main ───────────────────────────────────────────
def main():
    print("LAYER 3 — PPI Network Construction")
    print("-" * 40)

    # Load Layer 2 output
    gene_df   = pd.read_csv(INPUT_FILE)
    gene_list = gene_df["gene_symbol"].tolist()
    print(f"Loaded {len(gene_list)} genes: {', '.join(gene_list)}\n")

    # Fetch interactions
    interactions = fetch_string_interactions(gene_list)

    # Filter by score
    interactions = [i for i in interactions if i["score"] >= STRING_SCORE_MIN]
    print(f"\n── Interactions (score ≥ {STRING_SCORE_MIN}) ──────────────")
    for i in interactions:
        print(f"  {i['gene_a']:8} ↔ {i['gene_b']:8}  score={i['score']}")

    # Build graph
    G = build_graph(interactions, gene_df)
    print(f"\n── Network stats ──────────────────────────")
    print(f"  Nodes (genes)       : {G.number_of_nodes()}")
    print(f"  Edges (interactions): {G.number_of_edges()}")

    # Top connected genes
    degree = sorted(G.degree(), key=lambda x: x[1], reverse=True)
    print(f"\n── Most connected genes (hub candidates) ──")
    for gene, deg in degree[:5]:
        print(f"  {gene:10} — {deg} interactions")

    # Draw network
    draw_network(G, gene_df)

    # Save interaction table
    df_out = pd.DataFrame(interactions)
    df_out.to_csv(OUTPUT_CSV, index=False)
    print(f"\n✓ Saved interactions: {OUTPUT_CSV}")
    print("✓ Saved network image: outputs/layer3_network.png")
    print("\nLayer 3 COMPLETE. Run layer4_reactome.py next.")

    return G, interactions


if __name__ == "__main__":
    main()