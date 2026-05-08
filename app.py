"""
BBS Digenic Variant Explorer — Interactive Streamlit App
=========================================================
Run with: streamlit run app.py
"""

import streamlit as st
import pandas as pd
import networkx as nx
import os
import plotly.graph_objects as go
import plotly.express as px

st.set_page_config(
    page_title="BBS Digenic Explorer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Space+Mono:wght@400;700&family=DM+Sans:wght@300;400;500;600&display=swap');
:root{--bg:#0a0e1a;--bg2:#111827;--bg3:#1a2235;--accent1:#00f5c4;--accent2:#6366f1;--accent3:#f59e0b;--accent4:#ef4444;--text:#e2e8f0;--text2:#94a3b8;--border:rgba(255,255,255,0.08);}
html,body,[class*="css"]{font-family:'DM Sans',sans-serif;background-color:var(--bg)!important;color:var(--text)!important;}
.stApp{background-color:var(--bg)!important;}
section[data-testid="stSidebar"]{background:var(--bg2)!important;border-right:1px solid var(--border);}
.metric-box{background:var(--bg3);border:1px solid var(--border);border-radius:10px;padding:1rem;text-align:center;}
.metric-value{font-family:'Space Mono',monospace;font-size:2rem;font-weight:700;color:var(--accent1);line-height:1;}
.metric-label{font-size:0.75rem;color:var(--text2);margin-top:4px;letter-spacing:0.05em;text-transform:uppercase;}
.hero-title{font-family:'Space Mono',monospace;font-size:2.2rem;font-weight:700;background:linear-gradient(135deg,#00f5c4,#6366f1);-webkit-background-clip:text;-webkit-text-fill-color:transparent;line-height:1.2;}
.section-title{font-family:'Space Mono',monospace;font-size:0.85rem;color:var(--accent1);letter-spacing:0.1em;text-transform:uppercase;margin-bottom:0.5rem;}
.card{background:var(--bg3);border:1px solid var(--border);border-radius:12px;padding:1.2rem 1.5rem;margin-bottom:1rem;}
#MainMenu{visibility:hidden;}footer{visibility:hidden;}header{visibility:hidden;}
</style>
""", unsafe_allow_html=True)

# ══════════════════════════════════════════════════
# LOAD DATA (FIXED: NO outputs/)
# ══════════════════════════════════════════════════
@st.cache_data
def load_data():
    data = {}
    files = {
        "genes": "layer1_bbs_genes.csv",
        "constrained": "layer2_constrained_genes.csv",
        "interactions": "layer3_interactions.csv",
        "pathways": "layer4_pathway_scores.csv",
        "clinvar": "layer4b_clinvar_scores.csv",
        "scores": "layer5_digenic_scores.csv"
    }

    for key, path in files.items():
        if os.path.exists(path):
            data[key] = pd.read_csv(path, encoding="utf-8")
        else:
            data[key] = pd.DataFrame()

    return data

def load_ai_text():
    path = "layer6_ai_interpretations.txt"
    if os.path.exists(path):
        with open(path, "r", encoding="utf-8") as f:
            return f.read()
    return None

data = load_data()
ai_text = load_ai_text()

scores_df = data["scores"]
genes_df = data["constrained"]
inter_df = data["interactions"]

# ══════════════════════════════════════════════════
# SIDEBAR
# ══════════════════════════════════════════════════
with st.sidebar:
    st.markdown("🧬 **BBS Explorer**")
    page = st.radio(
        "",
        ["🏠  Overview","🔬  Gene Analysis","🕸️  Network","📊  Digenic Scores","🤖  AI Interpretation","⚙️  Pipeline"],
        label_visibility="collapsed"
    )

# ══════════════════════════════════════════════════
# OVERVIEW
# ══════════════════════════════════════════════════
if page == "🏠 Overview":

    st.title("BBS Digenic Explorer")

    col1, col2, col3 = st.columns(3)

    col1.metric("Genes", len(genes_df))
    col2.metric("Interactions", len(inter_df))
    col3.metric("Digenic Pairs", len(scores_df))

    st.subheader("Top Digenic Pairs")
    st.dataframe(scores_df.head(10), use_container_width=True)

# ══════════════════════════════════════════════════
# GENE ANALYSIS
# ══════════════════════════════════════════════════
elif page == "🔬 Gene Analysis":

    gene = st.selectbox("Select Gene", sorted(genes_df.iloc[:,0].dropna().unique()))

    result = scores_df[
        (scores_df["gene_a"] == gene) |
        (scores_df["gene_b"] == gene)
    ]

    st.dataframe(result, use_container_width=True)

# ══════════════════════════════════════════════════
# NETWORK
# ══════════════════════════════════════════════════
elif page == "🕸️ Network":

    G = nx.Graph()

    for _, row in inter_df.iterrows():
        if len(row) >= 2:
            G.add_edge(row.iloc[0], row.iloc[1])

    fig, ax = plt.subplots(figsize=(10, 8))
    pos = nx.spring_layout(G, seed=42)
    nx.draw(G, pos, with_labels=True, node_size=2000, ax=ax)

    st.pyplot(fig)

# ══════════════════════════════════════════════════
# SCORES
# ══════════════════════════════════════════════════
elif page == "📊 Digenic Scores":

    st.dataframe(
        scores_df.sort_values("digenic_score_normalized", ascending=False),
        use_container_width=True
    )

# ══════════════════════════════════════════════════
# AI INTERPRETATION
# ══════════════════════════════════════════════════
elif page == "🤖 AI Interpretation":

    st.subheader("AI Interpretation")

    if ai_text:
        st.text(ai_text)
    else:
        st.warning("No AI file found")

# ══════════════════════════════════════════════════
# PIPELINE
# ══════════════════════════════════════════════════
elif page == "⚙️ Pipeline":

    st.code("""
Layer 1 → Gene collection
Layer 2 → Constraint filtering
Layer 3 → STRING interactions
Layer 4 → Pathway scoring
Layer 4B → ClinVar scoring
Layer 5 → Digenic ranking
Layer 6 → AI interpretation
""")
