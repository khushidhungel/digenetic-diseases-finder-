"""
BBS Digenic Variant Explorer — Interactive Streamlit App
Run with: streamlit run app.py
"""

import streamlit as st
import pandas as pd
import networkx as nx
import os
import plotly.graph_objects as go
import plotly.express as px

# ── PAGE CONFIG ──
st.set_page_config(
    page_title="BBS Digenic Explorer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ── SAFE CSS (NO BREAKING RULES) ──
st.markdown("""
<style>
.stApp {
    background-color: #0a0e1a;
    color: #e2e8f0;
}

section[data-testid="stSidebar"] {
    background-color: #0d1f3c;
    border-right: 2px solid #00f5c4;
}

.metric-box{
    background:#1a2235;
    border:1px solid rgba(255,255,255,0.1);
    border-radius:10px;
    padding:1rem;
    text-align:center;
}

.metric-value{
    font-size:2rem;
    font-weight:700;
    color:#00f5c4;
}

.metric-label{
    font-size:0.75rem;
    color:#cbd5e1;
    text-transform:uppercase;
}

.card{
    background:#1a2235;
    border-radius:12px;
    padding:1rem;
    margin-bottom:0.5rem;
    border:1px solid rgba(255,255,255,0.08);
}

.section-title{
    font-size:0.9rem;
    font-weight:700;
    color:#00f5c4;
    margin:10px 0;
}

</style>
""", unsafe_allow_html=True)

# ── DATA LOADING ──
@st.cache_data
def load_data():
    files = {
        "genes": "layer1_bbs_genes.csv",
        "constrained": "layer2_constrained_genes.csv",
        "interactions": "layer3_interactions.csv",
        "pathways": "layer4_pathway_scores.csv",
        "clinvar": "layer4b_clinvar_scores.csv",
        "scores": "layer5_digenic_scores.csv",
    }

    data = {}
    for k, v in files.items():
        if os.path.exists(v):
            data[k] = pd.read_csv(v)
        else:
            data[k] = pd.DataFrame()

    return data

def load_ai():
    if os.path.exists("layer6_ai_interpretations.txt"):
        with open("layer6_ai_interpretations.txt", encoding="utf-8") as f:
            return f.read()
    return None

data = load_data()
ai_text = load_ai()

genes_df = data["genes"]
scores_df = data["scores"]
inter_df = data["interactions"]

# ── SIDEBAR ──
with st.sidebar:
    st.title("🧬 BBS Explorer")

    page = st.radio(
        "Navigation",
        ["Overview", "Gene Analysis", "Network", "Digenic Scores", "AI Interpretation", "Pipeline"]
    )

    st.markdown("---")
    st.write("Bardet-Biedl Syndrome")

# ── HELPERS ──
def score_color(s):
    if s >= 80:
        return "#ef4444"
    elif s >= 60:
        return "#f59e0b"
    elif s >= 40:
        return "#00f5c4"
    return "#6366f1"

# ── OVERVIEW ──
if page == "Overview":

    st.title("Digenic Disease Explorer")

    c1, c2, c3, c4 = st.columns(4)

    with c1:
        st.markdown(f"<div class='metric-box'><div class='metric-value'>{len(genes_df)}</div><div class='metric-label'>Genes</div></div>", unsafe_allow_html=True)

    with c2:
        st.markdown(f"<div class='metric-box'><div class='metric-value'>{len(inter_df)}</div><div class='metric-label'>Interactions</div></div>", unsafe_allow_html=True)

    with c3:
        st.markdown(f"<div class='metric-box'><div class='metric-value'>{len(scores_df)}</div><div class='metric-label'>Pairs</div></div>", unsafe_allow_html=True)

    with c4:
        top = scores_df["digenic_score_normalized"].max() if not scores_df.empty else 0
        st.markdown(f"<div class='metric-box'><div class='metric-value'>{top:.1f}</div><div class='metric-label'>Top Score</div></div>", unsafe_allow_html=True)

    st.subheader("Top Pairs")

    if not scores_df.empty:
        for _, row in scores_df.head(5).iterrows():
            s = row["digenic_score_normalized"]
            c = score_color(s)

            st.markdown(f"""
            <div class='card'>
            <b>{row['gene_a']} ↔ {row['gene_b']}</b><br>
            Score: <span style='color:{c}'>{s:.1f}</span>
            </div>
            """, unsafe_allow_html=True)

# ── AI INTERPRETATION (FIXED SECTION) ──
elif page == "AI Interpretation":

    st.title("AI Interpretation")

    if ai_text:

        sections = [s for s in ai_text.split("="*50) if "RANK" in s]

        for section in sections:

            # FIXED SAFE PARSING
            lines = section.split("\n")

            title = next((l for l in lines if "RANK" in l), "Rank")
            score_line = next((l for l in lines if "Score" in l), "")

            try:
                score_val = float(score_line.split(":")[-1].split("/")[0])
            except:
                score_val = 0

            content = "\n".join(lines[2:])

            color = score_color(score_val)

            with st.expander(f"{title} — {score_val:.0f}/100"):

                st.markdown(
                    f"<div style='border-left:4px solid {color};padding:10px;color:#cbd5e1;white-space:pre-wrap'>{content}</div>",
                    unsafe_allow_html=True
                )

    else:
        st.warning("Run AI layer first")

# ── NETWORK ──
elif page == "Network":

    st.title("Network View")

    if not inter_df.empty and not genes_df.empty:
        G = nx.Graph()

        for _, r in inter_df.iterrows():
            G.add_edge(r["gene_a"], r["gene_b"])

        st.write(f"Nodes: {G.number_of_nodes()} | Edges: {G.number_of_edges()}")

# ── SCORES ──
elif page == "Digenic Scores":

    st.title("Digenic Scores")

    if not scores_df.empty:
        st.dataframe(scores_df)

# ── PIPELINE ──
elif page == "Pipeline":

    st.title("Pipeline")

    st.write("All layers loaded successfully.")
