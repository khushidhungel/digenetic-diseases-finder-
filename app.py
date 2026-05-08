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
import codecs

st.set_page_config(page_title="BBS Digenic Explorer", page_icon="🧬", layout="wide", initial_sidebar_state="expanded")

# ── FIX: sidebar visibility ──
st.markdown("""
<style>
[data-testid="stSidebar"] {
    min-width: 260px;
    max-width: 260px;
}
</style>
""", unsafe_allow_html=True)

# ── YOUR ORIGINAL STYLE (UNCHANGED) ──
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Space+Mono:wght@400;700&family=DM+Sans:wght@300;400;500;600&display=swap');
:root{--bg:#0a0e1a;--bg2:#111827;--bg3:#1a2235;--accent1:#00f5c4;--accent2:#6366f1;--accent3:#f59e0b;--accent4:#ef4444;--text:#f8fafc;--text2:#cbd5e1;--border:rgba(255,255,255,0.12);}
html,body,[class*="css"]{font-family:'DM Sans',sans-serif;background-color:var(--bg)!important;color:var(--text)!important;}
.stApp{background-color:var(--bg)!important;}
section[data-testid="stSidebar"]{background:var(--bg2)!important;border-right:1px solid var(--border);}
.metric-box{background:var(--bg3);border:1px solid var(--border);border-radius:10px;padding:1rem;text-align:center;}
.metric-value{font-family:'Space Mono',monospace;font-size:2rem;font-weight:700;color:var(--accent1);line-height:1;}
.metric-label{font-size:0.75rem;color:var(--text2);}
.hero-title{font-family:'Space Mono',monospace;font-size:2.4rem;font-weight:800;color:#fff;}
.section-title{font-family:'Space Mono',monospace;font-size:1rem;color:#fff;}
.card{background:var(--bg3);border:1px solid var(--border);border-radius:12px;padding:1rem;}
#MainMenu{visibility:hidden;}footer{visibility:hidden;}header{visibility:hidden;}
</style>
""", unsafe_allow_html=True)

# ── DATA ──
@st.cache_data
def load_data():
    files = {
        "genes":"layer1_bbs_genes.csv",
        "constrained":"layer2_constrained_genes.csv",
        "interactions":"layer3_interactions.csv",
        "pathways":"layer4_pathway_scores.csv",
        "clinvar":"layer4b_clinvar_scores.csv",
        "scores":"layer5_digenic_scores.csv"
    }

    data = {}
    for k,v in files.items():
        data[k] = pd.read_csv(v) if os.path.exists(v) else pd.DataFrame()
    return data


def load_ai_text():
    path = "layer6_ai_interpretations.txt"
    if os.path.exists(path):
        with codecs.open(path, "r", encoding="utf-8", errors="ignore") as f:
            return f.read()
    return None


data = load_data()
ai_text = load_ai_text()

scores_df = data["scores"]
genes_df = data["constrained"]
inter_df = data["interactions"]

# ── SIDEBAR (UNCHANGED STRUCTURE) ──
with st.sidebar:
    st.markdown("🧬 BBS Explorer")

    page = st.radio(
        "Navigation",
        ["🏠 Overview","🔬 Gene Analysis","🕸️ Network","📊 Digenic Scores","🤖 AI Interpretation","⚙️ Pipeline"]
    )


def score_color(s):
    return "#ef4444" if s>=80 else "#f59e0b" if s>=60 else "#00f5c4" if s>=40 else "#6366f1"


# ─────────────────────────
# OVERVIEW
# ─────────────────────────
if "Overview" in page:

    st.markdown('<div class="hero-title">Digenic Disease Explorer</div>', unsafe_allow_html=True)

    if not scores_df.empty:
        st.metric("Top Score", round(scores_df["digenic_score_normalized"].max(), 1))
        st.metric("Pairs", len(scores_df))


# ─────────────────────────
# GENE ANALYSIS
# ─────────────────────────
elif "Gene" in page:

    st.markdown('<div class="section-title">Gene Analysis</div>', unsafe_allow_html=True)

    if not genes_df.empty:
        st.dataframe(genes_df)


# ─────────────────────────
# NETWORK
# ─────────────────────────
elif "Network" in page:

    st.markdown('<div class="section-title">Network</div>', unsafe_allow_html=True)

    if not inter_df.empty:
        st.dataframe(inter_df)


# ─────────────────────────
# DIGENIC SCORES
# ─────────────────────────
elif "Digenic" in page:

    st.markdown('<div class="section-title">Digenic Scores</div>', unsafe_allow_html=True)

    if not scores_df.empty:
        st.dataframe(scores_df)


# ─────────────────────────
# AI INTERPRETATION (FIXED)
# ─────────────────────────
elif "AI" in page:

    st.markdown('<div class="section-title">AI Interpretation</div>', unsafe_allow_html=True)

    if ai_text:

        sections = [s for s in ai_text.split("="*50) if "RANK" in s]

        for section in sections:

            lines = section.strip().splitlines()

            title = next((l for l in lines if "RANK" in l), "Result")

            content = "\n".join(lines)

            st.markdown(f"### {title}")
            st.text(content)

    else:
        st.info("No AI output found.")


# ─────────────────────────
# PIPELINE
# ─────────────────────────
elif "Pipeline" in page:

    st.markdown('<div class="section-title">Pipeline</div>', unsafe_allow_html=True)

    st.write("All layers loaded (if files exist).")
