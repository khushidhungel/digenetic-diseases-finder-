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
import codecs  # FIX for encoding issues

st.set_page_config(page_title="BBS Digenic Explorer", page_icon="🧬", layout="wide", initial_sidebar_state="expanded")

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

# ── LOAD DATA ──
@st.cache_data
def load_data():
    data = {}
    files = {
        "genes":"layer1_bbs_genes.csv",
        "constrained":"layer2_constrained_genes.csv",
        "interactions":"layer3_interactions.csv",
        "pathways":"layer4_pathway_scores.csv",
        "clinvar":"layer4b_clinvar_scores.csv",
        "scores":"layer5_digenic_scores.csv"
    }

    for key, path in files.items():
        data[key] = pd.read_csv(path) if os.path.exists(path) else pd.DataFrame()

    return data

def load_ai_text():
    path = "layer6_ai_interpretations.txt"
    if os.path.exists(path):
        return codecs.open(path, encoding="utf-8", errors="ignore").read()
    return None

data = load_data()
ai_text = load_ai_text()

scores_df = data["scores"]
genes_df = data["constrained"]
inter_df = data["interactions"]

# ── SIDEBAR ──
with st.sidebar:
    st.markdown("🧬 BBS Explorer")
    page = st.radio("Navigation", ["Overview","Gene Analysis","Network","Digenic Scores","AI Interpretation","Pipeline"])

# ── COLOR FUNCTION ──
def score_color(s):
    if s >= 80:
        return "#ef4444"
    elif s >= 60:
        return "#f59e0b"
    elif s >= 40:
        return "#00f5c4"
    else:
        return "#6366f1"

# ── OVERVIEW ──
if page == "Overview":
    st.title("Digenic Disease Explorer")

# ── AI INTERPRETATION (FIXED ONLY HERE) ──
elif page == "AI Interpretation":

    st.title("AI Interpretation")

    if ai_text:

        sections = [s for s in ai_text.split("="*50) if "RANK" in s]

        for section in sections:

            # FIXED: correct newline handling
            lines = section.strip().split("\n")

            title = next((l for l in lines if "RANK" in l), "Rank")
            score_line = next((l for l in lines if "Score" in l), "")

            content = "\n".join(lines[3:]).strip()

            score_val = 0
            try:
                score_val = float(score_line.split(":")[1].split("/")[0])
            except:
                pass

            color = score_color(score_val)

            with st.expander(f"{title} — {score_val:.0f}/100"):
                st.markdown(
                    f"<div style='border-left:4px solid {color};padding:10px;color:#cbd5e1;white-space:pre-wrap'>{content}</div>",
                    unsafe_allow_html=True
                )

    else:
        st.info("Run layer6_ai.py first")

# ── NETWORK ──
elif page == "Network":
    st.title("Network View")

# ── SCORES ──
elif page == "Digenic Scores":
    st.title("Scores")

# ── PIPELINE ──
elif page == "Pipeline":
    st.title("Pipeline")
