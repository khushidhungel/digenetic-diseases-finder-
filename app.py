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
import codecs   # ✅ FIX ADDED

st.set_page_config(
    page_title="BBS Digenic Explorer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ✅ FIX: sidebar visibility + width
st.markdown("""
<style>
[data-testid="stSidebar"] {
    min-width: 260px;
    max-width: 260px;
}
</style>
""", unsafe_allow_html=True)

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Space+Mono:wght@400;700&family=DM+Sans:wght@300;400;500;600&display=swap');
:root{--bg:#0a0e1a;--bg2:#111827;--bg3:#1a2235;--accent1:#00f5c4;--accent2:#6366f1;--accent3:#f59e0b;--accent4:#ef4444;--text:#f8fafc;--text2:#cbd5e1;--border:rgba(255,255,255,0.12);}
html,body,[class*="css"]{font-family:'DM Sans',sans-serif;background-color:var(--bg)!important;color:var(--text)!important;}
.stApp{background-color:var(--bg)!important;}
section[data-testid="stSidebar"]{background:var(--bg2)!important;border-right:1px solid var(--border);}
section[data-testid="stSidebar"] label, section[data-testid="stSidebar"] div, section[data-testid="stSidebar"] span, section[data-testid="stSidebar"] button { color: var(--text) !important; }
section[data-testid="stSidebar"] .stRadio, section[data-testid="stSidebar"] .stRadio label, section[data-testid="stSidebar"] .stRadio span { color: var(--text) !important; font-weight: 700 !important; }
.metric-box{background:var(--bg3);border:1px solid var(--border);border-radius:10px;padding:1rem;text-align:center;}
.metric-value{font-family:'Space Mono',monospace;font-size:2rem;font-weight:700;color:var(--accent1);line-height:1;}
.metric-label{font-size:0.75rem;color:var(--text);letter-spacing:0.05em;text-transform:uppercase;}
.hero-title{font-family:'Space Mono',monospace;font-size:2.4rem;font-weight:800;color:#ffffff;line-height:1.1;}
.section-title{font-family:'Space Mono',monospace;font-size:1rem;color:#ffffff;letter-spacing:0.16em;text-transform:uppercase;margin-bottom:0.75rem;}
.card{background:var(--bg3);border:1px solid rgba(255,255,255,0.12);border-radius:12px;padding:1.2rem 1.5rem;margin-bottom:1rem;}
[data-baseweb="tab-list"]{background:var(--bg2)!important;border-radius:12px;padding:4px;}
[data-baseweb="tab"]{background:transparent !important;color:var(--text) !important;font-weight:700 !important;border-radius:10px !important;}
[data-baseweb="tab"][aria-selected="true"]{background:#111827 !important;color:var(--accent1) !important;}
.stRadio [class*="label"]{color:var(--text) !important;}
#MainMenu{visibility:hidden;}footer{visibility:hidden;}header{visibility:hidden;}
</style>
""", unsafe_allow_html=True)

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

# ✅ FIX: safe UTF-8 reading
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

with st.sidebar:
    st.markdown("🧬 BBS Explorer")
    page = st.radio(
        "Navigation",   # ✅ FIX: label visible
        ["🏠 Overview","🔬 Gene Analysis","🕸️ Network","📊 Digenic Scores","🤖 AI Interpretation","⚙️ Pipeline"]
    )

def score_color(s):
    return "#ef4444" if s>=80 else "#f59e0b" if s>=60 else "#00f5c4" if s>=40 else "#6366f1"


# ── AI FIXED SECTION ──
elif "AI" in page:
    st.markdown("### AI Interpretation")

    if ai_text:
        sections = [s.strip() for s in ai_text.split("="*55) if s.strip() and "RANK" in s]

        for section in sections:
            lines = section.strip().splitlines()   # ✅ FIX HERE

            title = next((l for l in lines if "RANK" in l), "")
            score_ln = next((l for l in lines if "Score" in l or "Digenic Score" in l), "")

            content = "\n".join(lines[3:]).strip()

            score_val = 0
            try:
                score_val = float(score_ln.split(":")[1].split("/")[0])
            except:
                pass

            color = score_color(score_val)

            with st.expander(f"{title} — {score_val}/100", expanded=score_val>=90):
                st.markdown(
                    f"<div style='background:#111827;border-left:3px solid {color};padding:1rem;color:#cbd5e1;white-space:pre-wrap;'>{content}</div>",
                    unsafe_allow_html=True
                )
    else:
        st.info("Run AI pipeline first.")
