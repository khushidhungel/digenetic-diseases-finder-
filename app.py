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

st.set_page_config(page_title="BBS Digenic Explorer", page_icon="🧬", layout="wide", initial_sidebar_state="auto")

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Space+Mono:wght@400;700&family=DM+Sans:wght@300;400;500;600&display=swap');
:root{--bg:#0a0e1a;--bg2:#0d1f3c;--bg3:#1a2235;--accent1:#00f5c4;--accent2:#6366f1;--accent3:#f59e0b;--accent4:#ef4444;--text:#f8fafc;--text2:#cbd5e1;--border:rgba(255,255,255,0.12);}
html,body,[class*="css"]{font-family:'DM Sans',sans-serif;background-color:var(--bg)!important;color:var(--text)!important;}
.stApp{background-color:var(--bg)!important;}
section[data-testid="stSidebar"]{background:#0d1f3c !important;border-right:2px solid #00f5c4 !important;}
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
#MainMenu{visibility:hidden;}footer{visibility:hidden;}
header{visibility:hidden;}
[data-testid="collapsedControl"]{visibility:visible !important;display:flex !important;}
</style>
""", unsafe_allow_html=True)

@st.cache_data
def load_data():
    data = {}
    files = {
        "genes":       "layer1_bbs_genes.csv",
        "constrained": "layer2_constrained_genes.csv",
        "interactions":"layer3_interactions.csv",
        "pathways":    "layer4_pathway_scores.csv",
        "clinvar":     "layer4b_clinvar_scores.csv",
        "scores":      "layer5_digenic_scores.csv",
    }
    for key, path in files.items():
        data[key] = pd.read_csv(path) if os.path.exists(path) else pd.DataFrame()
    return data

def load_ai_text():
    path = "layer6_ai_interpretations.txt"
    if os.path.exists(path):
        with open(path, encoding="utf-8") as f:
            return f.read()
    return None

data    = load_data()
ai_text = load_ai_text()
scores_df = data["scores"]
genes_df  = data["constrained"]
inter_df  = data["interactions"]

with st.sidebar:
    st.markdown("<div style='padding:1rem 0;'><div style='font-family:Space Mono,monospace;font-size:1.1rem;color:#00f5c4;font-weight:700;'>🧬 BBS Explorer</div><div style='font-size:0.75rem;color:#94a3b8;margin-top:4px;'>Digenic Variant Architecture</div></div>", unsafe_allow_html=True)
    st.markdown("---")
    page = st.radio("", ["🏠  Overview","🔬  Gene Analysis","🕸️  Network","📊  Digenic Scores","🤖  AI Interpretation","⚙️  Pipeline"], label_visibility="collapsed")
    st.markdown("---")
    st.markdown("<div class='card' style='padding:0.8rem;'><div style='font-size:0.7rem;color:#94a3b8;text-transform:uppercase;letter-spacing:0.05em;margin-bottom:8px;'>Disease</div><div style='font-size:0.9rem;color:#e2e8f0;font-weight:500;'>Bardet-Biedl Syndrome</div><div style='font-size:0.75rem;color:#6366f1;margin-top:4px;'>MONDO:0015229</div></div>", unsafe_allow_html=True)
    st.markdown("<div style='font-size:0.7rem;color:#475569;text-align:center;margin-top:1rem;'>Biohackathon 2026 · KU Bioinformatics</div>", unsafe_allow_html=True)

def score_color(s):
    return "#ef4444" if s>=80 else "#f59e0b" if s>=60 else "#00f5c4" if s>=40 else "#6366f1"

# ── OVERVIEW ──
if "Overview" in page:
    st.markdown('<div class="hero-title">Digenic Disease<br>Architecture Explorer</div>', unsafe_allow_html=True)
    st.markdown('<div style="color:#94a3b8;font-size:1rem;margin:0.5rem 0 1.5rem;">Cross-database variant prioritization · Bardet-Biedl Syndrome</div>', unsafe_allow_html=True)
    c1,c2,c3,c4 = st.columns(4)
    with c1: st.markdown(f"<div class='metric-box'><div class='metric-value'>{len(genes_df) if not genes_df.empty else 9}</div><div class='metric-label'>BBS Genes</div></div>", unsafe_allow_html=True)
    with c2: st.markdown(f"<div class='metric-box'><div class='metric-value' style='color:#6366f1;'>{len(inter_df) if not inter_df.empty else 16}</div><div class='metric-label'>Interactions</div></div>", unsafe_allow_html=True)
    with c3: st.markdown(f"<div class='metric-box'><div class='metric-value' style='color:#f59e0b;'>{len(scores_df) if not scores_df.empty else 15}</div><div class='metric-label'>Digenic Pairs</div></div>", unsafe_allow_html=True)
    with c4:
        top = round(scores_df["digenic_score_normalized"].max(),1) if not scores_df.empty else 100.0
        st.markdown(f"<div class='metric-box'><div class='metric-value' style='color:#ef4444;'>{top}</div><div class='metric-label'>Top Score</div></div>", unsafe_allow_html=True)
    st.markdown("<br>", unsafe_allow_html=True)
    col_a,col_b = st.columns([1.2,1])
    with col_a:
        st.markdown('<div class="section-title">Top Digenic Candidates</div>', unsafe_allow_html=True)
        if not scores_df.empty:
            for i,row in scores_df.head(5).iterrows():
                s=row["digenic_score_normalized"]; c=score_color(s)
                st.markdown(f"<div class='card' style='padding:0.8rem 1rem;margin-bottom:0.5rem;'><div style='display:flex;justify-content:space-between;align-items:center;margin-bottom:6px;'><span style='font-family:Space Mono,monospace;font-size:0.9rem;color:#e2e8f0;'>{row['gene_a']} <span style='color:#475569;'>↔</span> {row['gene_b']}</span><span style='font-family:Space Mono,monospace;font-size:0.85rem;color:{c};font-weight:700;'>{s:.1f}/100</span></div><div style='background:#1e293b;border-radius:4px;height:4px;'><div style='background:{c};width:{int(s)}%;height:100%;border-radius:4px;'></div></div><div style='display:flex;gap:12px;margin-top:6px;'><span style='font-size:0.7rem;color:#64748b;'>STRING: {row['string_score']}</span><span style='font-size:0.7rem;color:#64748b;'>Pathways: {row['shared_pathways']}</span><span style='font-size:0.7rem;color:#64748b;'>ClinVar: {row['clinvar_pair']:.2f}</span></div></div>", unsafe_allow_html=True)
    with col_b:
        st.markdown('<div class="section-title">Databases Used</div>', unsafe_allow_html=True)
        for icon,name,desc,color in [("🎯","Open Targets","Gene harvesting","#00f5c4"),("🧬","gnomAD v4","Constraint filtering","#6366f1"),("🕸️","STRING DB","PPI network","#f59e0b"),("🔬","Reactome","Pathway co-membership","#00f5c4"),("🏥","ClinVar","Clinical evidence","#ef4444"),("🤖","Gemini AI","Clinical interpretation","#6366f1")]:
            st.markdown(f"<div style='display:flex;align-items:center;gap:10px;padding:8px 10px;background:#111827;border-radius:8px;margin-bottom:4px;border:1px solid rgba(255,255,255,0.05);'><span style='font-size:1.1rem;'>{icon}</span><div><div style='font-size:0.85rem;color:#e2e8f0;font-weight:500;'>{name}</div><div style='font-size:0.7rem;color:#64748b;'>{desc}</div></div><div style='margin-left:auto;width:6px;height:6px;border-radius:50%;background:{color};'></div></div>", unsafe_allow_html=True)

# ── GENE ANALYSIS ──
elif "Gene" in page:
    st.markdown('<div class="section-title">Gene Analysis</div>', unsafe_allow_html=True)
    st.markdown("### BBS Gene Constraint Profile")
    if not genes_df.empty:
        col1,col2 = st.columns([1.5,1])
        with col1:
            colors = ["#ef4444" if float(p)>=0.99 else "#f59e0b" if float(p)>=0.95 else "#00f5c4" for p in genes_df["pLI"]]
            fig = go.Figure(go.Bar(x=genes_df["pLI"].astype(float),y=genes_df["gene_symbol"],orientation="h",marker_color=colors,text=[f"{float(p):.3f}" for p in genes_df["pLI"]],textposition="outside",hovertemplate="<b>%{y}</b><br>pLI: %{x:.4f}<extra></extra>"))
            fig.add_vline(x=0.9,line_dash="dash",line_color="#475569",annotation_text="cutoff 0.9",annotation_font_color="#94a3b8")
            fig.update_layout(plot_bgcolor="#111827",paper_bgcolor="#0a0e1a",font_color="#94a3b8",height=380,xaxis=dict(range=[0,1.08],gridcolor="#1e293b",title="pLI Score"),yaxis=dict(autorange="reversed",gridcolor="#1e293b"),title=dict(text="gnomAD pLI Constraint Scores",font_color="#e2e8f0"),margin=dict(l=10,r=60,t=40,b=10))
            st.plotly_chart(fig,use_container_width=True)
        with col2:
            st.markdown('<div class="section-title">Gene Details</div>', unsafe_allow_html=True)
            for _,row in genes_df.iterrows():
                pli=float(row["pLI"]); color="#ef4444" if pli>=0.99 else "#f59e0b" if pli>=0.95 else "#00f5c4"
                st.markdown(f"<div class='card' style='padding:0.7rem;margin-bottom:0.4rem;'><div style='display:flex;justify-content:space-between;'><span style='font-family:Space Mono,monospace;font-size:0.85rem;color:#e2e8f0;'>{row['gene_symbol']}</span><span style='font-family:Space Mono,monospace;font-size:0.8rem;color:{color};'>pLI={pli:.3f}</span></div><div style='font-size:0.7rem;color:#64748b;margin-top:2px;'>OT score: {row['ot_score']:.3f}</div></div>", unsafe_allow_html=True)

# ── NETWORK ──
elif "Network" in page:
    st.markdown('<div class="section-title">PPI Network</div>', unsafe_allow_html=True)
    st.markdown("### Interactive Protein-Protein Interaction Network")
    if not inter_df.empty and not genes_df.empty:
        all_genes = sorted(genes_df["gene_symbol"].tolist())
        selected  = st.multiselect("🔍 Type a gene name to filter interactions:", options=all_genes, default=all_genes)
        if not selected: selected = all_genes
        filtered_inter = inter_df[inter_df["gene_a"].isin(selected) & inter_df["gene_b"].isin(selected)]
        G = nx.Graph()
        pli_map = dict(zip(genes_df["gene_symbol"],genes_df["pLI"].astype(float)))
        for gene in selected: G.add_node(gene,pLI=pli_map.get(gene,0.5))
        for _,row in filtered_inter.iterrows(): G.add_edge(row["gene_a"],row["gene_b"],weight=float(row["score"]))
        pos = nx.spring_layout(G,seed=42,k=2.5)
        edge_x,edge_y=[],[]
        for u,v in G.edges():
            x0,y0=pos[u]; x1,y1=pos[v]
            edge_x+=[x0,x1,None]; edge_y+=[y0,y1,None]
        node_color=["#ef4444" if pli_map.get(n,0)>=0.99 else "#f59e0b" if pli_map.get(n,0)>=0.95 else "#00f5c4" for n in G.nodes()]
        node_size=[20+G.degree(n)*6 for n in G.nodes()]
        fig=go.Figure()
        fig.add_trace(go.Scatter(x=edge_x,y=edge_y,mode="lines",line=dict(width=1.5,color="#6366f1"),opacity=0.5,hoverinfo="none"))
        fig.add_trace(go.Scatter(x=[pos[n][0] for n in G.nodes()],y=[pos[n][1] for n in G.nodes()],mode="markers+text",marker=dict(size=node_size,color=node_color,line=dict(width=1.5,color="#0a0e1a")),text=list(G.nodes()),textposition="top center",textfont=dict(color="#e2e8f0",size=11),hovertext=[f"<b>{n}</b><br>pLI: {pli_map.get(n,0):.3f}<br>Connections: {G.degree(n)}" for n in G.nodes()],hoverinfo="text"))
        fig.update_layout(plot_bgcolor="#0a0e1a",paper_bgcolor="#0a0e1a",height=550,showlegend=False,xaxis=dict(showgrid=False,zeroline=False,showticklabels=False),yaxis=dict(showgrid=False,zeroline=False,showticklabels=False),margin=dict(l=20,r=20,t=40,b=20),title=dict(text=f"BBS PPI Network — {len(selected)} genes, {G.number_of_edges()} interactions",font_color="#e2e8f0"))
        st.plotly_chart(fig,use_container_width=True)
        st.markdown('<div class="section-title">Hub Genes</div>', unsafe_allow_html=True)
        top_hubs = sorted(G.degree(),key=lambda x:x[1],reverse=True)[:min(5,len(selected))]
        cols = st.columns(len(top_hubs))
        for i,(gene,deg) in enumerate(top_hubs):
            with cols[i]: st.markdown(f"<div class='metric-box'><div class='metric-value' style='font-size:1.2rem;color:#00f5c4;'>{gene}</div><div class='metric-label'>{deg} connections</div></div>", unsafe_allow_html=True)

# ── DIGENIC SCORES ──
elif "Digenic" in page:
    st.markdown('<div class="section-title">Digenic Scores</div>', unsafe_allow_html=True)
    st.markdown("### Final Ranked Digenic Candidate Pairs")
    if not scores_df.empty:
        col_f1,col_f2 = st.columns(2)
        with col_f1:
            all_genes   = sorted(set(scores_df["gene_a"].tolist()+scores_df["gene_b"].tolist()))
            gene_filter = st.multiselect("🔍 Filter by gene:",all_genes,default=[])
        with col_f2:
            min_score = st.slider("Minimum score:",0,100,0,5)
        filtered = scores_df.copy()
        if gene_filter: filtered = filtered[filtered["gene_a"].isin(gene_filter)|filtered["gene_b"].isin(gene_filter)]
        filtered = filtered[filtered["digenic_score_normalized"]>=min_score]
        st.markdown(f"<div style='font-size:0.8rem;color:#64748b;margin-bottom:12px;'>Showing {len(filtered)} pairs</div>", unsafe_allow_html=True)
        tab1,tab2,tab3 = st.tabs(["📊 Bar Chart","🔥 Heatmap","📋 Table"])
        with tab1:
            top=filtered.head(15).copy(); top["label"]=top["gene_a"]+" ↔ "+top["gene_b"]
            colors=[score_color(s) for s in top["digenic_score_normalized"]]
            fig=go.Figure(go.Bar(x=top["digenic_score_normalized"],y=top["label"],orientation="h",marker_color=colors,text=[f"{s:.1f}" for s in top["digenic_score_normalized"]],textposition="outside",hovertemplate="<b>%{y}</b><br>Score: %{x:.1f}/100<extra></extra>"))
            fig.update_layout(plot_bgcolor="#111827",paper_bgcolor="#0a0e1a",font_color="#94a3b8",height=max(400,len(top)*40),xaxis=dict(range=[0,115],gridcolor="#1e293b",title="Digenic Score"),yaxis=dict(autorange="reversed",gridcolor="#1e293b"),title=dict(text="Digenic Score Ranking",font_color="#e2e8f0"),margin=dict(l=10,r=60,t=40,b=10))
            st.plotly_chart(fig,use_container_width=True)
        with tab2:
            gene_list=genes_df["gene_symbol"].tolist() if not genes_df.empty else []
            matrix=pd.DataFrame(0.0,index=gene_list,columns=gene_list)
            for _,row in filtered.iterrows():
                a,b,s=row["gene_a"],row["gene_b"],row["digenic_score_normalized"]
                if a in matrix.index and b in matrix.columns: matrix.loc[a,b]=s; matrix.loc[b,a]=s
            fig=px.imshow(matrix,color_continuous_scale="YlOrRd",text_auto=".0f",aspect="auto",title="Digenic Score Heatmap")
            fig.update_layout(plot_bgcolor="#111827",paper_bgcolor="#0a0e1a",font_color="#94a3b8",height=450,title_font_color="#e2e8f0")
            fig.update_traces(textfont_size=10)
            st.plotly_chart(fig,use_container_width=True)
        with tab3:
            display=filtered[["gene_a","gene_b","digenic_score_normalized","shared_pathways","string_score","clinvar_pair"]].copy()
            display.columns=["Gene A","Gene B","Score","Pathways","STRING","ClinVar"]
            display["Score"]=display["Score"].round(1)
            st.dataframe(display.style.background_gradient(subset=["Score"],cmap="YlOrRd"),use_container_width=True,height=400)

# ── AI INTERPRETATION ──
elif "AI" in page:
    st.markdown('<div class="section-title">AI Interpretation</div>', unsafe_allow_html=True)
    st.markdown("### Gemini-Generated Clinical Briefs")
    if ai_text:
        sections=[s.strip() for s in ai_text.split("="*55) if s.strip() and "RANK" in s]
        for section in sections:
            lines = section.strip().split("\n")
            title = next((l for l in lines if "RANK" in l), "")
            score_ln = next((l for l in lines if "Digenic Score" in l), "")
            content = "\n".join(lines[3:]).strip()
            score_val=0
            try: score_val=float(score_ln.split(":")[1].strip().split("/")[0])
            except: pass
            color=score_color(score_val)
            with st.expander(f"🧬 {title} — Score: {score_val:.0f}/100",expanded=score_val>=90):
                st.markdown(f"<div style='background:#111827;border-left:3px solid {color};padding:1rem;border-radius:0 8px 8px 0;font-size:0.875rem;line-height:1.7;color:#cbd5e1;white-space:pre-wrap;'>{content}</div>", unsafe_allow_html=True)
    else:
        st.info("Run layer6_ai.py first to generate AI interpretations.")

# ── PIPELINE ──
elif "Pipeline" in page:
    st.markdown('<div class="section-title">Pipeline</div>', unsafe_allow_html=True)
    st.markdown("### 6-Layer Analysis Pipeline")
    for num,name,source,outfile,color in [
        ("1",  "Gene Harvesting",    "Open Targets GraphQL",   "layer1_bbs_genes.csv",        "#00f5c4"),
        ("2",  "Constraint Filtering","gnomAD pLI scores",     "layer2_constrained_genes.csv","#6366f1"),
        ("3",  "PPI Network",         "STRING DB interactions", "layer3_interactions.csv",     "#f59e0b"),
        ("4",  "Pathway Analysis",    "Reactome co-membership", "layer4_pathway_scores.csv",   "#00f5c4"),
        ("4b", "Clinical Evidence",   "ClinVar pathogenicity",  "layer4b_clinvar_scores.csv",  "#ef4444"),
        ("5",  "Digenic Scoring",     "DiVaS-inspired formula", "layer5_digenic_scores.csv",   "#f59e0b"),
        ("6",  "AI Interpretation",   "Gemini 2.0 Flash",       "layer6_ai_summary.csv",       "#6366f1"),
    ]:
        exists="✅" if os.path.exists(outfile) else "⏳"
        st.markdown(f"<div class='card' style='display:flex;align-items:center;gap:16px;padding:0.9rem 1.2rem;margin-bottom:0.5rem;border-left:3px solid {color};'><div style='font-family:Space Mono,monospace;font-size:1rem;color:{color};font-weight:700;min-width:32px;'>L{num}</div><div style='flex:1;'><div style='font-size:0.9rem;color:#e2e8f0;font-weight:500;'>{name}</div><div style='font-size:0.75rem;color:#64748b;margin-top:2px;'>{source}</div></div><div style='font-size:1.2rem;'>{exists}</div></div>", unsafe_allow_html=True)
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown('<div class="section-title">Digenic Score Formula</div>', unsafe_allow_html=True)
    st.code("DS = (string_conf × pLI_A × pLI_B × pathway_overlap × clinvar_pair)
     ÷ (AF_A × AF_B × 1e6)

Normalized to 0–100  |  Inspired by DiVaS algorithm", language="text")
