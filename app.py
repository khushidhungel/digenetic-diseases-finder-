"""
BBS Digenic Variant Explorer — Complete Final App
Run: streamlit run app_final.py
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

/* Force sidebar always open */
section[data-testid="stSidebar"] {
    min-width: 260px !important;
    max-width: 260px !important;
    transform: none !important;
    visibility: visible !important;
}
[data-testid="collapsedControl"] { display: none !important; }

html, body, [class*="css"] {
    font-family: 'DM Sans', sans-serif;
    background-color: #0a0e1a !important;
    color: #e2e8f0 !important;
}
.stApp { background-color: #0a0e1a !important; }
section[data-testid="stSidebar"] {
    background: #111827 !important;
    border-right: 1px solid rgba(255,255,255,0.08);
}
.metric-box {
    background: #1a2235;
    border: 1px solid rgba(255,255,255,0.08);
    border-radius: 10px;
    padding: 1rem;
    text-align: center;
    margin-bottom: 0.5rem;
}
.metric-value {
    font-family: 'Space Mono', monospace;
    font-size: 2rem;
    font-weight: 700;
    color: #00f5c4;
    line-height: 1;
}
.metric-label {
    font-size: 0.7rem;
    color: #94a3b8;
    margin-top: 4px;
    letter-spacing: 0.08em;
    text-transform: uppercase;
}
.hero-title {
    font-family: 'Space Mono', monospace;
    font-size: 2rem;
    font-weight: 700;
    background: linear-gradient(135deg, #00f5c4, #6366f1);
    -webkit-background-clip: text;
    -webkit-text-fill-color: transparent;
    line-height: 1.3;
    margin-bottom: 0.3rem;
}
.section-title {
    font-family: 'Space Mono', monospace;
    font-size: 0.8rem;
    color: #00f5c4;
    letter-spacing: 0.12em;
    text-transform: uppercase;
    margin-bottom: 0.6rem;
}
.card {
    background: #1a2235;
    border: 1px solid rgba(255,255,255,0.08);
    border-radius: 10px;
    padding: 1rem 1.2rem;
    margin-bottom: 0.6rem;
}
.stTabs [data-baseweb="tab-list"] {
    background: #111827 !important;
    border-radius: 8px;
    padding: 4px;
    gap: 4px;
}
.stTabs [data-baseweb="tab"] {
    background: transparent !important;
    color: #94a3b8 !important;
    border-radius: 6px !important;
}
.stTabs [aria-selected="true"] {
    background: #1a2235 !important;
    color: #00f5c4 !important;
}
div[data-testid="stRadio"] label {
    color: #94a3b8 !important;
    font-size: 0.9rem;
    padding: 6px 10px;
    border-radius: 6px;
    transition: all 0.15s;
}
div[data-testid="stRadio"] label:hover { background: #1a2235 !important; color: #e2e8f0 !important; }
#MainMenu { visibility: hidden; }
footer { visibility: hidden; }
header { visibility: hidden; }
</style>
""", unsafe_allow_html=True)

# ── Data loading ───────────────────────────────────
@st.cache_data
def load_data():
    base  = os.path.dirname(os.path.abspath(__file__))
    files = {
        "genes":       "outputs/layer1_bbs_genes.csv",
        "constrained": "outputs/layer2_constrained_genes.csv",
        "interactions":"outputs/layer3_interactions.csv",
        "pathways":    "outputs/layer4_pathway_scores.csv",
        "clinvar":     "outputs/layer4b_clinvar_scores.csv",
        "scores":      "outputs/layer5_digenic_scores.csv",
    }
    data = {}
    for key, rel in files.items():
        path = os.path.join(base, rel)
        data[key] = pd.read_csv(path) if os.path.exists(path) else pd.DataFrame()
    return data

@st.cache_data
def load_ai():
    base = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(base, "outputs/layer6_ai_interpretations.txt")
    return open(path, encoding="utf-8").read() if os.path.exists(path) else None

data      = load_data()
ai_text   = load_ai()
scores_df = data["scores"]
genes_df  = data["constrained"]
inter_df  = data["interactions"]

def score_color(s):
    return "#ef4444" if s>=80 else "#f59e0b" if s>=60 else "#00f5c4" if s>=40 else "#6366f1"

# ── Sidebar ────────────────────────────────────────
with st.sidebar:
    st.markdown("""
    <div style='padding:0.8rem 0 1rem;'>
        <div style='font-family:Space Mono,monospace;font-size:1.1rem;color:#00f5c4;font-weight:700;letter-spacing:0.05em;'>🧬 BBS Explorer</div>
        <div style='font-size:0.72rem;color:#64748b;margin-top:3px;'>Digenic Variant Architecture</div>
    </div>
    """, unsafe_allow_html=True)

    st.markdown("---")
    page = st.radio("Navigation", [
        "🏠  Overview",
        "🔬  Gene Analysis",
        "🕸️  Network",
        "📊  Digenic Scores",
        "🤖  AI Interpretation",
        "⚙️  Pipeline"
    ], label_visibility="collapsed")

    st.markdown("---")

    # Quick stats in sidebar
    if not scores_df.empty:
        top_pair = scores_df.iloc[0]
        st.markdown(f"""
        <div class='card' style='padding:0.7rem;'>
            <div style='font-size:0.65rem;color:#64748b;text-transform:uppercase;letter-spacing:0.08em;margin-bottom:6px;'>⭐ Top Finding</div>
            <div style='font-family:Space Mono,monospace;font-size:0.85rem;color:#00f5c4;'>{top_pair['gene_a']} ↔ {top_pair['gene_b']}</div>
            <div style='font-size:0.75rem;color:#ef4444;margin-top:2px;font-weight:600;'>Score: {top_pair['digenic_score_normalized']:.1f}/100</div>
        </div>
        """, unsafe_allow_html=True)

    st.markdown("""
    <div class='card' style='padding:0.7rem;'>
        <div style='font-size:0.65rem;color:#64748b;text-transform:uppercase;letter-spacing:0.08em;margin-bottom:6px;'>Disease</div>
        <div style='font-size:0.85rem;color:#e2e8f0;font-weight:500;'>Bardet-Biedl Syndrome</div>
        <div style='font-size:0.7rem;color:#6366f1;margin-top:2px;'>MONDO:0015229</div>
    </div>
    """, unsafe_allow_html=True)

    st.markdown("<div style='font-size:0.65rem;color:#334155;text-align:center;padding-top:0.5rem;'>Biohackathon 2026 · KU Bioinformatics</div>", unsafe_allow_html=True)

# ══════════════════════════════════════════════════
# OVERVIEW
# ══════════════════════════════════════════════════
if "Overview" in page:
    st.markdown('<div class="hero-title">Digenic Disease<br>Architecture Explorer</div>', unsafe_allow_html=True)
    st.markdown('<div style="color:#64748b;font-size:0.95rem;margin-bottom:1.5rem;">Cross-database variant prioritization · Bardet-Biedl Syndrome</div>', unsafe_allow_html=True)

    c1,c2,c3,c4 = st.columns(4)
    ng = len(genes_df) if not genes_df.empty else 9
    ni = len(inter_df) if not inter_df.empty else 16
    np_ = len(scores_df) if not scores_df.empty else 15
    ts = round(scores_df["digenic_score_normalized"].max(),1) if not scores_df.empty else 100.0

    with c1: st.markdown(f"<div class='metric-box'><div class='metric-value'>{ng}</div><div class='metric-label'>BBS Genes</div></div>", unsafe_allow_html=True)
    with c2: st.markdown(f"<div class='metric-box'><div class='metric-value' style='color:#6366f1;'>{ni}</div><div class='metric-label'>Interactions</div></div>", unsafe_allow_html=True)
    with c3: st.markdown(f"<div class='metric-box'><div class='metric-value' style='color:#f59e0b;'>{np_}</div><div class='metric-label'>Digenic Pairs</div></div>", unsafe_allow_html=True)
    with c4: st.markdown(f"<div class='metric-box'><div class='metric-value' style='color:#ef4444;'>{ts}</div><div class='metric-label'>Top Score</div></div>", unsafe_allow_html=True)

    st.markdown("<br>", unsafe_allow_html=True)
    col_a, col_b = st.columns([1.3, 1])

    with col_a:
        st.markdown('<div class="section-title">Top Digenic Candidates</div>', unsafe_allow_html=True)
        if not scores_df.empty:
            for _, row in scores_df.head(5).iterrows():
                s = float(row["digenic_score_normalized"])
                c = score_color(s)
                st.markdown(f"""
                <div class='card' style='padding:0.8rem 1rem;margin-bottom:0.4rem;'>
                    <div style='display:flex;justify-content:space-between;align-items:center;margin-bottom:6px;'>
                        <span style='font-family:Space Mono,monospace;font-size:0.88rem;color:#e2e8f0;'>{row['gene_a']} <span style='color:#334155;'>↔</span> {row['gene_b']}</span>
                        <span style='font-family:Space Mono,monospace;font-size:0.82rem;color:{c};font-weight:700;'>{s:.1f}/100</span>
                    </div>
                    <div style='background:#0d1320;border-radius:4px;height:5px;'>
                        <div style='background:{c};width:{int(s)}%;height:100%;border-radius:4px;'></div>
                    </div>
                    <div style='display:flex;gap:14px;margin-top:6px;'>
                        <span style='font-size:0.7rem;color:#475569;'>STRING: {row["string_score"]}</span>
                        <span style='font-size:0.7rem;color:#475569;'>Pathways: {row["shared_pathways"]}</span>
                        <span style='font-size:0.7rem;color:#475569;'>ClinVar: {float(row["clinvar_pair"]):.2f}</span>
                    </div>
                </div>""", unsafe_allow_html=True)
        else:
            st.info("Run the pipeline first to see results")

    with col_b:
        st.markdown('<div class="section-title">Databases Used</div>', unsafe_allow_html=True)
        for icon, name, desc, color in [
            ("🎯","Open Targets","Gene harvesting","#00f5c4"),
            ("🧬","gnomAD v4","Constraint filtering","#6366f1"),
            ("🕸️","STRING DB","PPI network","#f59e0b"),
            ("🔬","Reactome","Pathway co-membership","#00f5c4"),
            ("🏥","ClinVar","Clinical evidence","#ef4444"),
            ("🤖","Gemini AI","Clinical interpretation","#6366f1"),
        ]:
            st.markdown(f"""
            <div style='display:flex;align-items:center;gap:10px;padding:8px 10px;background:#111827;border-radius:8px;margin-bottom:4px;border:1px solid rgba(255,255,255,0.05);'>
                <span style='font-size:1rem;'>{icon}</span>
                <div style='flex:1;'>
                    <div style='font-size:0.82rem;color:#e2e8f0;font-weight:500;'>{name}</div>
                    <div style='font-size:0.68rem;color:#475569;'>{desc}</div>
                </div>
                <div style='width:6px;height:6px;border-radius:50%;background:{color};flex-shrink:0;'></div>
            </div>""", unsafe_allow_html=True)

# ══════════════════════════════════════════════════
# GENE ANALYSIS
# ══════════════════════════════════════════════════
elif "Gene" in page:
    st.markdown('<div class="section-title">Gene Analysis</div>', unsafe_allow_html=True)
    st.markdown("### BBS Gene Constraint Profile")

    if not genes_df.empty:
        col1, col2 = st.columns([1.6, 1])
        with col1:
            pli_vals = genes_df["pLI"].astype(float).tolist()
            colors   = ["#ef4444" if p>=0.99 else "#f59e0b" if p>=0.95 else "#00f5c4" for p in pli_vals]
            fig = go.Figure(go.Bar(
                x=pli_vals, y=genes_df["gene_symbol"].tolist(),
                orientation="h", marker_color=colors,
                text=[f"{p:.3f}" for p in pli_vals], textposition="outside",
                hovertemplate="<b>%{y}</b><br>pLI: %{x:.4f}<br>Higher = more constrained<extra></extra>"
            ))
            fig.add_vline(x=0.9, line_dash="dash", line_color="#475569",
                         annotation_text="cutoff", annotation_font_color="#64748b")
            fig.update_layout(
                plot_bgcolor="#111827", paper_bgcolor="#0a0e1a",
                font_color="#94a3b8", height=400,
                xaxis=dict(range=[0,1.1], gridcolor="#1e293b", title="pLI Score", color="#64748b"),
                yaxis=dict(autorange="reversed", gridcolor="#1e293b", color="#94a3b8"),
                title=dict(text="gnomAD pLI Scores — gene intolerance to mutations", font_color="#e2e8f0", font_size=13),
                margin=dict(l=10, r=70, t=45, b=10),
                hoverlabel=dict(bgcolor="#1a2235", font_color="#e2e8f0")
            )
            st.plotly_chart(fig, use_container_width=True)

            st.markdown("""
            <div style='display:flex;gap:16px;margin-top:4px;flex-wrap:wrap;'>
                <span style='font-size:0.72rem;color:#94a3b8;'><span style='color:#ef4444;'>■</span> pLI≥0.99 extremely constrained</span>
                <span style='font-size:0.72rem;color:#94a3b8;'><span style='color:#f59e0b;'>■</span> pLI≥0.95 highly constrained</span>
                <span style='font-size:0.72rem;color:#94a3b8;'><span style='color:#00f5c4;'>■</span> pLI<0.95 moderate</span>
            </div>""", unsafe_allow_html=True)

        with col2:
            st.markdown('<div class="section-title">Gene Details</div>', unsafe_allow_html=True)
            for _, row in genes_df.iterrows():
                pli   = float(row["pLI"])
                color = "#ef4444" if pli>=0.99 else "#f59e0b" if pli>=0.95 else "#00f5c4"
                st.markdown(f"""
                <div class='card' style='padding:0.65rem 0.9rem;margin-bottom:0.35rem;'>
                    <div style='display:flex;justify-content:space-between;align-items:center;'>
                        <span style='font-family:Space Mono,monospace;font-size:0.82rem;color:#e2e8f0;'>{row['gene_symbol']}</span>
                        <span style='font-size:0.75rem;color:{color};font-weight:600;'>pLI={pli:.3f}</span>
                    </div>
                    <div style='display:flex;justify-content:space-between;margin-top:3px;'>
                        <span style='font-size:0.68rem;color:#475569;'>OT: {float(row['ot_score']):.3f}</span>
                        <span style='font-size:0.68rem;color:{color};'>{"★★★" if pli>=0.99 else "★★" if pli>=0.95 else "★"}</span>
                    </div>
                </div>""", unsafe_allow_html=True)
    else:
        st.warning("Run layer2_gnomad.py first")

# ══════════════════════════════════════════════════
# NETWORK
# ══════════════════════════════════════════════════
elif "Network" in page:
    st.markdown('<div class="section-title">PPI Network</div>', unsafe_allow_html=True)
    st.markdown("### Interactive Protein-Protein Interaction Network")

    if not inter_df.empty and not genes_df.empty:
        all_genes = sorted(genes_df["gene_symbol"].tolist())

        col_f1, col_f2 = st.columns([2, 1])
        with col_f1:
            selected = st.multiselect(
                "🔍 Type gene name to filter network:",
                options=all_genes, default=all_genes,
                help="Remove genes to focus on specific interactions"
            )
        with col_f2:
            min_edge = st.slider("Min interaction score:", 500, 999, 700, 50)

        if not selected: selected = all_genes

        filtered_inter = inter_df[
            inter_df["gene_a"].isin(selected) &
            inter_df["gene_b"].isin(selected) &
            (inter_df["score"] >= min_edge)
        ]

        G       = nx.Graph()
        pli_map = dict(zip(genes_df["gene_symbol"], genes_df["pLI"].astype(float)))
        for gene in selected:
            G.add_node(gene, pLI=pli_map.get(gene, 0.5))
        for _, row in filtered_inter.iterrows():
            G.add_edge(row["gene_a"], row["gene_b"], weight=float(row["score"]))

        pos = nx.spring_layout(G, seed=42, k=2.8)

        # Edges
        edge_x, edge_y, edge_hover = [], [], []
        for u, v, d in G.edges(data=True):
            x0,y0=pos[u]; x1,y1=pos[v]
            edge_x+=[x0,x1,None]; edge_y+=[y0,y1,None]
            edge_hover.append(f"{u} ↔ {v}: {d['weight']}")

        node_x     = [pos[n][0] for n in G.nodes()]
        node_y     = [pos[n][1] for n in G.nodes()]
        node_names = list(G.nodes())
        node_color = ["#ef4444" if pli_map.get(n,0)>=0.99 else "#f59e0b" if pli_map.get(n,0)>=0.95 else "#00f5c4" for n in G.nodes()]
        node_size  = [22+G.degree(n)*7 for n in G.nodes()]
        node_hover = [f"<b>{n}</b><br>pLI: {pli_map.get(n,0):.3f}<br>Connections: {G.degree(n)}<br>Click to explore" for n in G.nodes()]

        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=edge_x, y=edge_y, mode="lines",
            line=dict(width=1.8, color="#6366f1"), opacity=0.45, hoverinfo="none"
        ))
        fig.add_trace(go.Scatter(
            x=node_x, y=node_y, mode="markers+text",
            marker=dict(size=node_size, color=node_color, line=dict(width=2, color="#0a0e1a"),
                       opacity=0.95),
            text=node_names, textposition="top center",
            textfont=dict(color="#e2e8f0", size=11, family="Space Mono"),
            hovertext=node_hover, hoverinfo="text",
            hoverlabel=dict(bgcolor="#1a2235", font_color="#e2e8f0")
        ))
        fig.update_layout(
            plot_bgcolor="#0a0e1a", paper_bgcolor="#0a0e1a",
            height=520, showlegend=False,
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            margin=dict(l=20, r=20, t=45, b=20),
            title=dict(
                text=f"BBS PPI Network — {len(selected)} genes · {G.number_of_edges()} interactions · score≥{min_edge}",
                font_color="#e2e8f0", font_size=13
            )
        )
        st.plotly_chart(fig, use_container_width=True)

        # Hub genes
        if G.number_of_nodes() > 0:
            st.markdown('<div class="section-title">Hub Genes (most connected)</div>', unsafe_allow_html=True)
            top_hubs = sorted(G.degree(), key=lambda x: x[1], reverse=True)[:min(5, len(selected))]
            cols     = st.columns(len(top_hubs))
            for i, (gene, deg) in enumerate(top_hubs):
                with cols[i]:
                    pli = pli_map.get(gene, 0)
                    c   = "#ef4444" if pli>=0.99 else "#f59e0b" if pli>=0.95 else "#00f5c4"
                    st.markdown(f"""
                    <div class='metric-box'>
                        <div class='metric-value' style='font-size:1.1rem;color:{c};'>{gene}</div>
                        <div class='metric-label'>{deg} connections</div>
                        <div style='font-size:0.65rem;color:#475569;margin-top:2px;'>pLI={pli:.3f}</div>
                    </div>""", unsafe_allow_html=True)
    else:
        st.warning("Run layer3_string.py first")

# ══════════════════════════════════════════════════
# DIGENIC SCORES
# ══════════════════════════════════════════════════
elif "Digenic" in page:
    st.markdown('<div class="section-title">Digenic Scores</div>', unsafe_allow_html=True)
    st.markdown("### Final Ranked Digenic Candidate Pairs")

    if not scores_df.empty:
        col_f1, col_f2 = st.columns(2)
        with col_f1:
            all_genes   = sorted(set(scores_df["gene_a"].tolist()+scores_df["gene_b"].tolist()))
            gene_filter = st.multiselect("🔍 Filter by gene:", all_genes, default=[],
                                        placeholder="Type gene name...")
        with col_f2:
            min_score = st.slider("Minimum digenic score:", 0, 100, 0, 5)

        filtered = scores_df.copy()
        if gene_filter:
            filtered = filtered[filtered["gene_a"].isin(gene_filter)|filtered["gene_b"].isin(gene_filter)]
        filtered = filtered[filtered["digenic_score_normalized"]>=min_score]

        st.markdown(f"<div style='font-size:0.78rem;color:#475569;margin-bottom:10px;'>Showing {len(filtered)} of {len(scores_df)} pairs</div>", unsafe_allow_html=True)

        tab1, tab2, tab3 = st.tabs(["📊  Bar Chart", "🔥  Heatmap", "📋  Full Table"])

        with tab1:
            top = filtered.head(15).copy()
            top["label"] = top["gene_a"] + " ↔ " + top["gene_b"]
            colors = [score_color(float(s)) for s in top["digenic_score_normalized"]]
            fig = go.Figure(go.Bar(
                x=top["digenic_score_normalized"].astype(float),
                y=top["label"],
                orientation="h", marker_color=colors,
                text=[f"{float(s):.1f}" for s in top["digenic_score_normalized"]],
                textposition="outside",
                customdata=top[["string_score","shared_pathways","clinvar_pair"]].values,
                hovertemplate="<b>%{y}</b><br>Score: %{x:.1f}/100<br>STRING: %{customdata[0]}<br>Pathways: %{customdata[1]}<br>ClinVar: %{customdata[2]:.2f}<extra></extra>"
            ))
            fig.update_layout(
                plot_bgcolor="#111827", paper_bgcolor="#0a0e1a",
                font_color="#94a3b8", height=max(380, len(top)*38),
                xaxis=dict(range=[0,118], gridcolor="#1e293b", title="Digenic Score (0-100)", color="#64748b"),
                yaxis=dict(autorange="reversed", gridcolor="#1e293b", color="#94a3b8"),
                title=dict(text="Digenic Score Ranking — hover for details", font_color="#e2e8f0", font_size=13),
                margin=dict(l=10, r=60, t=45, b=10),
                hoverlabel=dict(bgcolor="#1a2235", font_color="#e2e8f0")
            )
            st.plotly_chart(fig, use_container_width=True)

        with tab2:
            gene_list = genes_df["gene_symbol"].tolist() if not genes_df.empty else []
            matrix    = pd.DataFrame(0.0, index=gene_list, columns=gene_list)
            for _, row in filtered.iterrows():
                a,b,s = row["gene_a"], row["gene_b"], float(row["digenic_score_normalized"])
                if a in matrix.index and b in matrix.columns:
                    matrix.loc[a,b] = s
                    matrix.loc[b,a] = s
            fig = px.imshow(
                matrix, color_continuous_scale="YlOrRd",
                text_auto=".0f", aspect="auto",
                title="Digenic Score Heatmap — darker = higher score"
            )
            fig.update_layout(
                plot_bgcolor="#111827", paper_bgcolor="#0a0e1a",
                font_color="#94a3b8", height=460,
                title_font_color="#e2e8f0", title_font_size=13
            )
            fig.update_traces(textfont_size=10)
            st.plotly_chart(fig, use_container_width=True)

        with tab3:
            display = filtered[["gene_a","gene_b","digenic_score_normalized",
                               "shared_pathways","string_score","clinvar_pair"]].copy()
            display.columns = ["Gene A","Gene B","Score","Pathways","STRING","ClinVar"]
            display["Score"] = display["Score"].round(1)
            st.dataframe(
                display.style.background_gradient(subset=["Score"], cmap="YlOrRd")
                             .format({"Score":"{:.1f}","ClinVar":"{:.2f}"}),
                use_container_width=True, height=420
            )
    else:
        st.warning("Run layer5_scoring.py first")

# ══════════════════════════════════════════════════
# AI INTERPRETATION
# ══════════════════════════════════════════════════
elif "AI" in page:
    st.markdown('<div class="section-title">AI + Database Interpretation</div>', unsafe_allow_html=True)
    st.markdown("### Clinical Briefs — Gemini AI + MyGene.info + ClinVar")

    if ai_text:
        sections = [s.strip() for s in ai_text.split("="*55) if s.strip() and "RANK" in s]

        if not sections:
            st.info("No interpretations parsed — check layer6_ai_interpretations.txt format")
        else:
            for section in sections:
                lines    = section.strip().split("\n")
                title    = next((l for l in lines if "RANK" in l), "Unknown")
                score_ln = next((l for l in lines if "Digenic Score" in l), "")
                source_ln= next((l for l in lines if "Source" in l), "")
                content  = "\n".join(l for l in lines if l not in [title, score_ln]).strip()

                score_val = 0.0
                try:
                    score_val = float(score_ln.split(":")[1].strip().split("/")[0].split("|")[0].strip())
                except:
                    pass

                source = "Real Data" if "MyGene" in source_ln else "AI" if "Gemini" in source_ln else "Database"
                color  = score_color(score_val)
                badge  = "🤖 Gemini AI" if "Gemini" in source_ln else "🔬 Real Database"

                with st.expander(f"🧬 {title}  |  Score: {score_val:.0f}/100  |  {badge}", expanded=score_val>=90):
                    st.markdown(f"""
                    <div style='background:#111827;border-left:3px solid {color};padding:1rem 1.2rem;
                               border-radius:0 8px 8px 0;font-size:0.85rem;line-height:1.8;
                               color:#cbd5e1;white-space:pre-wrap;font-family:DM Sans,sans-serif;'>
{content}
                    </div>""", unsafe_allow_html=True)
    else:
        st.info("Run layer6_ai.py first to generate interpretations")

# ══════════════════════════════════════════════════
# PIPELINE
# ══════════════════════════════════════════════════
elif "Pipeline" in page:
    st.markdown('<div class="section-title">Analysis Pipeline</div>', unsafe_allow_html=True)
    st.markdown("### 6-Layer Cross-Database Pipeline")

    base = os.path.dirname(os.path.abspath(__file__))
    for num, name, source, outfile, color, desc in [
        ("1",  "Gene Harvesting",      "Open Targets GraphQL",  "outputs/layer1_bbs_genes.csv",          "#00f5c4", "Disease-associated genes filtered by score>0.3"),
        ("2",  "Constraint Filtering", "gnomAD v4",             "outputs/layer2_constrained_genes.csv",  "#6366f1", "Genes with pLI>0.9 — intolerant to mutations"),
        ("3",  "PPI Network",          "STRING DB",             "outputs/layer3_interactions.csv",       "#f59e0b", "Protein interactions with confidence>700"),
        ("4",  "Pathway Analysis",     "Reactome",              "outputs/layer4_pathway_scores.csv",     "#00f5c4", "Shared biological pathway co-membership"),
        ("4b", "Clinical Evidence",    "ClinVar + NCBI",        "outputs/layer4b_clinvar_scores.csv",    "#ef4444", "Confirmed pathogenic variants in patients"),
        ("5",  "Digenic Scoring",      "DiVaS-inspired",        "outputs/layer5_digenic_scores.csv",     "#f59e0b", "Combined score ranking all gene pairs 0-100"),
        ("6",  "AI Interpretation",    "Gemini + MyGene+ClinVar","outputs/layer6_ai_summary.csv",        "#6366f1", "Clinical briefs with validation experiments"),
    ]:
        exists = os.path.exists(os.path.join(base, outfile))
        status = "✅" if exists else "⏳"
        rows   = 0
        if exists:
            try:
                rows = len(pd.read_csv(os.path.join(base, outfile)))
            except:
                rows = 0
        row_txt = f"{rows} rows" if rows > 0 else ""
        st.markdown(f"""
        <div class='card' style='display:flex;align-items:center;gap:14px;padding:0.85rem 1.1rem;margin-bottom:0.4rem;border-left:3px solid {color};'>
            <div style='font-family:Space Mono,monospace;font-size:0.95rem;color:{color};font-weight:700;min-width:34px;'>L{num}</div>
            <div style='flex:1;'>
                <div style='font-size:0.88rem;color:#e2e8f0;font-weight:500;'>{name}</div>
                <div style='font-size:0.7rem;color:#475569;margin-top:2px;'>{source} · {desc}</div>
            </div>
            <div style='font-size:0.7rem;color:#475569;'>{row_txt}</div>
            <div style='font-size:1.1rem;'>{status}</div>
        </div>""", unsafe_allow_html=True)

    st.markdown("<br>", unsafe_allow_html=True)
    col1, col2 = st.columns(2)
    with col1:
        st.markdown('<div class="section-title">Digenic Score Formula</div>', unsafe_allow_html=True)
        st.code("""DS = (string_conf × pLI_A × pLI_B
       × pathway_overlap × clinvar_pair)
     ÷ (AF_A × AF_B × 1e6)

Normalized to 0–100 scale
Inspired by DiVaS algorithm (Gazzo et al. 2016)""", language="text")
    with col2:
        st.markdown('<div class="section-title">Score Interpretation</div>', unsafe_allow_html=True)
        for rng, label, color in [
            ("80–100", "Top priority — strong multi-database evidence",   "#ef4444"),
            ("60–79",  "High priority — validated interaction + pathways","#f59e0b"),
            ("40–59",  "Moderate — good interaction confidence",          "#00f5c4"),
            ("0–39",   "Low — limited evidence",                          "#6366f1"),
        ]:
            st.markdown(f"""
            <div style='display:flex;align-items:center;gap:10px;padding:6px 10px;background:#111827;border-radius:6px;margin-bottom:4px;'>
                <span style='font-family:Space Mono,monospace;font-size:0.8rem;color:{color};min-width:50px;'>{rng}</span>
                <span style='font-size:0.78rem;color:#94a3b8;'>{label}</span>
            </div>""", unsafe_allow_html=True)
