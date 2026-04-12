import streamlit as st
import pandas as pd
import networkx as nx
import plotly.graph_objects as go

st.title("Experimental IBD Clustering for Ancient DNA")

st.sidebar.header("Upload IBD pairs")
file = st.sidebar.file_uploader("CSV or TSV with id1,id2,total_cM", type=["csv", "tsv"])

if file is not None:
    sep = "\t" if file.name.endswith(".tsv") else ","
    df = pd.read_csv(file, sep=sep)

    cols = {c.lower(): c for c in df.columns}
    id1 = cols.get("sample1") or cols.get("id1") or list(df.columns)[0]
    id2 = cols.get("sample2") or cols.get("id2") or list(df.columns)[1]
    tot = cols.get("total_cm") or cols.get("length_cm") or list(df.columns)[2]

    def clean_id(x):
        s = str(x)
        if s.startswith("[") and "](" in s:
            return s.split("[", 1)[1].split("]", 1)[0]
        return s

    df["sample1"] = df[id1].apply(clean_id)
    df["sample2"] = df[id2].apply(clean_id)
    df["total_cM"] = pd.to_numeric(df[tot], errors="coerce")
    df = df.dropna(subset=["total_cM"])

    st.subheader("Pairwise IBD (input)")
    st.dataframe(
        df[["sample1", "sample2", "total_cM"]]
        .sort_values("total_cM", ascending=False)
    )

    G = nx.Graph()
    for _, r in df.iterrows():
        G.add_edge(r["sample1"], r["sample2"], weight=r["total_cM"])

    components = list(nx.connected_components(G))
    cluster_map = {}
    for i, comp in enumerate(components, start=1):
        for node in comp:
            cluster_map[node] = f"Cluster {i}"

    st.subheader("Sample assignments")
    rows = []
    for node in G.nodes():
        partners = list(G.neighbors(node))
        totals = [G[node][p]["weight"] for p in partners]
        rows.append(
            {
                "sample": node,
                "cluster": cluster_map[node],
                "partners": len(partners),
                "max_pairwise_cM": max(totals) if totals else 0.0,
            }
        )
    df_samples = pd.DataFrame(rows).sort_values(["cluster", "sample"])
    st.dataframe(df_samples)

    st.sidebar.header("Cluster view")
    clusters = sorted(df_samples["cluster"].unique())
    selected = st.sidebar.selectbox("Select cluster", ["All"] + clusters)

    if selected != "All":
        nodes = [s for s in G.nodes() if cluster_map[s] == selected]
        H = G.subgraph(nodes)
    else:
        H = G

    pos = nx.spring_layout(H, weight="weight", seed=42)
    x, y, text = [], [], []
    for n, (px, py) in pos.items():
        x.append(px)
        y.append(py)
        text.append(n)

    edge_x, edge_y = [], []
    for u, v in H.edges():
        edge_x += [pos[u][0], pos[v][0], None]
        edge_y += [pos[u][1], pos[v][1], None]

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=edge_x,
            y=edge_y,
            mode="lines",
            line=dict(color="rgba(150,150,150,0.5)", width=1),
            hoverinfo="none",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=x,
            y=y,
            mode="markers+text",
            text=text,
            textposition="top center",
            marker=dict(size=10, color="steelblue"),
            hovertext=text,
            hoverinfo="text",
        )
    )
    fig.update_layout(
        title=f"IBD network ({selected})",
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        showlegend=False,
        height=600,
    )
    st.subheader("IBD network")
    st.plotly_chart(fig, use_container_width=True)
else:
    st.info("Upload a CSV/TSV file with id1,id2,total_cM to begin.")
