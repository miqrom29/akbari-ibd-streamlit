import streamlit as st
import pandas as pd
import networkx as nx
import plotly.graph_objects as go

st.set_page_config(page_title="Akbari IBD clustering", layout="wide")

st.title("Experimental IBD Clustering for Ancient DNA")

# ───────────────────────── Sidebar: uploads ─────────────────────────

st.sidebar.header("Upload IBD pairs")
ibd_file = st.sidebar.file_uploader(
    "CSV or TSV with id1,id2,total_cM", type=["csv", "tsv"], key="ibd"
)

st.sidebar.header("Optional sample metadata")
meta_file = st.sidebar.file_uploader(
    "Metadata CSV (sample,haplogroup_mt,...)", type=["csv"], key="meta"
)

# ───────────────────────── Helpers ─────────────────────────


def clean_id(x: str) -> str:
    """Convert '[ID](url)' → 'ID', otherwise return the string."""
    s = str(x)
    if s.startswith("[") and "](" in s:
        return s.split("[", 1)[1].split("]", 1)[0]
    return s


# ───────────────────────── Main logic ─────────────────────────

if ibd_file is None:
    st.info("Upload a CSV/TSV file with id1,id2,total_cM to begin.")
else:
    # Detect separator
    sep = "\t" if ibd_file.name.endswith(".tsv") else ","
    df = pd.read_csv(ibd_file, sep=sep)

    # Normalise column names
    cols = {c.lower(): c for c in df.columns}
    id1_col = cols.get("sample1") or cols.get("id1") or list(df.columns)[0]
    id2_col = cols.get("sample2") or cols.get("id2") or list(df.columns)[1]
    tot_col = (
        cols.get("total_cm")
        or cols.get("length_cm")
        or cols.get("tot_cm")
        or list(df.columns)[2]
    )

    # Clean IDs and parse total_cM
    df["sample1"] = df[id1_col].apply(clean_id)
    df["sample2"] = df[id2_col].apply(clean_id)
    df["total_cM"] = pd.to_numeric(df[tot_col], errors="coerce")
    df = df.dropna(subset=["total_cM"]).copy()

    st.subheader("Pairwise IBD (input)")
    st.dataframe(
        df[["sample1", "sample2", "total_cM"]]
        .sort_values("total_cM", ascending=False),
        use_container_width=True,
    )

    # ── Load metadata if provided ──
    meta = None
if meta_file is not None:
    meta = pd.read_csv(meta_file)
    meta_cols = {c.lower(): c for c in meta.columns}

    # identify sample column
    sid_col = meta_cols.get("sample") or list(meta.columns)[0]

    # clean sample IDs ([ID](url) -> ID)
    meta["sample_clean"] = meta[sid_col].apply(clean_id)
    meta = meta.set_index("sample_clean")

    # normalise haplogroup column names
    # support both mt_haplogroup / y_haplogroup and haplogroup_mt / haplogroup_y
    if "haplogroup_mt" not in meta.columns:
        col_mt = meta_cols.get("mt_haplogroup")
        if col_mt:
            meta.rename(columns={col_mt: "haplogroup_mt"}, inplace=True)

    if "haplogroup_y" not in meta.columns:
        col_y = meta_cols.get("y_haplogroup")
        if col_y:
            meta.rename(columns={col_y: "haplogroup_y"}, inplace=True)

    # ── Build graph ──
    G = nx.Graph()
    for _, r in df.iterrows():
        G.add_edge(r["sample1"], r["sample2"], weight=r["total_cM"])

    # Connected components as clusters
    components = list(nx.connected_components(G))
    cluster_map = {}
    for i, comp in enumerate(components, start=1):
        label = f"Cluster {i}"
        for node in comp:
            cluster_map[node] = label

    # ── Sample-level table with optional metadata ──
    rows = []
    for node in G.nodes():
        partners = list(G.neighbors(node))
        totals = [G[node][p]["weight"] for p in partners]
        row = {
            "sample": node,
            "cluster": cluster_map[node],
            "partners": len(partners),
            "max_pairwise_cM": max(totals) if totals else 0.0,
        }
        if meta is not None and node in meta.index:
            for col in meta.columns:
                row[col] = meta.loc[node, col]
        rows.append(row)

    df_samples = pd.DataFrame(rows).sort_values(["cluster", "sample"])
    st.subheader("Sample assignments")
    st.dataframe(df_samples, use_container_width=True)

    # ───────────────────── Cluster selector & favorites ─────────────────────

    if "favorites" not in st.session_state:
        st.session_state["favorites"] = set()

    st.sidebar.header("Cluster view")
    clusters = sorted(df_samples["cluster"].unique())
    selected = st.sidebar.selectbox("Select cluster", ["All"] + clusters)

    col_fav_add, col_fav_clear = st.sidebar.columns(2)
    with col_fav_add:
        if st.button("★ Add to favorites") and selected != "All":
            st.session_state["favorites"].add(selected)
    with col_fav_clear:
        if st.button("Clear"):
            st.session_state["favorites"] = set()

    if st.session_state["favorites"]:
        st.sidebar.markdown("**Favorite clusters:**")
        for c in sorted(st.session_state["favorites"]):
            st.sidebar.write(f"- {c}")

    # Download favorites table
    if st.session_state["favorites"]:
        fav_df = df_samples[df_samples["cluster"].isin(st.session_state["favorites"])]
        st.download_button(
            "Download favorites table",
            fav_df.to_csv(index=False).encode("utf-8"),
            file_name="ibd_favorite_clusters.csv",
            mime="text/csv",
        )

    # ───────────────────── Network graph ─────────────────────

    # Subgraph by cluster
    if selected != "All":
        nodes = [s for s in G.nodes() if cluster_map[s] == selected]
        H = G.subgraph(nodes)
        title = f"IBD network ({selected})"
    else:
        H = G
        title = "IBD network (All clusters – may be large)"

    if H.number_of_nodes() == 0:
        st.warning("No nodes to display for this selection.")
    else:
        # Layout
        pos = nx.spring_layout(H, weight="weight", seed=42)

        node_x, node_y, text_labels, hovertext = [], [], [], []
        for n, (px, py) in pos.items():
            node_x.append(px)
            node_y.append(py)
            text_labels.append(n)
            label = n
            if meta is not None and n in meta.index:
                if "haplogroup_mt" in meta.columns:
                    label += f" | mt: {meta.loc[n, 'haplogroup_mt']}"
                if "haplogroup_y" in meta.columns:
                    label += f" | Y: {meta.loc[n, 'haplogroup_y']}"
            hovertext.append(label)

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
                line=dict(color="rgba(150,150,150,0.4)", width=1),
                hoverinfo="none",
            )
        )
        fig.add_trace(
            go.Scatter(
                x=node_x,
                y=node_y,
                mode="markers+text",
                text=text_labels if selected != "All" else None,
                textposition="top center",
                marker=dict(size=10, color="steelblue"),
                hovertext=hovertext,
                hoverinfo="text",
            )
        )
        fig.update_layout(
            title=title,
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            showlegend=False,
            height=650,
            margin=dict(l=10, r=10, t=50, b=10),
        )

        st.subheader("IBD network")
        st.plotly_chart(fig, use_container_width=True)
