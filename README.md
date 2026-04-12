# Akbari IBD Streamlit

Experimental web app for exploring **pairwise IBD sharing and clustering** in the Akbari et al. ancient DNA dataset.

This is a lightweight, fully client‑side prototype built with **Streamlit + NetworkX + Plotly**.

---

## Features

- Upload IBD pair files:
  - CSV or TSV with columns like `id1,id2,total_cM`  
  - Also accepts headers `sample1,sample2,length_cM`
  - Cleans Markdown links such as `[FRA021.SG](http://FRA021.SG)` → `FRA021.SG`
- Build an IBD **network graph**:
  - Nodes = samples
  - Edges = IBD sharing (weight = total cM)
- Infer **clusters** as connected components in the graph
- Per‑sample summary table:
  - cluster ID
  - number of IBD partners
  - maximum pairwise IBD (cM)
- Interactive 2D graph:
  - Cluster selector (“All” or a specific cluster)
  - Layout based on edge weights

---

## Input format

Minimal expected input (CSV or TSV):

```text
id1,id2,total_cM
FRA021.SG,FRA025.SG,1711.4
FRA021.SG,FRA020.SG,1706.04
FRA022.SG,FRA025.SG,3526.57
...
```

Alternative headers that also work:

- `sample1,sample2,total_cM`
- `sample1,sample2,length_cM`

IDs can optionally be Markdown links:

```text
[FRA021.SG](http://FRA021.SG)    [FRA025.SG](http://FRA025.SG)    1711.4
```

These are automatically cleaned to `FRA021.SG` and `FRA025.SG`.

---

## How to run locally

```bash
# create and activate a virtual environment (optional but recommended)
python -m venv .venv
source .venv/bin/activate  # on Windows: .venv\Scripts\activate

# install dependencies
pip install -r requirements.txt

# run the app
streamlit run app.py
```

Then open the URL shown in the terminal (default: `http://localhost:8501`).

---

## Deploying on Streamlit Community Cloud

1. Push this repository to GitHub (already done for `miqrom29/akbari-ibd-streamlit`).
2. Go to https://share.streamlit.io and sign in with GitHub.
3. Create a new app:
   - Repo: `miqrom29/akbari-ibd-streamlit`
   - Branch: `main`
   - File: `app.py`
4. Wait for the build to finish; you’ll get a public URL to share.

---

## Notes and limitations

- This is an **experimental** tool for exploratory analysis, not a replacement for ancIBD or KING.
- “Clusters” are currently defined as **connected components** in the IBD graph.  
  For serious kinship inference, use:
  - segment‑level IBD with length thresholds (e.g. > 12 cM)
  - methods benchmarked for ancient DNA.

---

## License

You are free to reuse and modify this code for research and educational purposes.
