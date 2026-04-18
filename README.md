# Akbari IBD Streamlit

Experimental web app for exploring pairwise IBD sharing and clustering in ancient DNA datasets (Akbari et al. 2026, AADR v66, and custom IBD data).
This is a lightweight, fully client‑side prototype built with Streamlit + NetworkX + Plotly.
IBD analysis (core)
- Upload IBD pair files:
  - CSV or TSV with columns like `id1,id2,total_cM`
  - Also accepts headers `sample1,sample2,length_cM`
  - Cleans Markdown links such as `[FRA021.SG](http://FRA021.SG)` → `FRA021.SG`
- Build an IBD network graph:
  - Nodes = samples
  - Edges = IBD sharing (weight = total cM)
  - Infer clusters as connected components in the graph
- Per‑sample summary table:
  - cluster ID
  - number of IBD partners
  - maximum pairwise IBD (cM)
- Interactive 2D graph:
  - Cluster selector (“All” or a specific cluster)
  - Layout based on edge weights

Metadata extraction (NEW)
- Extract sample metadata from:
  - **AADR v66 `.anno` files** (Harvard’s Allen Ancient DNA Resource)
  - **Akbari et al. 2026 supplementary tables** (`.xlsx` format)
  - CSV/TSV files with custom column names
- **Interactive column mapping**: The app asks you to select which columns correspond to:
  - Sample ID (mandatory)
  - mtDNA and Y haplogroups
  - Date (mean BP)
  - Geographic region, country, locality, site
  - ROH segments
  - Family relations (parses `1.0d daughter-mother:AV1:AV2` format)
- Automatically merges multiple metadata files
- Displays haplogroup counts and date distributions per cluster

### ID linking tool (NEW)
- When IBD sample IDs don't match metadata IDs exactly (e.g., `Loschbour` vs `Loschbour.AG`), the app provides a manual linking interface:
  - Upload a CSV mapping file (`ibd_id, metadata_id`)
  - Assign IDs one by one with search & select
  - Mappings persist during the session

## ⚠️ Important: Memory usage and deployment

### Local installation (RECOMMENDED for real analysis)
- The app can process **full AADR v66 metadata** (~23,000 samples, 14 MB `.anno` file) when run **locally**.
- Memory usage: **500 MB – 1.5 GB** depending on IBD file size.
- **Always pre-filter your IBD pairs** before loading:
  - **Use a threshold >40 cM** (e.g., `awk '$3 > 40' ibd.tsv > ibd_filtered.tsv`)
  - This reduces the number of edges and keeps memory under control.
- For large IBD files (>10,000 pairs), the app may become slow; consider downsampling or increasing your local machine’s RAM.

### Streamlit Community Cloud (DEMO ONLY)
- The public demo at `https://akbari-ibd-streamlit.streamlit.app` is **only for testing with very small files**.
- Streamlit Cloud free tier has **200 MB memory limit** – the app will crash with full AADR metadata or large IBD files.
- **Do not upload large `.anno` files (>5 MB) or IBD files with >2,000 pairs** to the cloud version.
- The cloud version is intended for:
  - Quick feature previews
  - Testing with tiny subsets (e.g., first 100 lines of an `.anno` file)
  - Sharing visualisations of small datasets

## 📥 Input format

### Minimal expected input (CSV or TSV):

```csv
id1,id2,total_cM
FRA021.SG,FRA025.SG,1711.4
FRA021.SG,FRA020.SG,1706.04
FRA022.SG,FRA025.SG,3526.57
Alternative headers that also work:

sample1,sample2,total_cM

sample1,sample2,length_cM

IDs can optionally be Markdown links:

text
[FRA021.SG](http://FRA021.SG)    [FRA025.SG](http://FRA025.SG)    1711.4
These are automatically cleaned to FRA021.SG and FRA025.SG.

Metadata input examples
AADR v66 .anno file – first row is the header (long descriptions), second row onward contains data.

Akbari 2026 .xlsx – automatically detects columns, but you can manually remap them.

Any CSV/TSV with sample IDs and additional columns.

🖥️ How to run locally
bash
# create and activate a virtual environment (optional but recommended)
python -m venv .venv
source .venv/bin/activate  # on Windows: .venv\Scripts\activate

# install dependencies
pip install -r requirements.txt

# run the app
streamlit run app.py
Then open the URL shown in the terminal (default: http://localhost:8501).

☁️ Deploying on Streamlit Community Cloud (DEMO only)
Push this repository to GitHub (already done for miqrom29/akbari-ibd-streamlit).

Go to https://share.streamlit.io and sign in with GitHub.

Create a new app:

Repo: miqrom29/akbari-ibd-streamlit

Branch: main

File: app.py

Wait for the build to finish; you’ll get a public URL to share.

⚠️ Remember: The cloud version has only 200 MB RAM. Use it only for small test files.

📚 Notes and limitations
This is an experimental tool for exploratory analysis, not a replacement for ancIBD or KING.

“Clusters” are currently defined as connected components in the IBD graph (single-linkage clustering).

For serious kinship inference, use:

segment‑level IBD with length thresholds (e.g. >12 cM)

methods benchmarked for ancient DNA

Family relations are parsed but not yet converted to expected cM values; the raw relations are displayed in the metadata table.

The app does not perform IBD calling – it only visualises pre‑computed IBD pairs.

🧬 Citation
If you use this tool, please cite:

Akbari, A. et al. (2026) – [[paper details]](https://www.nature.com/articles/s41586-026-10358-1)

The Allen Ancient DNA Resource (AADR) v66 https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/FFIDCW

This repository: https://github.com/miqrom29/akbari-ibd-streamlit

ancient DNA sheets and developing in genarchivist.net
