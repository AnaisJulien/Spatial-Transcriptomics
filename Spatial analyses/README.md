# Spatial Transcriptomics Pipeline — BANKSY + scVI

A two-stage pipeline for analysing **MERSCOPE** (or similar targeted spatial transcriptomics) data:

1. **Per-section spatial embedding** with [BANKSY](https://github.com/prabhakarlab/Banksy)
2. **Cross-section batch integration** with [scVI](https://scvi-tools.org/)

---

## Overview

```
Raw .h5ad sections
       │
       ▼
┌─────────────────────────────┐
│  01_per_section_banksy      │  ← run once per tissue section
│  • QC: auto-select k_geom   │
│  • Normalize + log1p        │
│  • BANKSY spatial embedding │
│  → *_banksy_full.h5ad       │  1500 features (genes + neighbourhood)
│  → *_banksy_genes.h5ad      │  genes only, all layers
└─────────────────────────────┘
       │
       ├─────────────────────────────────────────┐
       │ (*_banksy_genes.h5ad)                   │ (*_banksy_full.h5ad)
       ▼                                         ▼
┌──────────────────────────────┐   ┌──────────────────────────────────┐
│  02_global_integration_scvi  │   │  04_banksy_aware_integration     │
│  • Concatenate sections      │   │  • Concatenate sections          │
│  • QC filtering              │   │  • Regress out total_counts      │
│  • scVI batch correction     │   │  • PCA on 1500 features          │
│    (raw counts, 500 genes)   │   │  • Harmony batch correction      │
│  • Leiden clustering + UMAP  │   │  • Leiden clustering + UMAP      │
│  → integrated.h5ad           │   │  → banksy_aware_integrated.h5ad  │
└──────────────────────────────┘   └──────────────────────────────────┘
       │
       ▼  (optional)
┌────────────────────────────────────┐
│  02_1_global_integration_scvi      │  sub-population re-integration
│  Variant: retrain scVI on a subset │
│  (e.g. immune cells) with Apple    │
│  Silicon MPS acceleration support  │
└────────────────────────────────────┘
```

The two integration branches serve different purposes:
- **Notebook 02 (scVI)** — integrates on gene expression only. Better for cell type identification.
- **Notebook 04 (Harmony)** — integrates on BANKSY features (genes + spatial neighbourhood). Better for identifying spatially-defined tissue domains.

---

## Notebooks

### `01_per_section_banksy.ipynb`

Processes each tissue section independently.

**Steps:**
1. **k_geom QC** — computes median k-NN distances across candidate values (6, 8, 10, 12) to suggest the best spatial neighbourhood size per section. Override manually in `K_GEOM_PER_SECTION` if needed.
2. **Preprocessing** — normalises to 10 000 counts per cell, log1p-transforms, and stores all layers (`counts`, `normalized`, `transformed`, `scaled`). Raw counts are frozen in `adata.raw` for downstream DE tools.
3. **BANKSY** — builds a 1500-feature matrix combining self-expression (genes) and spatial neighbourhood expression (λ = 0.8, `scaled_gaussian` weights).
4. **Save** — writes two files per section:
   - `*_banksy_full.h5ad` — 1500 features, used for spatial label transfer
   - `*_banksy_genes.h5ad` — genes only, all layers, used as input to scVI

**Key parameters (edit in cell 2):**

| Parameter | Default | Description |
|---|---|---|
| `DEFAULT_K_GEOM` | 8 | Spatial neighbourhood size |
| `lambda_list` | `[0.8]` | BANKSY mixing weight (0 = no neighbourhood) |
| `max_m` | 1 | BANKSY moment order |
| `nbr_weight_decay` | `scaled_gaussian` | Neighbour weight function |

---

### `02_global_integration_scvi.ipynb`

Integrates all sections into a single embedding.

**Steps:**
1. **Load & concatenate** — scans the BANKSY output directory for `*_banksy_genes.h5ad` files, makes `obs_names` unique (`<section>_<barcode>`), sets `adata.X` to raw counts for scVI.
2. **Exclude sparse sections** — drops sections with fewer than `min_cells_section` cells (default 200), which lack enough cells for scVI to learn their batch effect.
3. **QC filtering** — filters on gene count and total count thresholds. No mitochondrial filtering (targeted panel). Plots distributions before filtering.
4. **scVI** — trains on a stratified subsample (~`n_train_cells` cells, evenly drawn per section), then embeds **all** cells. Batch key = `section`.
5. **Cluster & UMAP** — builds a k-NN graph on `X_scVI`, computes UMAP, and runs Leiden at multiple resolutions.
6. **Save** — writes the integrated `.h5ad` with all layers, `X_scVI`, `X_umap_scVI`, and `leiden_*` columns.

**Key parameters (edit in cell 2):**

| Parameter | Default | Description |
|---|---|---|
| `min_genes` / `max_genes` | 17 / 500 | Gene count QC bounds (set max to panel size) |
| `min_counts` / `max_counts` | 5 / 2500 | Total count QC bounds |
| `min_cells_section` | 200 | Minimum cells to keep a section |
| `n_latent` | 20 | scVI latent dimensions |
| `n_layers` | 2 | scVI encoder/decoder depth |
| `max_epochs` | 300 | scVI training cap (early stopping enabled) |
| `n_train_cells` | 150 000 | Cells used for training (stratified by section) |
| `n_neighbors` | 30 | Neighbours for graph / UMAP |
| `leiden_resolutions` | `[0.3, 0.5, 0.8, 1.0, 1.5]` | Clustering resolutions |

> **Why scVI instead of Harmony?** The dataset spans many batches (multiple injury models × multiple timepoints). scVI's generative model handles this better than linear correction methods.

> **Why not use BANKSY features in scVI?** BANKSY neighbourhood features are continuous spatial averages, not count data. scVI models raw counts via a negative binomial likelihood — mixing in BANKSY features would violate that assumption.

---

### `02_1_global_integration_scvi.ipynb`

A variant of notebook 02 for **re-integrating a cell subset** (e.g. immune cells) with refined parameters. Key differences:

- Accepts a pre-filtered `.h5ad` as direct input (skips load/concatenate steps)
- Supports **Apple Silicon MPS acceleration** via a Lightning patch (falls back to CPU automatically)
- `n_train_cells = None` trains on **all** cells (appropriate for smaller subsets)
- Slightly different defaults: `n_latent = 15`, `max_epochs = 500`, `batch_size = 512`

---

## Repository structure

```
.
├── 01_per_section_banksy.ipynb          # Stage 1: per-section BANKSY
├── 02_global_integration_scvi.ipynb     # Stage 2: global scVI integration
├── 02_1_global_integration_scvi.ipynb   # Stage 2 variant: subset re-integration
└── README.md
```

---

## Requirements

### Conda environments

Two separate environments are needed because BANKSY and scVI have conflicting dependencies:

**`banksy` environment** — used for notebook 01 only:
```
scanpy
banksy          # https://github.com/prabhakarlab/Banksy
scipy
scikit-learn
matplotlib
seaborn
```

**`scvi` environment** — used for notebooks 02 and 02_1:
```
scvi-tools
scanpy
pytorch
lightning
```

### Input data format

Each input `.h5ad` file must have:
- `adata.X` — raw count matrix (integers)
- `adata.obs["x"]`, `adata.obs["y"]` — spatial coordinates

---

## Quick start

**1. Edit paths in notebook 01:**
```python
PROJ_DIR = "/path/to/your/input/h5ad/files/"
SAVE_DIR = "/path/to/your/output/directory/"
```

**2. Run Step 1** (k_geom QC) — inspect the distance table, override `K_GEOM_PER_SECTION` for any section if the default doesn't fit.

**3. Run Step 2** — processes all sections sequentially.

**4. Edit paths in notebook 02:**
```python
BANKSY_DIR = "/path/to/banksy/output/"          # = SAVE_DIR from notebook 01
SAVE_PATH  = "/path/to/output/integrated.h5ad"
```

**5. Run notebook 02** step by step (switch kernel before the scVI cell).

---

## Output layers & embeddings

| Key | Location | Description |
|---|---|---|
| `counts` | `adata.layers` | Raw integer counts |
| `normalized` | `adata.layers` | CPM-normalized (target = 10 000) |
| `transformed` | `adata.layers` | log1p of normalized |
| `X_scVI` | `adata.obsm` | scVI latent embedding (20-dim) |
| `X_umap_scVI` | `adata.obsm` | UMAP of scVI latent |
| `leiden_<res>` | `adata.obs` | Leiden cluster labels per resolution |
