# Spatial Transcriptomics Pipeline

## Pipeline overview

```
┌─────────────────────────────────────────────────────┐
│               Raw .h5ad sections                    │
└─────────────────────────┬───────────────────────────┘
                          │
            ┌─────────────┴─────────────┐
            │ (optional, manual)        │
            ▼                           │
┌────────────────────────────┐          │
│ Manual_section_cluster      │          │
│ · Napari ROI annotation     │          │
│ · Split slide into sections │          │
└─────────────┬────────────────┘         │
              │                          │
              └─────────────┬────────────┘
                            ▼
          ┌───────────────────────────────┐
          │      01_per_section_banksy    │
          │  · QC · auto-select k_geom   │
          │  · Normalize + log1p         │
          │  · BANKSY spatial embedding  │
          └──────────────┬────────────────┘
                         │
           ┌─────────────┴──────────────┐
           │                            │
    *_banksy_genes               *_banksy_full
           │                            │
           ▼                            ▼
┌──────────────────────┐   ┌────────────────────────────┐
│ 02_global_integration│   │ 04_banksy_aware_integration│
│        _scvi         │   │ · Regress out total_counts │
│ · Concatenate        │   │ · PCA on 1500 features     │
│ · QC filter          │   │ · Harmony batch correction │
│ · scVI batch correct │   │ · Leiden + UMAP            │
│ · Leiden + UMAP      │   └────────────────────────────┘
└──────────┬───────────┘
           │
           ▼
┌──────────────────────────────────────┐
│  Cell type annotation                │
│  → subset population of interest     │
└──────────┬───────────────────────────┘
           │
           │ expression-based          
           ▼
     ┌──────────┐    
     │  02_1    │   
     │  scVI    │   
     │  subset  │    
     │  re-integ│    
     └────┬─────┘   
          │
          │
          ▼
┌─────────────────────────┐
│  Annotated spatial      │
│      object             │
└──────────┬──────────────┘
```

## Notebooks

### `Manual_section_cluster.ipynb`
Optional preprocessing step, used only when a slide needs manual intervention before entering the main pipeline.

Interactive [Napari](https://napari.org/)-based tool to draw regions of interest (ROIs) directly on the spatial point cloud of a section, then split the slide into separate `.h5ad` files per ROI. Use this when automated clustering (BANKSY, Leiden) doesn't cleanly separate regions you can recognise by eye — for example splitting a single slide that contains multiple physically distinct tissue sections, or manually delimiting a lesion/scar boundary.

Two annotation modes:
1. **Without gene expression** — draw polygons directly on the spatial coordinates.
2. **With gene expression** — highlight candidate cells based on a marker gene signature (e.g. fibrotic or immune markers) before drawing, making it easier to trace a region defined by expression rather than shape alone.

**Input:** merged spatial `.h5ad` with `obs["sample"]`, `obs["x"]`, `obs["y"]`
**Output:** one `.h5ad` per manually annotated ROI, ready to feed into `01_per_section_banksy.ipynb`

> **Note:** requires a graphical environment — Napari opens a desktop window and will not run headless or on a remote/server-only Jupyter session.


Processes each tissue section independently.

**Steps:**
1. **k_geom QC** — computes median k-NN distances to suggest the best spatial neighbourhood size per section
2. **Preprocessing** — normalises to 10 000 counts per cell, log1p-transforms, stores all layers (`counts`, `normalized`, `transformed`, `scaled`)
3. **BANKSY** — builds a 1500-feature matrix combining self-expression and spatial neighbourhood expression (λ = 0.8)
4. **Save** — writes two files per section:
   - `*_banksy_full.h5ad` — 1500 features, used for spatial-aware integration (notebook 04)
   - `*_banksy_genes.h5ad` — genes only, all layers, used as input to scVI (notebook 02)

**Key parameters:**

| Parameter | Default | Description |
|---|---|---|
| `DEFAULT_K_GEOM` | 8 | Spatial neighbourhood size |
| `lambda_list` | `[0.8]` | BANKSY mixing weight |
| `max_m` | 1 | BANKSY moment order |
| `nbr_weight_decay` | `scaled_gaussian` | Neighbour weight function |

---

### `02_global_integration_scvi.ipynb`
Integrates all sections on gene expression using scVI.

**Steps:**
1. Load & concatenate all `*_banksy_genes.h5ad` files
2. Exclude sparse sections (< `min_cells_section` cells)
3. QC filtering on gene count and total count thresholds
4. scVI training on a stratified subsample, then embed all cells
5. Leiden clustering at multiple resolutions + UMAP
6. Save integrated `.h5ad`

**Key parameters:**

| Parameter | Default | Description |
|---|---|---|
| `min_genes` / `max_genes` | 17 / 500 | Gene count QC bounds |
| `min_counts` / `max_counts` | 5 / 2500 | Total count QC bounds |
| `min_cells_section` | 200 | Minimum cells to keep a section |
| `n_latent` | 20 | scVI latent dimensions |
| `n_layers` | 2 | scVI encoder/decoder depth |
| `max_epochs` | 300 | Training cap (early stopping enabled) |
| `n_train_cells` | 150 000 | Cells used for training |
| `leiden_resolutions` | `[0.3, 0.5, 0.8, 1.0, 1.5]` | Clustering resolutions |

> **Why scVI over Harmony?** The dataset spans many batches (multiple injury models × multiple timepoints). scVI's generative model handles this better than linear correction.

> **Why not use BANKSY features in scVI?** BANKSY neighbourhood features are continuous spatial averages, not counts. scVI models raw counts via a negative binomial likelihood — feeding BANKSY features would violate that assumption.

---

### `02_1a_global_integration_scvi.ipynb`
Re-integrates a cell subset (e.g. immune cells) with refined scVI parameters.

**Differences from 02:**

| | `02` | `02_1a` |
|---|---|---|
| Input | Scans BANKSY dir, loads all sections | Pre-filtered `.h5ad` directly |
| Training cells | Stratified subsample ~150k | All cells (`n_train_cells = None`) |
| `n_latent` | 20 | 15 |
| `max_epochs` | 300 | 500 |
| Hardware | CPU only | CPU + Apple Silicon MPS (auto-fallback) |

---

### `02_1b_tangram.ipynb`
Reference-based label transfer using [Tangram](https://github.com/broadinstitute/Tangram) (via `scvi.external.Tangram`).

Maps a multiome scRNA-seq reference (immune cells + fibroblasts) onto spatial cells to assign fine-grained cell type identities directly from a well-annotated reference atlas. Run separately for each cell population, then labels are merged back into the spatial object.

**Input:** integrated spatial `.h5ad` + scRNA-seq multiome reference `.h5ad`
**Output:** spatial object with transferred cell type labels

---

### `04_banksy_aware_integration.ipynb`
Spatially-aware integration using the full 1500-feature BANKSY objects.

Unlike notebook 02 (genes only, scVI), this integrates cells using both their own expression and their spatial neighbourhood signal. Harmony is used instead of scVI because BANKSY features are continuous, not counts.

**Steps:**
1. Load & concatenate all `*_banksy_full.h5ad` files
2. Exclude sparse sections
3. QC
4. Regress out `total_counts` across all 1500 features
5. PCA (50 components) + Harmony batch correction
6. Leiden clustering + UMAP

**Key parameters:**

| Parameter | Default | Description |
|---|---|---|
| `n_pcs` | 50 | PCA components (more than usual — 1500 features have more variance) |
| `max_iter_harmony` | 50 | Harmony iterations |
| `leiden_resolutions` | `[3.5]` | Clustering resolution |
| `min_cells_section` | 100 | Minimum cells to keep a section |

---

## Repository structure

```
.
├── Manual_section_cluster.ipynb          # Stage 0 (optional): manual ROI annotation
├── 01_per_section_banksy.ipynb           # Stage 1: per-section BANKSY
├── 02_global_integration_scvi.ipynb      # Stage 2: global scVI integration (all cells)
├── 02_1a_global_integration_scvi.ipynb   # Stage 2a: subset scVI re-integration
├── 02_1b_tangram.ipynb                   # Stage 2b: reference-based label transfer
├── 04_banksy_aware_integration.ipynb     # Stage 2 parallel: spatially-aware Harmony
└── README.md
```

---

## Requirements

### Conda environments

| Environment | Used for |
|---|---|
| `banksy` | Notebook 01 |
| `scvi` | Notebooks 02, 02_1a, 02_1b, 04, Augur, COMMOT |

**`Manual_section_cluster.ipynb`** has its own lightweight requirements and can run in either environment, or a separate one: `napari`, `shapely`, `magicgui` (in addition to `scanpy`, `numpy`, `matplotlib`). It requires a graphical display.

**`banksy`:**
```
scanpy · banksy · scipy · scikit-learn · matplotlib · seaborn
```

**`scvi`:**
```
scvi-tools · scanpy · pytorch · lightning · pertpy · commot · tangram-sc
```

### Input format
Each input `.h5ad` must have:
- `adata.X` — raw count matrix (integers)
- `adata.obs["x"]`, `adata.obs["y"]` — spatial coordinates

---

## Output layers & embeddings

| Key | Location | Description |
|---|---|---|
| `counts` | `adata.layers` | Raw integer counts |
| `normalized` | `adata.layers` | CPM-normalized (target = 10 000) |
| `transformed` | `adata.layers` | log1p of normalized |
| `X_scVI` | `adata.obsm` | scVI latent embedding |
| `X_umap_scVI` | `adata.obsm` | UMAP of scVI latent |
| `X_pca_harmony` | `adata.obsm` | Harmony-corrected PCA (notebook 04) |
| `X_umap_banksy` | `adata.obsm` | UMAP of Harmony embedding (notebook 04) |
| `leiden_<res>` | `adata.obs` | Leiden cluster labels per resolution |
