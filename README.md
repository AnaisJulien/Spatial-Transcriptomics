[README_repo.md](https://github.com/user-attachments/files/29557747/README_repo.md)
# Spatial Transcriptomics Pipeline

Analysis pipeline for MERSCOPE spatial transcriptomics data, covering spatial embedding, batch integration, cell type annotation, and downstream analyses.
---

## Repository structure

## Folders

### `Global-analyses/`

Core pipeline — runs on every dataset.

| Notebook | Description |
|---|---|
| `Manual_slide_cluster.ipynb` | Optional Stage 0. Interactive Napari-based tool to draw ROIs on a spatial section and split it into separate `.h5ad` files. Use when one slide contains multiple distinct tissue sections. |
| `01_per_section_banksy.ipynb` | Preprocesses each section (normalize → log1p) and runs BANKSY spatial embedding (λ=0.8). Saves two files per section: `*_banksy_genes.h5ad` (genes only) and `*_banksy_full.h5ad` (1500 features including neighbourhood). |
| `02_global_integration_scvi.ipynb` | Concatenates all `*_banksy_genes.h5ad` files, runs QC, trains scVI for batch correction, and clusters with Leiden at multiple resolutions. |
| `02_1_global_integration_scvi.ipynb` | Re-integrates a cell subset (e.g. immune cells) with refined scVI parameters. Supports Apple Silicon MPS acceleration. |
| `04_banksy_aware_integration.ipynb` | Integrates all `*_banksy_full.h5ad` files using PCA + Harmony on the full 1500-feature BANKSY matrix. Better for identifying spatially-defined tissue domains. |

---

### `Spatial-neighbourhood/`

Identifies tissue domains that recur across conditions.

| Notebook | Description |
|---|---|
| `Banksy_domain_analysis.ipynb` | Runs BANKSY + Harmony per group (one injury type × timepoint), extracts domain expression profiles, then clusters domains cross-condition into **meta-domains** via hierarchical correlation clustering. |

---

### `Add-on-analyses/`

Cell type annotation and downstream analyses on the integrated object.

| Notebook | Description |
|---|---|
| `Tangram.ipynb` | Maps a multiome scRNA-seq reference (immune cells + fibroblasts) onto spatial cells using Tangram to transfer fine-grained cell type labels. |
| `Augur.ipynb` | Ranks immune cell types by transcriptional perturbation (AUC) across injury models and phases using Augur (via `pertpy`). |
| `Commot.ipynb` | Infers spatial cell-cell communication using COMMOT and the CellChat ligand-receptor database (mouse, 50 µm threshold). |

---

### `Differential-gene-expression/`

Pseudobulk DEG analysis between injury models.

| Notebook | Description |
|---|---|
| `DEG.ipynb` | Pseudobulk differential expression between two injury conditions, stratified by phase. Three methods provided: Wilcoxon (quick exploration), edgeR (recommended), DESeq2 (Python-native alternative). Includes automatic bias correction for transcript detection efficiency differences across slides. |

---

### `Converstion_R_Python/`

Utilities for moving data between R and Python.

| File | Description |
|---|---|
| `Data_extraction_for_Python_conve…` | Exports count matrices and metadata from a Seurat object to files readable by Python/AnnData. |
| `Seurat_object_conversion_part2…` | Converts a Seurat object to AnnData format. |
| `Tutorial` | Step-by-step guide for the conversion workflow. |

---

## Requirements

Two conda environments cover the full pipeline:

| Environment | Used for |
|---|---|
| `banksy` | `01_per_section_banksy`, `Manual_slide_cluster` |
| `scvi` | All other notebooks |

**`banksy`:** `scanpy · banksy · scipy · scikit-learn · napari · shapely · magicgui`

**`scvi`:** `scvi-tools · scanpy · pytorch · lightning · pertpy · commot · tangram-sc · harmonypy · pydeseq2 (optional)`

For DEG with edgeR: R installation with the `edgeR` Bioconductor package.

## Input format

Each input `.h5ad` must have:
- `adata.X` — raw count matrix (integers)
- `adata.obs["x"]`, `adata.obs["y"]` — spatial coordinates
- `adata.obs["batch"]` or `adata.obs["sample"]` — slide/section ID
