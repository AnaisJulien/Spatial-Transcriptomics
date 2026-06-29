Spatial Transcriptomics Pipeline

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
       ▼
  (cell type annotation, then subset a population of interest)
       │
       ├─────────────────────────────────────────┐
       │ 02_1a — expression-based                │ 02_1b — reference-based
       ▼                                         ▼
┌──────────────────────────────┐   ┌──────────────────────────────────┐
│  02_1a_global_integration    │   │  02_1b_tangram                   │
│       _scvi                  │   │  • Map scRNA-seq reference        │
│  • Re-integrate subset with  │   │    (multiome) onto spatial data  │
│    scVI (e.g. immune cells)  │   │  • Transfer fine-grained cell    │
│  • MPS acceleration support  │   │    type labels to spatial cells  │
│  • Leiden clustering + UMAP  │   │  • Deconvolve mixed spots        │
│  → subset_integrated.h5ad    │   │  → spatial + cell type labels    │
└──────────────────────────────┘   └──────────────────────────────────┘
       │                                         │
       └──────────────┬──────────────────────────┘
                      │  (annotated spatial object)
                      │
       ┌──────────────┴──────────────────────────┐
       │  Downstream analyses                    │
       ▼                                         ▼
┌──────────────────────────┐   ┌──────────────────────────────────────┐
│  Augur                   │   │  COMMOT                              │
│  • Rank cell types by    │   │  • Spatial cell-cell communication   │
│    perturbation score    │   │  • CellChat LR database (mouse)      │
│  • Compare injury types  │   │  • 50 µm distance threshold          │
│                          │   │  • Focused on a specific domain      │
│  • Compare time points   │   │                                      │
│                          │   │  → per-slide commot .h5ad objects    │
│  → AUC scores per cell   │   │  → sender/receiver scores per LR     │
│    type & condition      │   │    pair and pathway                  │
└──────────────────────────┘   └──────────────────────────────────────┘
