devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
library(spacexr)
library(SeuratObject)
library(SeuratWrappers, lib.loc = "/Library/Frameworks/R.framework/Versions/4.4-x86_64/Resources/library")
library(SeuratDisk)
library(Seurat)
library(SummarizedExperiment)
library(SpatialExperiment)
library(ggplot2)

#Use Margherita's data as reference
MZ <- readRDS("~/Desktop/multiome/MZ data/multiome_MZ_fragm.rds")
DefaultAssay(MZ) <- "RNA"

Idents(MZ) <- "RNA.weight_res.0.3"

new_cluster_ids <- c(
  '0' = 'Oligodendroglia',
  '1' = 'Macrophages',
  '2' = 'Microglia',
  '3' = 'Astrocytes',
  '4' = 'Neurons_V',
  '5' = 'Ependymal',
  '6' = 'Oligodendroglia',
  '7' = 'Macrophages',
  '8' = 'Neurons_D',
  '9' = 'Oligodendroglia',
  '10' = 'Neurons_D',
  '11' = 'Microglia',
  '12' = 'Neurons_D',
  '13' = 'Ependymal',
  '14' = 'Pericytes',
  '15' = 'Neurons_D',
  '16' = 'Oligodendroglia',
  '17' = 'Endothelial_cells',
  '18' = 'Neurons_D',
  '19' = 'Fibro',
  '20' = 'Fibro'
)

MZ <- RenameIdents(MZ, new_cluster_ids)
MZ[["cluster_RCTD"]] <- Idents(object = MZ)

Idents(MZ) <- "cluster_RCTD"
counts <- GetAssayData(MZ, assay = "RNA", slot = "counts")
cluster <- as.factor(MZ$cluster_RCTD)
names(cluster) <- colnames(MZ)
nUMI <- MZ$nCount_RNA
names(nUMI) <- colnames(MZ)
nUMI <- colSums(counts)
levels(cluster) <- gsub("/", "-", levels(cluster))
reference <- Reference(counts, cluster, nUMI)

#Convert h5ad file from spatial
#Extract coords and counts from python direcly

# Build SpatialExperiment

coords <- read.csv("~/Desktop/Goritz_Lab/Spatial_transcriptomic/Analyse_3/For_RTCD/Slide_12_coords.csv")
rownames(coords) <- coords$cell
coords$cell <- NULL
coords$cell.1 <- NULL
coords <- as.matrix(coords)
mode(coords) <- "numeric"
coords <- as.data.frame(coords)

query.counts <- read.csv("~/Desktop/Goritz_Lab/Spatial_transcriptomic/Analyse_3/For_RTCD/Slide_12_counts.csv")
rownames(query.counts) <- query.counts$cell
query.counts$cell <- NULL
query.counts <- t(query.counts)
new_row <- matrix(100, nrow = 1, ncol = ncol(query.counts))
rownames(new_row) <- "Pseudo_gene"  # optional, give it a name
query.counts <- rbind(query.counts, new_row)

query <- SpatialRNA(coords, query.counts, colSums(query.counts))

RTCD <- create.RCTD(
  spatialRNA = query,
  reference,
  max_cores = 8,
  test_mode = FALSE,
  gene_cutoff = 1e-5,
  fc_cutoff = 0.1,
  gene_cutoff_reg = 1e-5,
  fc_cutoff_reg = 0.25,
  UMI_min = 10,
  UMI_max = 2e7,
  counts_MIN = 2,
  UMI_min_sigma = 100,
  class_df = NULL,
  CELL_MIN_INSTANCE = 10,
  cell_type_names = NULL,
  MAX_MULTI_TYPES = 4,
  keep_reference = FALSE,
  cell_type_profiles = NULL,
  CONFIDENCE_THRESHOLD = 3,
  DOUBLET_THRESHOLD = 25
)

RTCD <- run.RCTD(RTCD, doublet_mode = "doublet")

annotations.df <- RTCD@results$results_df
annotations <- annotations.df$first_type
names(annotations) <- rownames(annotations.df)
write.csv2(annotations.df, "~/Desktop/Goritz_Lab/Spatial_transcriptomic/Analyse_3/For_RTCD/annotations_df.csv")

#How to plot results
resultsdir <- '~/Desktop/Goritz_Lab/Spatial_transcriptomic/Analyse_3/For_RTCD' ## you may change this to a more accessible directory on your computer.
dir.create(resultsdir)

results <- RTCD@results
norm_weights = normalize_weights(results$weights) 
cell_type_names <- RTCD@cell_type_info$info[[2]] #list of cell type names
spatialRNA <- RTCD@spatialRNA

# Plots the confident weights for each cell type as in full_mode (saved as 'results/cell_type_weights_unthreshold.pdf')
plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights)

# Plot proportions of each cell type

barplot(prop.table(table(annotations.df$spot_class)),
        main = "Proportion of Spot Classes",
        ylab = "Proportion",
        xlab = "Spot Class", ylim = c(0, 1),
        col = "steelblue")

# Plot the combination of cell type that are in doublets
par(mar = c(3, 12, 4, 2))
# Subset to doublet_certain
doublets <- subset(annotations.df, spot_class == "doublet_certain")

# Build combinations from first_type + second_type
doublets$cell_combinations <- paste(doublets$first_type, doublets$second_type, sep = "+")

# Proportions
prop_doublets <- prop.table(table(doublets$cell_combinations))

# Sort by proportion
prop_doublets <- sort(prop_doublets)

# Horizontal barplot
bp <- barplot(
  prop_doublets,
  horiz = TRUE,
  las = 1,   # horizontal y-axis labels
  col = "steelblue",
  ylab = "Doublet Composition",
  main = "Proportion of Doublet Compositions",
  xlim = c(0, 0.1)   # scale proportion axis 0–1
)

# Extract the weights to use them as metadata in python 

W <- RTCD@results$weights
weights <- as.data.frame(as.matrix(W))
colnames(weights) <- paste0(colnames(weights), "_rtcd_weights")
write.csv2(weights, "~/Desktop/Goritz_Lab/Spatial_transcriptomic/Analyse_3/For_RTCD/weights.csv")



