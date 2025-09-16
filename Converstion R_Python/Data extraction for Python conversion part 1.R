counts_matrix <- GetAssayData(seu, slot = "counts") #raw counts
norm_matrix <- GetAssayData(seu, slot = "data") #norm counts

cell_metadata <- seu@meta.data #meta data

gene_metadata <- data.frame(gene = rownames(counts_matrix))

pca_embeddings <- Embeddings(seurat_obj[["pca"]])
umap_embeddings <- Embeddings(seurat_obj[["umap"]])

# Save count matrix (sparse)
writeMM(counts_matrix, "/Users/anajul/Desktop/Goritz_Lab/Spatial_transcriptomic/Tangram/counts.mtx")
writeMM(norm_matrix, "/Users/anajul/Desktop/Goritz_Lab/Spatial_transcriptomic/Tangram/data.mtx")

# Save metadata
write.csv(cell_metadata, "/Users/anajul/Desktop/Goritz_Lab/Spatial_transcriptomic/Tangram/cell_metadata.csv", quote=FALSE)
write.csv(gene_metadata, "/Users/anajul/Desktop/Goritz_Lab/Spatial_transcriptomic/Tangram/gene_metadata.csv", quote=FALSE)

# Save embeddings
write.csv(pca_embeddings, "/Users/anajul/Desktop/Goritz_Lab/Spatial_transcriptomic/Tangram/pca_embeddings.csv", quote=FALSE)
write.csv(umap_embeddings, "/Users/anajul/Desktop/Goritz_Lab/Spatial_transcriptomic/Tangram/umap_embeddings.csv", quote=FALSE)
