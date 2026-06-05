library(dplyr)
library(ggplot2)
library(here)
library(pheatmap)
library(Seurat)
library(scAnalysisR)

source(here("src/common_setup.R"))
source(here("src/01_dataset_visualization_helpers.R"))

# Load in the Seurat object after integratation
merged_object <- readRDS(file.path(data_dir, "02_merged_integrated_object.rds"))

# Subset to just include the ferret cells
ferret_knn_object <- subset(merged_object, subset = batch == "Ferret")

# Process the object to get in an easier to visualize projection
ferret_knn_object <- NormalizeData(ferret_knn_object) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(npcs = 15) %>%
  FindNeighbors(dims = 1:15) %>%
  FindClusters(resolution = 0.8) %>%
  RunUMAP(reduction = "pca", dims = 1:15)

# Visualize the knn cell labels from the human dataset integration
ferret_knn_umap <- plotDimRed(ferret_knn_object,
  col_by = "knn_celltype",
  color = human_palette,
  plot_type = "umap",
  ggrastr = TRUE
)[[1]]

pdf(
  file = file.path(image_dir, "03_ferret_knn_umap.pdf"),
  width = 9,
  height = 6
)
print(ferret_knn_umap)
dev.off()

# Calculate the resolution 7 clustering
ferret_res_7_object <- NormalizeData(ferret_knn_object) %>%
  FindClusters(resolution = 7) %>%
  RunUMAP(reduction = "pca", dims = 1:15)

# For each res 7 cluster, find the majority knn_celltype and use it to rename
cluster_labels <- ferret_res_7_object@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    majority_celltype = names(sort(table(knn_celltype), decreasing = TRUE))[1],
    .groups = "drop"
  )

# Create a mapping that will be easily applied to the Seurat object
cluster_map <- setNames(
  cluster_labels$majority_celltype,
  cluster_labels$seurat_clusters
)

# Assign the new majority celltype names to the object with the mapping
ferret_res_7_object@meta.data$majority_celltype <- cluster_map[
  as.character(ferret_res_7_object$seurat_clusters)
]

# Visualize the resolution 7 seurat clusters used for majority renaming
ferret_sc_majority <- plotDimRed(ferret_res_7_object,
  col_by = "RNA_snn_res.7",
  plot_type = "umap",
  ggrastr = TRUE
)[[1]]

pdf(
  file = file.path(image_dir, "03_ferret_seurat_clusters_majority.pdf"),
  width = 8,
  height = 6
)
print(ferret_sc_majority)
dev.off()

# Visualize the majority cell type renaming before manual steps
ferret_majority_umap <- plotDimRed(ferret_res_7_object,
  col_by = "majority_celltype",
  color = human_palette,
  plot_type = "umap",
  ggrastr = TRUE
)[[1]]

pdf(
  file = file.path(image_dir, "03_ferret_majority_umap.pdf"),
  width = 8,
  height = 6
)
print(ferret_majority_umap)
dev.off()

# Create another celltype naming scheme so we have all steps documented
ferret_res_7_object$final_renamed_celltype <-
  ferret_res_7_object$majority_celltype

# Retrieve the 0.8 resolution plot to perform final renaming
ferret_final_majority <- NormalizeData(ferret_res_7_object) %>%
  FindClusters(resolution = 0.8) %>%
  RunUMAP(reduction = "pca", dims = 1:15)

# Visualize the Seurat clusters that gave rise to this manual renaming
ferret_seurat_clusters <- plotDimRed(ferret_final_majority,
  col_by = "RNA_snn_res.0.8",
  plot_type = "umap", ggrastr = TRUE
)[[1]]

pdf(
  file = file.path(image_dir, "03_ferret_seurat_clusters.pdf"),
  width = 7,
  height = 6
)
print(ferret_seurat_clusters)
dev.off()

# Relabel cells in cluster 0 manually to be ductal immature
cluster0_cells <- WhichCells(ferret_final_majority, idents = "0")
ferret_final_majority$final_renamed_celltype[cluster0_cells] <- "Ductal Immature"

# Relabel cells in cluster 13 manually to be delta
cluster13_cells <- WhichCells(ferret_final_majority, idents = "13")
ferret_final_majority$final_renamed_celltype[cluster13_cells] <- "Delta"

# Relabel cells in cluster 7 manually to be alpha
cluster7_cells <- WhichCells(ferret_final_majority, idents = "7")
ferret_final_majority$final_renamed_celltype[cluster7_cells] <- "Alpha"

# Relabel cells in cluster 2 manually to be beta
cluster2_cells <- WhichCells(ferret_final_majority, idents = "2")
ferret_final_majority$final_renamed_celltype[cluster2_cells] <- "Beta"

# Visualize the resulting cluster renaming here
ferret_final_majority_umap <- plotDimRed(ferret_final_majority,
  col_by = "final_renamed_celltype",
  color = human_palette,
  plot_type = "umap", ggrastr = TRUE
)[[1]]

pdf(
  file = file.path(image_dir, "03_ferret_final_majority_umap.pdf"),
  width = 8,
  height = 6
)
print(ferret_final_majority_umap)
dev.off()

# Visualize barplot of the ferret majority cell types across the 4 timepoints
ferret_majority_barplot <- create_barplot(
  "final_renamed_celltype",
  ferret_final_majority,
  "day",
  human_palette
)

pdf(
  file = file.path(image_dir, "03_ferret_majority_barplot.pdf"),
  width = 6,
  height = 6
)
print(ferret_majority_barplot)
dev.off()

# Get DE heatmap of newly labeled majority cell type
ferret_majority_heatmap <-
  get_de_heatmap(
    seurat_object = ferret_final_majority,
    group_by = "final_renamed_celltype",
    n = 100,
    upregulated = TRUE,
    palette = blue_yellow,
    custom_order = c(
      "T lymphocyte", "Delta", "Beta", "Alpha", "B lymphocyte",
      "Ductal Immature", "Mesenchyme Subtype 1", "Ductal Proliferating",
      "Acinar Proliferating", "Mesenchyme Subtype 2", "Endothelial Lymphatic",
      "Acinar Immature", "Acinar Maturing", "Acinar Progenitor",
      "Pancreatic Progenitor"
    )
  )

pdf(
  file = file.path(image_dir, "03_ferret_majority_heatmap.pdf"),
  width = 15,
  height = 8
)
print(ferret_majority_heatmap)
dev.off()

# Save the updated object to an output Seurat file
saveRDS(
  file = file.path(data_dir, "03_combined_ferret_res_7.rds"),
  ferret_final_majority
)
