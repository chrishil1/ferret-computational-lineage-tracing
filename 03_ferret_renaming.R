library(dplyr)
library(ggplot2)
library(here)
library(Seurat)
library(scAnalysisR)

source(here("03_ferret_renaming_helpers.R"))

# Useful directories
data_dir <- here("rds_objects")
image_dir <- here("output_images")

# Load in knn relabeled ferret object
ferret_knn_object <- readRDS(file.path(
  data_dir,
  "03_ferret_relabelled_clusters.rds"
))
ferret_knn_object$knn_celltype <- ferret_knn_object$final_combined_celltype

ferret_knn_umap <- plotDimRed(ferret_knn_object,
  col_by = "knn_celltype",
  color = human_palette,
  plot_type = "umap",
  ggrastr = TRUE
)[[1]]

pdf(
  file = file.path(image_dir, "ferret_knn_umap.pdf"),
  width = 9,
  height = 6
)
print(ferret_knn_umap)
dev.off()

# TODO: Find the code that does the majority cell type labeling and put it here

# Load in resolution 7 seurat cluster ferret object
ferret_majority_object <- readRDS(file.path(
  data_dir,
  "03_combined_ferret_res_7.rds"
))

# Visualize the resolution 7 seurat clusters used for majority renaming
ferret_sc_majority <- plotDimRed(ferret_majority_object,
  col_by = "seurat_clusters",
  plot_type = "umap",
  ggrastr = TRUE
)[[1]]

pdf(
  file = file.path(image_dir, "ferret_seurat_clusters_majority.pdf"),
  width = 8,
  height = 6
)
print(ferret_sc_majority)
dev.off()

# Visualize the majority cell type renaming before manual steps
ferret_majority_umap <- plotDimRed(ferret_majority_object,
  col_by = "majority_celltype",
  color = human_palette,
  plot_type = "umap",
  ggrastr = TRUE
)[[1]]

pdf(
  file = file.path(image_dir, "ferret_majority_umap.pdf"),
  width = 8,
  height = 6
)
print(ferret_majority_umap)
dev.off()

# TODO: Find the code that does the reclustering with correct resolution

# Load in ferret object after majority cell type renaming
ferret_final_majority <-
  readRDS(file.path(data_dir, "03_ferret_relabelled_clusters_pre.rds"))

# Visualize the Seurat clusters that gave rise to this manual renaming
ferret_seurat_clusters <- plotDimRed(ferret_final_majority,
  col_by = "seurat_clusters",
  plot_type = "umap", ggrastr = TRUE
)[[1]]

pdf(
  file = file.path(image_dir, "ferret_seurat_clusters.pdf"),
  width = 7,
  height = 6
)
print(ferret_seurat_clusters)
dev.off()

# Relabel cells in cluster 0 manually to be ductal immature
cluster0_cells <- WhichCells(ferret_final_majority, idents = "0")
ferret_final_majority$majority_celltype[cluster0_cells] <- "Ductal Immature"

# Relabel cells in cluster 13 manually to be delta
cluster13_cells <- WhichCells(ferret_final_majority, idents = "13")
ferret_final_majority$majority_celltype[cluster13_cells] <- "Delta"

# Relabel cells in cluster 7 manually to be alpha
cluster7_cells <- WhichCells(ferret_final_majority, idents = "7")
ferret_final_majority$majority_celltype[cluster7_cells] <- "Alpha"

# Relabel cells in cluster 2 manually to be beta
cluster2_cells <- WhichCells(ferret_final_majority, idents = "2")
ferret_final_majority$majority_celltype[cluster2_cells] <- "Beta"

# Visualize the resulting cluster renaming here
ferret_final_majority_umap <- plotDimRed(ferret_final_majority,
  col_by = "majority_celltype",
  color = human_palette,
  plot_type = "umap", ggrastr = TRUE
)[[1]]

pdf(
  file = file.path(image_dir, "ferret_final_majority_umap.pdf"),
  width = 8,
  height = 6
)
print(ferret_final_majority_umap)
dev.off()

# Visualize barplot of the ferret majority cell types across the 4 timepoints
ferret_majority_barplot <- create_barplot(
  "majority_celltype",
  ferret_majority_object,
  "day",
  human_palette
)

pdf(
  file = file.path(image_dir, "ferret_majority_barplot.pdf"),
  width = 6,
  height = 6
)
print(ferret_majority_barplot)
dev.off()
