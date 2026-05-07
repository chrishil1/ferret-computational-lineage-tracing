library(dplyr)
library(ggplot2)
library(FNN)
library(here)
library(scAnalysisR)
library(scCustomize)
library(Seurat)
library(stringr)

source(here("src/common_setup.R"))
source(here("src/02_integration_helpers.R"))

# Load in Seurat objects
ferret_object <- readRDS(file.path(
  data_dir,
  "01_ferret_processed_object.rds"
))
human_object <- readRDS(here(data_dir, "human-fetal-panc-allcells.rds"))
human_object$Celltype <- droplevels(human_object$Celltype)

# Convert human Seurat object to V5
human_object <- Convert_Assay(human_object,
  assay = c("RNA", "SCT"),
  convert_to = "V5"
)

# Perform integration on human and ferret datasets
merged_seurat <- calculate_human_integration(
  ferret_object = ferret_object,
  human_object = human_object,
  mapping_file = mapping_file
)

# Visualize the merged Seurat object by batch (human/ferret)
human_ferret_batch_umap <- plotDimRed(merged_seurat,
  col_by = "batch",
  color = c("#4285F4", "#EA4335"),
  plot_type = "umap", ggrastr = TRUE, reorder_cells = TRUE
)[[1]]

pdf(
  file = file.path(image_dir, "02_human_ferret_batch_umap.pdf"),
  width = 14,
  height = 12
)
print(human_ferret_batch_umap)
dev.off()

# Calculate k-nearest neighbor cell type mapping from merged seurat object
knn_merged_object <- calculate_human_knn_mapping(
  merged_seurat = merged_seurat,
  mapping_file = mapping_file
)

# Tidy up the cell type naming to look at shared names in UMAP space
knn_merged_object$final_combined_celltype <- gsub(
  "_human$",
  "",
  knn_merged_object$final_combined_celltype
)

# Plot the combined species UMAP with all cells labeled with human cell types
human_ferret_celltype_umap <- plotDimRed(knn_merged_object,
  col_by = "final_combined_celltype",
  color = human_palette,
  plot_type = "umap", ggrastr = TRUE
)[[1]]

pdf(
  file = file.path(image_dir, "02_human_ferret_celltype_umap.pdf"),
  width = 17,
  height = 12
)
print(human_ferret_celltype_umap)
dev.off()

# Subset and plot the ferret only cells in integrated UMAP space
merged_ferret <- subset(knn_merged_object, subset = batch == "Ferret")

ferret_integrated_umap <- plotDimRed(merged_ferret,
  col_by = "final_combined_celltype",
  color = human_palette,
  plot_type = "umap", ggrastr = TRUE
)[[1]]

pdf(
  file = file.path(image_dir, "02_ferret_integrated_umap.pdf"),
  width = 17,
  height = 12
)
print(ferret_integrated_umap)
dev.off()

# Subset and plot the human only cells in integrated UMAP space
merged_human <- subset(knn_merged_object,
  subset = batch == "Human"
)

human_integrated_umap <- plotDimRed(merged_human,
  col_by = "final_combined_celltype",
  color = human_palette,
  plot_type = "umap", ggrastr = TRUE
)[[1]]

pdf(
  file = file.path(image_dir, "02_human_integrated_umap.pdf"),
  width = 17,
  height = 12
)
print(human_integrated_umap)
dev.off()

# Save Seurat objecet to output file
saveRDS(
  file = file.path(data_dir, "02_merged_integrated_object.rds"),
  merged_seurat
)
