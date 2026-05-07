library(dplyr)
library(ggplot2)
library(grid)
library(here)
library(pheatmap)
library(scAnalysisR)
library(Seurat)

source(here("src/common_setup.R"))
source(here("src/01_dataset_visualization_helpers.R"))

########################
# Ferret dataset plots #
########################
# ferret_object <- readRDS(file.path(
#   data_dir,
#   "01_ferret_processed_object.rds"
# ))
# 
# # TODO - figure out most effective way to plot this
# 
# new_ferret_object <- NormalizeData(new_ferret_object)
# new_ferret_object <- FindVariableFeatures(new_ferret_object)
# new_ferret_object <- ScaleData(new_ferret_object)
# new_ferret_object <- RunPCA(new_ferret_object, npcs = 15)
# new_ferret_object <- FindNeighbors(new_ferret_object, dims = 1:15)
# new_ferret_object <- FindClusters(new_ferret_object, resolution = 0.8)
# new_ferret_object <- RunUMAP(new_ferret_object, reduction = "pca", dims = 1:15)

# Load in Seurat object
ferret_object <- readRDS(file.path(
  data_dir,
  "01_initial_reclustered_ferret_object.rds"
))
ferret_object$day <- gsub("_b2$", "", ferret_object$sample)

# Visualize ferret day distribution in UMAP space
ferret_day_umap <- plotDimRed(ferret_object,
  col_by = "day",
  color = c("#4285F4", "#34A853", "#FBBC05", "#EA4335"),
  plot_type = "umap",
  ggrastr = TRUE
)[[1]]

pdf(
  file = file.path(image_dir, "01_ferret_day_umap.pdf"),
  width = 7,
  height = 6
)
print(ferret_day_umap)
dev.off()

# Visualize ferret initial celltype distribution in UMAP space
ferret_celltype_umap <- plotDimRed(ferret_object,
  col_by = "final_combined_celltype",
  color = ferret_palette,
  plot_type = "rna.umap",
  ggrastr = TRUE
)[[1]]

pdf(
  file = file.path(image_dir, "01_ferret_celltype_umap.pdf"),
  width = 7.5,
  height = 6
)
print(ferret_celltype_umap)
dev.off()

# Switch the order of the day variable so it plots correctly on a heatmap
ferret_object$fixed_day <- sub("^([0-9]+)day$", "day\\1", ferret_object$day)

# Visualize ferret DE genes by day in a heatmap
ferret_day_de_heatmap <- get_de_heatmap(
  seurat_object = ferret_object,
  group_by = "fixed_day",
  n = 100,
  custom_order = c("day22", "day24", "day28", "day35"),
  upregulated = TRUE,
  palette = plasma
)

pdf(
  file = file.path(image_dir, "01_ferret_day_de_heatmap.pdf"),
  width = 10,
  height = 8
)
print(ferret_day_de_heatmap)
dev.off()

# Plot ferret celltype distribution by day in barplot
ferret_day_barplot <- create_barplot(
  "final_combined_celltype",
  ferret_object,
  "day",
  ferret_palette
)

pdf(
  file = file.path(image_dir, "01_ferret_day_barplot.pdf"),
  width = 6,
  height = 6
)
print(ferret_day_barplot)
dev.off()

#######################
# Human dataset plots #
#######################
# Load in Seurat object
human_object <- readRDS(here(data_dir, "human-fetal-panc-allcells.rds"))
human_object$Celltype <- factor(human_object$Celltype,
  levels = sort(levels(human_object$Celltype))
)

DefaultAssay(human_object) <- "RNA"
human_object <- NormalizeData(human_object)

# Visualize human PCW distribution in UMAP space
human_pcw_umap <- plotDimRed(human_object,
  col_by = "PCW",
  color = c("#4285F4", "#FBBC05", "#EA4335"),
  plot_type = "umap",
  ggrastr = TRUE
)[[1]]

pdf(
  file = file.path(image_dir, "01_human_pcw_umap.pdf"),
  width = 15,
  height = 12
)
print(human_pcw_umap)
dev.off()

# Visualize human celltype distribution in UMAP space
human_celltype_umap <- plotDimRed(human_object,
  col_by = "Celltype",
  plot_type = "umap",
  color = human_palette,
  ggrastr = TRUE
)[[1]]

pdf(
  file = file.path(image_dir, "01_human_celltype_umap.pdf"),
  width = 17,
  height = 12
)
print(human_celltype_umap)
dev.off()

# Visualize human DE genes by celltype in a heatmap
human_celltype_de_heatmap <- get_de_heatmap(
  seurat_object = human_object,
  group_by = "Celltype",
  n = 100,
  upregulated = TRUE,
  custom_order = c(
    "Acinar Maturing", "Acinar Secretory", "Acinar Proliferating",
    "Acinar Immature", "Acinar Progenitor", "Pancreatic Progenitor",
    "Dendritic Cell", "Beta", "B lymphocyte", "T lymphocyte", "Mono-Macrophage",
    "Mesenchyme Subtype 3", "Mesenchyme Subtype 1", "Mesenchyme Subtype 4",
    "Mesenchyme Subtype 2", "Endothelial", "Mesenchyme Proliferating",
    "Schwann", "Delta", "Ductal Proliferating", "Ductal ECM", "Ductal Immature",
    "Ductal Maturing", "Ducto-Acinar", "Stromal", "Erythroblast",
    "Phagocytic Macrophage", "Endothelial Lymphatic", "Alpha", "Neuron"
  ),
  palette = blue_yellow
)

pdf(
  file = file.path(image_dir, "01_human_celltype_de_heatmap.pdf"),
  width = 15,
  height = 10
)
print(human_celltype_de_heatmap)
dev.off()

# Plot human celltype distribution by PCW in barplot
human_pcw_barplot <- create_barplot(
  "Celltype",
  human_object,
  "PCW",
  human_palette
)

pdf(
  file = file.path(image_dir, "01_human_pcw_barplot.pdf"),
  width = 6,
  height = 6
)
print(human_pcw_barplot)
dev.off()
