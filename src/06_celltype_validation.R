library(dplyr)
library(here)
library(JASPAR2020)
library(pheatmap)
library(Seurat)
library(TFBSTools)

source(here("src/common_setup.R"))
source(here("src/04_celltype_validation_helpers.R"))
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

# Define the cell types that are found in both datasets
shared_celltypes <- c(
  "Alpha", "Beta", "Delta", "Pancreatic Progenitor", "Acinar Immature",
  "Acinar Maturing", "Acinar Progenitor", "Acinar Proliferating",
  "Ductal Immature", "Ductal Proliferating"
)

# Load in Seurat objects for the two species
ferret_object <- readRDS(file.path(data_dir, "05_final_ferret.rds"))
human_object <- readRDS(file.path(data_dir, "05_final_human.rds"))
human_subset <- subset(human_object, subset = Celltype %in% shared_celltypes)

ferret_genes <- FindAllMarkers(ferret_object, group.by = "majority_celltype")
human_genes <- FindAllMarkers(human_object, group.by = "Celltype")

# Function to convert cell types into suffixes for variable naming ease
make_varname <- function(celltype) gsub(" ", "_", tolower(celltype))

# Loop over all cell types
for (celltype in shared_celltypes) {
  suffix <- make_varname(celltype)

  # Subset Seurat objects based on the current cell type of interest
  assign(paste0("ferret_", suffix), subset(
    ferret_object,
    majority_celltype == celltype
  ))
  assign(paste0("human_", suffix), subset(
    human_object,
    Celltype == celltype
  ))

  # Filter & rank ferret markers for this cell type
  ferret_celltype_genes <- subset(ferret_genes, cluster == celltype) %>%
    filter(abs(avg_log2FC) > 0.5, p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    arrange(cluster, desc(avg_log2FC)) %>%
    slice_head(n = 2000)

  # Filter & rank human markers for this cell type
  human_celltype_genes <- subset(human_genes, cluster == celltype) %>%
    filter(abs(avg_log2FC) > 0.5, p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    arrange(cluster, desc(avg_log2FC)) %>%
    slice_head(n = 2000)

  # Intersect, then keep only the top 100 genes in both datasets
  intersected <- intersect(
    ferret_celltype_genes$gene,
    human_celltype_genes$gene
  )
  top100 <- head(intersected, 100)

  assign(paste0("intersected_genes_", suffix), top100)
}

combined_gene_list <- lapply(shared_celltypes, function(celltype) {
  get(paste0("intersected_genes_", make_varname(celltype)))
})
names(combined_gene_list) <- make_varname(shared_celltypes)

combined_genes_unique <- unname(unlist(combined_gene_list))

# Calculate the initial ferret heatmap including all genes extracted
initial_ferret_heatmap <- get_validation_heatmap(
  ferret_object,
  custom_order = shared_celltypes,
  idents = "majority_celltype",
  genes_of_interest = combined_genes_unique,
  cluster_rows = FALSE,
  show_rownames = FALSE
)

# Calculate the initial human heatmap including all genes extracted
initial_human_heatmap <- get_validation_heatmap(
  human_subset,
  custom_order = shared_celltypes,
  idents = "Celltype",
  genes_of_interest = combined_genes_unique,
  cluster_rows = FALSE,
  show_rownames = FALSE
)

# Get expression data from the first heatmap expression
expr_ferret <- initial_ferret_heatmap$expression_data
expr_human <- initial_human_heatmap$expression_data

genes_per_cluster <- list()
all_genes_ordered <- c()

# For each celltype, only extract the genes that had high expression z-scores
# (> 2) in both species datasets
for (cluster in shared_celltypes) {
  if (cluster %in% colnames(expr_ferret) && cluster %in% colnames(expr_human)) {
    genes_ferret_high <- rownames(expr_ferret)[expr_ferret[, cluster] > 2]
    genes_human_high <- rownames(expr_human)[expr_human[, cluster] > 2]
    common_genes <- intersect(
      intersect(genes_ferret_high, genes_human_high),
      combined_genes_unique
    )
    genes_per_cluster[[cluster]] <- common_genes
  }
}

# Extract all of these genes per cluster into a final list
for (cluster in shared_celltypes) {
  if (cluster %in% names(genes_per_cluster)) {
    all_genes_ordered <- c(all_genes_ordered, genes_per_cluster[[cluster]])
  }
}

# Calculate a 2nd iteration of the ferret heatmap using the identified genes
final_ferret_heatmap <- get_validation_heatmap(
  ferret_object,
  custom_order = shared_celltypes,
  idents = "majority_celltype",
  genes_of_interest = all_genes_ordered,
  cluster_rows = FALSE,
  show_rownames = FALSE
)
pdf(file.path(image_dir, "04_final_ferret_heatmap.pdf"), height = 8, width = 3)
print(final_ferret_heatmap$heatmap)
dev.off()

# Calculate a 2nd iteration of the human heatmap using the identified genes
final_human_heatmap <- get_validation_heatmap(
  human_subset,
  custom_order = shared_celltypes,
  idents = "Celltype",
  genes_of_interest = all_genes_ordered,
  cluster_rows = FALSE,
  show_rownames = FALSE
)
pdf(file.path(image_dir, "04_final_human_heatmap.pdf"), height = 8, width = 3)
print(final_human_heatmap$heatmap)
dev.off()


# Save these genes and their expression data to output csv files
write.csv(
  final_ferret_heatmap$expression_data,
  file = file.path(csv_dir, "heatmap_ferret_expression_data.csv"),
  row.names = TRUE
)

write.csv(
  final_human_heatmap$expression_data,
  file = file.path(csv_dir, "heatmap_human_expression_data.csv"),
  row.names = TRUE
)
