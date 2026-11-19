library(dplyr)
library(here)
library(JASPAR2020)
library(pheatmap)
library(Seurat)
library(TFBSTools)

source(here("src/common_setup.R"))
source(here("src/06_celltype_validation_helpers.R"))
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

# Define the cell types that are found in both datasets
shared_celltypes <- c(
  "Acinar Immature", "Acinar Maturing", "Alpha", "Beta", "Delta",
  "Ductal Immature", "Ductal Proliferating", "Acinar Proliferating",
  "Pancreatic Progenitor"
)

# Load in Seurat objects for the two species
ferret_object <- readRDS(file.path(data_dir, "05_final_ferret.rds"))

human_object <- readRDS(file.path(data_dir, "05_final_human.rds"))
human_subset <- subset(human_object, subset = Celltype %in% shared_celltypes)

ferret_genes <- FindAllMarkers(ferret_object, group.by = "final_renamed_celltype")
human_genes <- FindAllMarkers(human_object, group.by = "Celltype")

# Function to convert cell types into suffixes for variable naming ease
make_varname <- function(celltype) gsub(" ", "_", tolower(celltype))

# Loop over all cell types
for (celltype in shared_celltypes) {
  suffix <- make_varname(celltype)

  # Subset Seurat objects based on the current cell type of interest
  assign(paste0("ferret_", suffix), subset(
    ferret_object,
    final_renamed_celltype == celltype
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

  intersected <- intersect(
    ferret_celltype_genes$gene,
    human_celltype_genes$gene
  )

  assign(
    paste0("intersected_genes_", suffix),
    intersected
  )
}

combined_gene_list <- lapply(shared_celltypes, function(celltype) {
  get(paste0("intersected_genes_", make_varname(celltype)))
})
names(combined_gene_list) <- make_varname(shared_celltypes)

combined_genes_unique <- unique(
  unlist(combined_gene_list, use.names = FALSE)
)

# Calculate the initial ferret heatmap including all genes extracted
initial_ferret_heatmap <- get_validation_heatmap(
  ferret_object,
  custom_order = shared_celltypes,
  idents = "final_renamed_celltype",
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

# Minimum z-score required in both species
z_threshold <- 1.5

# Maximum number of genes per cell type for final heatmaps
top_n <- 50

genes_per_cluster <- setNames(
  lapply(shared_celltypes, function(ct) {
    # Only consider shared DE markers for this cell type
    candidates <- combined_gene_list[[make_varname(ct)]]

    # Keep genes present in both expression matrices
    candidates <- intersect(
      candidates,
      intersect(rownames(expr_ferret), rownames(expr_human))
    )

    # Require enrichment in BOTH species
    keep <- is.finite(expr_ferret[candidates, ct]) &
      is.finite(expr_human[candidates, ct]) &
      expr_ferret[candidates, ct] > z_threshold &
      expr_human[candidates, ct] > z_threshold

    candidates[keep]
  }),
  shared_celltypes
)

# Rank passing genes by conserved log2FC
ranked_tables <- setNames(
  lapply(shared_celltypes, function(ct) {
    genes <- genes_per_cluster[[ct]]

    # Ferret differential expression
    f <- ferret_genes %>%
      filter(.data$cluster == ct, gene %in% genes) %>%
      distinct(gene, .keep_all = TRUE) %>%
      transmute(gene, Ferret_FC = avg_log2FC)

    # Human differential expression
    h <- human_genes %>%
      filter(.data$cluster == ct, gene %in% genes) %>%
      distinct(gene, .keep_all = TRUE) %>%
      transmute(gene, Human_FC = avg_log2FC)

    # Rank by minimum fold change across species
    inner_join(f, h, by = "gene") %>%
      filter(Ferret_FC > 0.5, Human_FC > 0.5) %>%
      mutate(Score = pmin(Ferret_FC, Human_FC)) %>%
      arrange(desc(Score), gene) %>%
      mutate(Rank = row_number())
  }),
  shared_celltypes
)

# Keep the top N genes per cell type
ranked_genes_per_cluster <- lapply(
  ranked_tables,
  function(tbl) head(tbl$gene, top_n)
)

# Combine selected genes for final heatmaps
all_genes_ordered <- unique(unlist(
  ranked_genes_per_cluster,
  use.names = FALSE
))

# Calculate a 2nd iteration of the ferret heatmap using the identified genes
final_ferret_heatmap <- get_validation_heatmap(
  ferret_object,
  custom_order = shared_celltypes,
  idents = "final_renamed_celltype",
  genes_of_interest = all_genes_ordered,
  cluster_rows = FALSE,
  show_rownames = FALSE
)
pdf(file.path(image_dir, "06_final_ferret_heatmap.pdf"), height = 8, width = 3)
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
pdf(file.path(image_dir, "06_final_human_heatmap.pdf"), height = 8, width = 3)
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
