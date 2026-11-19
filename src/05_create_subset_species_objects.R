library(biomaRt)
library(dplyr)
library(here)
library(Seurat)
library(org.Hs.eg.db)
library(rtracklayer)
library(dplyr)

source(here("src/common_setup.R"))
options(Seurat.object.assay.version = "v3")

# Establish mapping dictionary from ortholog file
mapping_gene_col <- "gene_id"
mapping_ortholog_col <- "Human.gene.name"
mapping_list <- setNames(
  mapping_file[[mapping_ortholog_col]],
  mapping_file[[mapping_gene_col]]
)

# Load in our Seurat objects
ferret_object <- readRDS(file.path(data_dir, "03_combined_ferret_res_7.rds"))
human_object <- readRDS(file.path(data_dir, "04_human_remade_seurat.rds"))
Idents(human_object) <- "Leiden_cluster"

# Clean up Seurat object after reclustering took place
ferret_object[["RNA"]]$data <- NULL
ferret_object[["RNA"]]$data.1 <- NULL
ferret_object[["RNA"]]$scale.data <- NULL
ferret_object[["RNA"]]$counts <- ferret_object[["RNA"]]$counts.1
ferret_object[["RNA"]]$counts.1 <- NULL

# Replace ferret gene names with human orthologs, if available
new_gene_names <- ifelse(rownames(ferret_object) %in% names(mapping_list) &
  mapping_list[rownames(ferret_object)] != "",
mapping_list[rownames(ferret_object)],
rownames(ferret_object)
)

# Remove the names from the new gene mapping
new_gene_names <- unname(new_gene_names)

# Identify and remove duplicated gene names after mapping
duplicated_gene_names <- new_gene_names[duplicated(new_gene_names)]
rows_to_keep <- !new_gene_names %in% duplicated_gene_names

# Subset the Seurat object to only keep unique gene names
ferret_object <- ferret_object[rows_to_keep, ]
new_gene_names <- new_gene_names[rows_to_keep]

# Extract the counts layer from the ferret object
count_data <- GetAssayData(
  object = ferret_object,
  assay = "RNA",
  layer = "counts"
)

# Rename the rows with the new mappings
rownames(count_data) <- new_gene_names

# Create a new object with all the metadata moved over
new_ferret_object <- CreateSeuratObject(
  counts = count_data,
  meta.data = ferret_object@meta.data
)

# Run standard preprocesing and dimensionality reduction on new ferret object
new_ferret_object <- NormalizeData(new_ferret_object)
new_ferret_object <- FindVariableFeatures(new_ferret_object)
new_ferret_object <- ScaleData(new_ferret_object)
new_ferret_object <- RunPCA(new_ferret_object, npcs = 20)
new_ferret_object <- FindNeighbors(new_ferret_object, dims = 1:20)
new_ferret_object <- FindClusters(new_ferret_object, resolution = 0.8)
new_ferret_object <- RunUMAP(new_ferret_object, reduction = "pca", dims = 1:20)

gtf_file <- here("Homo_sapiens.GRCh38.95.gtf")
gtf_data <- import(gtf_file)
gtf_df <- as.data.frame(gtf_data)

# Extract unique gene_id and gene_name (or gene_symbol) pairs
gene_mapping <- gtf_df %>%
  dplyr::select(gene_id, gene_name) %>%
  distinct() %>%
  filter(!is.na(gene_id) & !is.na(gene_name))

# Extract Ensembl IDs from the human Seurat object
ensg_ids <- rownames(human_object)

# Match Ensembl IDs to gene symbols
matched_gene_names <- gene_mapping$gene_name[
  match(ensg_ids, gene_mapping$gene_id)
]

# Use gene symbol if available; otherwise retain Ensembl ID
new_rownames <- ifelse(
  !is.na(matched_gene_names) & matched_gene_names != "",
  matched_gene_names,
  ensg_ids
)

# Remove names attribute, if present
new_rownames <- unname(new_rownames)

# Extract counts data from SCTransformed human data
count_data <- GetAssayData(
  object = human_object,
  assay = "SCT",
  layer = "counts"
)

# Apply the new rownames to the count matrix
rownames(count_data) <- new_rownames

# Extract embeddings from human seurat object
pca_embeddings <- Embeddings(human_object[["pca"]])
umap_embeddings <- Embeddings(human_object[["umap"]])

# Create a new object with these updated counts and metadata
new_human_object <- CreateSeuratObject(
  counts = count_data,
  meta.data = human_object@meta.data
)

# Add the PCA embeddings back into the object
new_human_object[["pca"]] <- CreateDimReducObject(
  embeddings = pca_embeddings,
  key = "PC_",
  assay = DefaultAssay(new_human_object)
)

# Add the UMAP embeddings back into the object
new_human_object[["umap"]] <- CreateDimReducObject(
  embeddings = umap_embeddings,
  key = "UMAP_",
  assay = DefaultAssay(new_human_object)
)

# Normalize and find variable features
new_human_object <- NormalizeData(new_human_object)
new_human_object <- FindVariableFeatures(new_human_object)

# Get the common genes for both and remove everything else
ferret_genes <- rownames(new_ferret_object)
human_genes <- rownames(new_human_object)
combined_genes <- intersect(ferret_genes, human_genes)

# Get the final ferret object to be used in downstream analyses
final_ferret <- subset(new_ferret_object, features = combined_genes)

# Convert the new ferret object back into the typical Seurat structure
final_ferret[["RNA3"]] <- as(object = final_ferret[["RNA"]], Class = "Assay")
DefaultAssay(final_ferret) <- "RNA3"
final_ferret[["RNA"]] <- NULL
final_ferret <- RenameAssays(object = final_ferret, RNA3 = "RNA")
final_ferret[["RNA"]]$data <- NULL
final_ferret[["RNA"]]$scale.data <- NULL

final_ferret <- NormalizeData(final_ferret)
final_ferret <- FindVariableFeatures(final_ferret)
final_ferret <- ScaleData(final_ferret)
final_ferret <- RunPCA(final_ferret, npcs = 20)
final_ferret <- FindNeighbors(final_ferret, dims = 1:20)
final_ferret <- FindClusters(final_ferret, resolution = 3.5)
final_ferret <- RunUMAP(final_ferret, reduction = "pca", dims = 1:20)

# Set clustering to RNA resolution 3.5
final_ferret$seurat_clusters <- final_ferret$RNA_snn_res.3.5
Idents(final_ferret) <- "seurat_clusters"

saveRDS(final_ferret, file.path(data_dir, "05_final_ferret.rds"))

# Get the final human object to be used in downstream analyses
final_human <- subset(new_human_object, features = combined_genes)
saveRDS(final_human, file.path(data_dir, "05_final_human.rds"))

# Determine the background genes by just taking all genes present in both
background_genes <- data.frame(
  Genes = rownames(final_human)
)

# Write the background genes to an output file
write.table(
  background_genes,
  file = here("csv_files/background_genes.csv"),
  quote = FALSE,
  row.names = FALSE,
  sep = ","
)
