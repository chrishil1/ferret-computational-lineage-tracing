library(biomaRt)
library(dplyr)
library(here)
library(Seurat)

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
new_ferret_object <- RunPCA(new_ferret_object, npcs = 15)
new_ferret_object <- FindNeighbors(new_ferret_object, dims = 1:15)
new_ferret_object <- FindClusters(new_ferret_object, resolution = 0.8)
new_ferret_object <- RunUMAP(new_ferret_object, reduction = "pca", dims = 1:15)

# Convert the human ENSG codes to gene IDs
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensg_ids <- rownames(human_object)

# Get the biomart mapping for the human ensembl and gene symbols
gene_mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ensg_ids,
  mart = mart
)

# Create a mapping of the match gene names for the human
matched_gene_names <-
  gene_mapping$hgnc_symbol[match(
    rownames(human_object),
    gene_mapping$ensembl_gene_id
  )]

# Apply the mapping if a non-empty entry exists
new_rownames <- ifelse(!is.na(matched_gene_names) & matched_gene_names != "",
  matched_gene_names,
  rownames(human_object)
)

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
