#' Update Seurat Object Gene Names with Ortholog Mappings
#'
#' Replaces gene names in a Seurat object with their ortholog equivalents using
#' a provided mapping file. Genes without a mapping retain their original names.
#' Genes that result in duplicate names after mapping are removed entirely.
#'
#' @param seurat_object A Seurat object whose RNA assay gene names will be
#'   updated.
#' @param mapping_file A data frame containing gene-to-ortholog mappings. Must
#'   include a column named gene_id with the original gene identifiers.
#' @param mapping_ortholog_col A string specifying the column name in
#'   mapping_file that contains the ortholog gene identifiers to map to.
#'
#' @return A Seurat object with updated gene names in the RNA assay counts
#'   layer. Genes that would result in duplicate names after mapping are
#'   dropped.
#'
#' @examples
#' mapping <- read.csv("ferret_to_human_orthologs.csv")
#' updated_seurat <- update_with_orthologs(
#'   seurat_object        = my_seurat,
#'   mapping_file         = mapping,
#'   mapping_ortholog_col = "human_gene_id"
#' )
update_with_orthologs <- function(seurat_object,
                                  mapping_file,
                                  mapping_ortholog_col) {
  # Define gene mapping columns
  mapping_gene_col <- "gene_id"

  # Create a mapping list of ferret genes to orthologs
  mapping_list <- setNames(
    mapping_file[[mapping_ortholog_col]],
    mapping_file[[mapping_gene_col]]
  )

  # Replace ferret gene names with orthologs if available
  new_gene_names <- ifelse(rownames(seurat_object) %in% names(mapping_list),
    mapping_list[rownames(seurat_object)],
    rownames(seurat_object)
  )

  # Identify and remove duplicated gene names after mapping
  duplicated_gene_names <- new_gene_names[duplicated(new_gene_names)]
  rows_to_keep <- !new_gene_names %in% duplicated_gene_names

  # Subset the Seurat object to only keep unique gene names
  seurat_object <- seurat_object[rows_to_keep, ]
  new_gene_names <- new_gene_names[rows_to_keep]

  # Update Seurat object with new gene names in the counts data
  count_data <- GetAssayData(
    object = seurat_object,
    assay = "RNA",
    layer = "counts"
  )
  rownames(count_data) <- new_gene_names
  new_assay <- CreateAssayObject(
    counts = count_data,
    min.cells = 0,
    min.features = 0
  )
  seurat_object[["RNA"]] <- new_assay

  return(seurat_object)
}

#' Subset Ferret and Human Seurat Objects by ANOVA-Based Gene Filtering
#'
#' Filters genes that show strong species-specific expression differences
#' between ferret and human Seurat objects. A one-way ANOVA is applied to each
#' gene across both species, and genes with an F statistic above the 90th
#' percentile are removed, retaining only genes with broadly comparable
#' expression patterns.
#'
#' @param ferret_subset A Seurat object containing ferret single-cell RNA data,
#'   with normalized expression in the RNA assay data layer and cell types
#'   stored in final_combined_celltype.
#' @param human_subset A Seurat object containing human single-cell RNA data,
#'   with cell types stored in a Celltype metadata column. Raw counts are
#'   expected, as normalization is applied within this function.
#'
#' @return A named list with two elements:
#'   ferret_anova_subset: The ferret Seurat object subsetted to
#'     genes passing the ANOVA filter.
#'   human_anova_subset: The human Seurat object subsetted to
#'     genes passing the ANOVA filter.
#' }
#'
#' @examples
#' result <- subset_anova_genes(
#'   ferret_subset = ferret_seurat,
#'   human_subset  = human_seurat
#' )
#' ferret_filtered <- result$ferret_anova_subset
#' human_filtered <- result$human_anova_subset
subset_anova_genes <- function(ferret_subset,
                               human_subset) {
  # Normalize the human subset
  human_subset <- NormalizeData(human_subset)

  # Rename the columns for correction step
  human_subset$final_combined_celltype <- human_subset$Celltype
  human_subset$final_combined_celltype <-
    paste0(human_subset$final_combined_celltype, "_human")

  # Assign batch labels for batch correction
  ferret_subset$batch <- "Ferret"
  human_subset$batch <- "Human"

  # Capitalize ferret cell types so they match comparison files
  capitalized_celltypes <- str_to_title(ferret_subset$final_combined_celltype)
  ferret_subset$final_combined_celltype <- capitalized_celltypes

  # Get expression matrices for both subsets
  ferret_expr_matrix <- LayerData(ferret_subset, layer = "data", assay = "RNA")
  human_expr_matrix <- LayerData(human_subset, layer = "data", assay = "RNA")

  # Combine expression matrices for comparison
  combined_expr_matrix <- cbind(ferret_expr_matrix, human_expr_matrix)

  # Create a species factor (matching the columns from ferret and human
  # expression data)
  species <- c(
    rep("ferret", ncol(ferret_expr_matrix)),
    rep("human", ncol(human_expr_matrix))
  )

  # Apply one-way ANOVA to each gene
  anova_results <- apply(combined_expr_matrix, 1, function(gene_expr) {
    species_factor <- factor(species)
    anova(lm(gene_expr ~ species_factor))$F[1] # Extract the F statistic
  })
  anova_results <- anova_results[!is.na(anova_results)]

  # Identify the 90th percentile cutoff for F statistic
  f_cutoff <- quantile(anova_results, 0.90)

  # Filter out genes with an F statistic above the 90th percentile
  # (genes strongly differing by species)
  genes_to_keep <- names(anova_results[anova_results <= f_cutoff])

  # Subset the Seurat objects to retain only the selected genes
  ferret_subset <- subset(ferret_subset, features = genes_to_keep)
  human_subset <- subset(human_subset, features = genes_to_keep)

  return(list(
    ferret_anova_subset = ferret_subset,
    human_anova_subset = human_subset
  ))
}

#' Integrate Ferret and Human Seurat Objects and Prepare for UMAP Visualization
#'
#' Orchestrates a cross-species integration pipeline for ferret and human
#' single-cell RNA data. Gene names in the ferret object are first mapped to
#' human orthologs, after which the two datasets are filtered to shared genes,
#' ANOVA-filtered to remove strongly species-biased genes, merged, and
#' integrated using Seurat's CCA method. The resulting object contains a UMAP
#' embedding ready for visualization.
#'
#' @param ferret_object A Seurat object containing ferret single-cell RNA data
#'   with gene names in ferret gene ID format.
#' @param human_object A Seurat object containing human single-cell RNA data
#'   with gene names in human gene symbol format.
#' @param mapping_file A data frame containing ferret-to-human ortholog mappings
#'   passed to update_with_orthologs. Must include a gene_id column
#'   and a Human.gene.name column.
#'
#' @return A merged and integrated Seurat object containing both ferret and
#'   human cells, with a CCA-based PCA reduction (integrated.cca) and a UMAP
#'   embedding (umap) computed across 50 dimensions. Cell type labels are
#'   stored in final_combined_celltype.
#'
#' @examples
#' mapping <- read.csv("ferret_to_human_orthologs.csv")
#' integrated_seurat <- plot_human_integration(
#'   ferret_object = ferret_seurat,
#'   human_object  = human_seurat,
#'   mapping_file  = mapping
#' )
calculate_human_integration <- function(ferret_object,
                                        human_object,
                                        mapping_file) {
  # Update the Seurat object with orthologs, if found
  ferret_object <- update_with_orthologs(
    seurat_object = ferret_object,
    mapping_file = mapping_file,
    mapping_ortholog_col = "Human.gene.name"
  )

  # Intersect the rownames, only keep those included in both
  ferret_mapped_genes <- rownames(ferret_object)
  human_genes <- rownames(human_object)
  intersected_genes <- intersect(ferret_mapped_genes, human_genes)

  # Subset both ferret and human Seurat objects to keep intersected genes
  ferret_subset <- ferret_object[intersected_genes, ]
  human_subset <- human_object[intersected_genes, ]

  # Subset out Seurat objects to remove highly biasing genes
  anova_subset_objects <- subset_anova_genes(ferret_subset, human_subset)

  # Merge ferret and human Seurat objects for integration
  merged_seurat <- merge(anova_subset_objects$ferret_anova_subset,
    y = anova_subset_objects$human_anova_subset,
    add.cell.ids = c("ferret", "human")
  )

  # Perform data normalization, variable feature selection, and scaling
  merged_seurat <- NormalizeData(merged_seurat)
  merged_seurat <- FindVariableFeatures(merged_seurat)
  merged_seurat <- ScaleData(merged_seurat)
  merged_seurat <- RunPCA(merged_seurat, npcs = 50)

  # Perform integration on merged object with Seurats CCAIntegration method
  merged_seurat <- IntegrateLayers(
    object = merged_seurat, method = CCAIntegration,
    orig.reduction = "pca", new.reduction = "integrated.cca",
    verbose = FALSE
  )
  merged_seurat <- RunUMAP(merged_seurat,
    reduction = "integrated.cca",
    dims = 1:50,
    reduction.name = "umap"
  )

  # Shorten cell type name for ease of reading on downstream plots
  merged_seurat$final_combined_celltype <-
    gsub(
      "Fev_endocrine_progenitor",
      "Fev EP",
      merged_seurat$final_combined_celltype
    )

  return(merged_seurat)
}

#' Reannotate Ferret Cells Using kNN Mapping from Human Cell Types
#'
#' Uses k-nearest neighbours (kNN) in UMAP space to transfer human cell type
#' labels onto ferret cells in an integrated Seurat object. For each ferret
#' cell, the majority cell type among its 10 nearest human neighbors is assigned
#' as the new label. Original cell type labels are preserved in a separate
#' metadata column before reassignment.
#'
#' @param merged_seurat A merged and integrated Seurat object produced by
#'  plot_human_integration, containing both ferret and human cells with a UMAP
#'  embedding and cell type labels in final_combined_celltype.
#' @param mapping_file A dataframe containing ferret-to-human ortholog mappings.
#'
#' @return The input Seurat object with updated final_combined_celltype
#'   metadata for ferret cells, relabeled based on kNN voting from neighboring
#'   human cells. The _human suffix is stripped from transferred labels.
#'   Original labels are retained in the original_labels metadata column.
#'
#' @examples
#' mapping <- read.csv("ferret_to_human_orthologs.csv")
#' annotated_seurat <- calculate_human_knn_mapping(
#'   merged_seurat = integrated_seurat,
#'   mapping_file  = mapping
#' )
calculate_human_knn_mapping <- function(merged_seurat,
                                        mapping_file) {
  # Save the original labels from integetion in case they are needed downstream
  merged_seurat$original_labels <- merged_seurat$final_combined_celltype

  # Extract UMAP embeddings and cell indices for ferret and human cells
  umap_coords <- Embeddings(merged_seurat, reduction = "umap")

  # Pull the ferret specific cells from the merged object
  ferret_pulled_cells <-
    WhichCells(merged_seurat,
      expression = final_combined_celltype %in% ferret_cells
    )

  # Pull the human specific cells from the merged object
  human_pulled_cells <-
    WhichCells(merged_seurat,
      expression = final_combined_celltype %in% human_cells
    )

  # Get UMAP coordinates for ferret and human cells
  ferret_umap_coords <- umap_coords[ferret_pulled_cells, ]
  human_umap_coords <- umap_coords[human_pulled_cells, ]

  # Apply k-nearest neighbors (kNN) to map ferret cells to nearest neighbors
  nn_result <- get.knnx(human_umap_coords, ferret_umap_coords, k = 10)

  # Initialize an empty vector to store new labels for ferret cells
  new_labels <- vector("character", length(ferret_pulled_cells))

  # For each ferret cell, find the majority cell type among nearest neighbors
  for (i in seq_along(ferret_pulled_cells)) {
    neighbor_indices <- nn_result$nn.index[i, ]
    neighbor_types <-
      merged_seurat@meta.data[
        human_pulled_cells[neighbor_indices],
        "final_combined_celltype"
      ]
    new_labels[i] <- names(sort(table(neighbor_types), decreasing = TRUE))[1]
  }

  # Initialise knn_celltype for all cells, then populate ferret cells with
  # their kNN-derived labels (with the _human suffix stripped)
  merged_seurat@meta.data$knn_celltype <- NA_character_
  merged_seurat@meta.data[ferret_pulled_cells, "knn_celltype"] <- gsub(
    "_human$",
    "",
    new_labels
  )

  return(merged_seurat)
}
