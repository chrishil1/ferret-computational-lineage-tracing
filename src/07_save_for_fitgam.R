options(Seurat.object.assay.version = "v3")

library(here)
library(Seurat)
library(slingshot)

source(here("src/common_setup.R"))

#' Run Slingshot trajectory inference on a subset of Seurat clusters
#'
#' This function subsets a Seurat object to specified clusters, converts
#' the subset to a SingleCellExperiment object, runs Slingshot trajectory
#' inference using PCA embeddings, and saves both the Slingshot object and
#' subsetted Seurat object to disk.
#'
#' @param seurat_obj A Seurat object containing PCA embeddings and cluster
#'   annotations.
#' @param cluster_col Character string specifying the metadata column used
#'   for cluster labels (e.g., "seurat_clusters" or "Leiden_cluster").
#' @param clusters Vector of cluster IDs to retain for trajectory analysis.
#' @param start_cluster Cluster ID used as the Slingshot starting point.
#' @param end_cluster Cluster ID used as the Slingshot endpoint.
#' @param output_prefix Character string used as the output filename prefix.
#'   Two files will be written:
#'   `{output_prefix}_sce.rds` and `{output_prefix}.rds`.
#' @param fitgam_dir Directory where output RDS files will be saved.
#' @param n_pcs Number of PCA dimensions to use for Slingshot. Default is 15.
#'
#' @examples
#' run_slingshot_subset(
#'   seurat_obj = ferret_subset,
#'   cluster_col = "seurat_clusters",
#'   clusters = c(9, 0, 5, 12, 17),
#'   start_cluster = 9,
#'   end_cluster = 0,
#'   output_prefix = "ferret_ductal_S9E0",
#'   fitgam_dir = fitgam_dir
#' )
run_slingshot_subset <- function(
  seurat_obj,
  cluster_col,
  clusters,
  start_cluster,
  end_cluster,
  output_prefix,
  fitgam_dir,
  n_pcs = 15
) {
  # Define the cells to keep from the object
  cells_keep <- rownames(
    seurat_obj@meta.data[
      seurat_obj@meta.data[[cluster_col]] %in% clusters,
      ,
      drop = FALSE
    ]
  )
  
  # Subset to desired clusters
  cluster_subset <- subset(
    seurat_obj,
    cells = cells_keep
  )

  # Convert to SingleCellExperiment
  sce <- as.SingleCellExperiment(cluster_subset)

  # Subset PCA dimensions
  pca_subset <- reducedDim(sce, "PCA")[, seq_len(n_pcs)]

  # Run slingshot
  slingshot_sce <- slingshot(
    sce,
    reducedDim = pca_subset,
    clusterLabels = colData(sce)[[cluster_col]],
    start.clus = start_cluster,
    end.clus = end_cluster
  )

  # Save outputs
  saveRDS(
    slingshot_sce,
    file.path(fitgam_dir, paste0(output_prefix, "_sce.rds"))
  )

  saveRDS(
    cluster_subset,
    file.path(fitgam_dir, paste0(output_prefix, ".rds"))
  )
}

#################
# Ferret saving #
#################
ferret_subset <- readRDS(
  file.path(data_dir, "05_final_ferret.rds")
)

run_slingshot_subset(
  seurat_obj = ferret_subset,
  cluster_col = "seurat_clusters",
  clusters = c(9, 0, 5, 12, 17),
  start_cluster = 9,
  end_cluster = 0,
  output_prefix = "ferret_ductal_S9E0",
  fitgam_dir = fitgam_dir
)

run_slingshot_subset(
  seurat_obj = ferret_subset,
  cluster_col = "seurat_clusters",
  clusters = c(9, 4, 14, 1, 16, 8, 11, 10),
  start_cluster = 9,
  end_cluster = 10,
  output_prefix = "ferret_acinar_S9E10",
  fitgam_dir = fitgam_dir
)

################
# Human saving #
################
human_object <- readRDS(
  file.path(data_dir, "04_human_remade_seurat.rds")
)
Idents(human_object) <- "Leiden_cluster"

run_slingshot_subset(
  seurat_obj = human_object,
  cluster_col = "Leiden_cluster",
  clusters = c(5, 8, 10, 19, 21, 27),
  start_cluster = 10,
  end_cluster = 19,
  output_prefix = "human_ductal_S10E19",
  fitgam_dir = fitgam_dir
)

run_slingshot_subset(
  seurat_obj = human_object,
  cluster_col = "Leiden_cluster",
  clusters = c(2, 3, 4, 6, 10, 11),
  start_cluster = 10,
  end_cluster = 2,
  output_prefix = "human_acinar_S10E2",
  fitgam_dir = fitgam_dir
)
