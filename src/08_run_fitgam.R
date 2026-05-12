library(here)
library(Seurat)
library(slingshot)
library(tradeSeq)
set.seed(0)

source(here("src/common_setup.R"))

#' Run the tradeSeq pipeline for a single dataset
#'
#' Loads a Seurat object and a Slingshot SCE object, isolates trajectory cells,
#' identifies variable features, evaluates optimal knot number via evaluateK,
#' and fits a GAM model via fitGAM. Saves the evaluateK plot, icMat object,
#' and fitted SCE object to disk, with filenames derived from sample_id.
#'
#' @param sample_id Character. Identifier for the sample, used to name all
#' output files (e.g. "ferret_ductal_S9E0" produces
#'   evaluateK_plots_ferret_ductal_S9E0.pdf,
#'   icmat_ferret_ductal_S9E0.rds, and
#'   fitgam_ferret_ductal_S9E0.rda).
#' @param nknots Integer. Number of knots to use in fitGAM. Default is 7.
#' @param nfeatures Integer. Number of variable features to identify via
#'   FindVariableFeatures. Default is 3000.
#' @param k_range Integer vector. Range of knot values to evaluate in
#'   evaluateK. Default is 3:10.
#' @param nGenes Integer. Number of genes to use in evaluateK. Default
#'   is 200.
#'
#' @examples
#' run_tradeseq_pipeline(
#'   slingshot_path = "ferret_ductal_S9E0_sce.rds",
#'   sample_id      = "ferret_ductal_S9E0"
#' )
run_tradeseq_pipeline <- function(sample_id,
                                  nknots = 7,
                                  nfeatures = 3000,
                                  k_range = 3:10,
                                  nGenes = 200) {
  # Load in the saved Seurat and slingshot objects
  seurat_object <- readRDS(file.path(
    fitgam_dir,
    paste0(sample_id, ".rds")
  ))
  slingshot_object <- readRDS(file.path(
    fitgam_dir,
    paste0(sample_id, "_sce.rds")
  ))

  # Isolate the cells involved in the slingshot trajectory
  slingshot_cells <- cellNames(slingshot_object$slingshot)
  seurat_object <- subset(seurat_object, cells = slingshot_cells)

  # Find variable features and pull out the counts associated with it
  seurat_object <- FindVariableFeatures(seurat_object, nfeatures = nfeatures)
  counts <- GetAssayData(seurat_object, slot = "data")
  counts <- counts[rownames(counts) %in% VariableFeatures(seurat_object), ]

  # EvaluateK between k_range for best fit
  pdf(file.path(image_dir, paste0("08_evaluateK_plots_", sample_id, ".pdf")))
  icMat <- evaluateK(
    counts = counts,
    sds = slingshot_object$slingshot,
    k = k_range,
    nGenes = nGenes,
    verbose = TRUE
  )
  dev.off()

  # Save icMat object to an output file
  saveRDS(icMat, file.path(fitgam_dir, paste0("icmat_", sample_id, ".rds")))

  # Extract the slingshot pseudotime and cell weights from the object
  pseudotime <- slingPseudotime(slingshot_object$slingshot, na = FALSE)
  cellWeights <- slingCurveWeights(slingshot_object$slingshot)

  # Run fitGAM with these counts, pseudotime, cell weights, and number of knots
  sce <- fitGAM(
    counts = counts,
    pseudotime = pseudotime,
    cellWeights = cellWeights,
    nknots = nknots,
    verbose = TRUE
  )

  # Save the output SingleCellExperiment object to a file
  saveRDS(sce, file.path(fitgam_dir, (paste0("fitgam_", sample_id, ".rda"))))
}

# Run for each ferret dataset
run_tradeseq_pipeline("ferret_ductal_S9E0")
run_tradeseq_pipeline("ferret_acinar_S9E10")

# Run for each human dataset
run_tradeseq_pipeline("human_ductal_S10E19")
run_tradeseq_pipeline("human_acinar_S10E2")
