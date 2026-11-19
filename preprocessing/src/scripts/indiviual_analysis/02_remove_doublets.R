library(DoubletFinder)
library(Seurat)
library(tidyverse)
library(here)
library(scAnalysisR)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

args <- commandArgs(trailingOnly = TRUE)


sample <- args[[1]]
sample <- gsub("__.*", "", sample)
#sample <- "combined_group_2"

sample_info <- args[[3]]
#sample_info <- here("files/sample_info.tsv")

results_dir <- args[[2]]
#results_dir <- here("results")

sample_info <- read.table(sample_info, fill = TRUE, header = TRUE)

sample_info <- sample_info[sample_info$sample == sample,]

HTO <- sample_info$HTO
ADT <- sample_info$ADT
hash_ident <- sample_info$hash_ident


if(normalization_method == "SCT"){
  SCT <- TRUE
  seurat_assay <- "SCT"
} else {
  SCT <- FALSE
  seurat_assay <- "RNA"
}

# Set directories
save_dir <- file.path(results_dir, "R_analysis", sample)

# Read in data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_start.rds"))

# Remove "negatives"
if(HTO){
  Idents(seurat_data) <- "HTO_classification.global"
  seurat_data <- subset(x = seurat_data, idents = c("Singlet", "Doublet"))
}

# PCA and UMAP
set.seed(0)
DefaultAssay(seurat_data) <- seurat_assay
seurat_data <- PCA_dimRed(seurat_data, assay = seurat_assay)
npcs <- 20

umap_data <- group_cells(seurat_data, assay = seurat_assay, nPCs = npcs)

seurat_data <- umap_data[[1]]

plotDimRed(seurat_data, col_by = "sample", plot_type = "rna.umap")

# Run Doublet finder -----------------------------------------------------------
set.seed(0)
sweep.res.seurat_data <- paramSweep_v3(seurat_data, PCs = 1:npcs, sct = SCT)

if(HTO){
  ## pK Identification (ground-truth)
  gt.calls <- seurat_data@meta.data[rownames(sweep.res.seurat_data[[1]]),
                                    "HTO_classification.global"]
  gt.calls <- as.factor(gt.calls)
  sweep.stats_sample <- summarizeSweep(sweep.res.seurat_data, GT = TRUE,
                                       GT.calls = gt.calls)
  pdf(file.path(save_dir, "images", "doublet_finder_peak.pdf"))
  
  bcmvn_sample <- find.pK(sweep.stats_sample)
  dev.off()
  max_AUC <- max(bcmvn_sample$MeanAUC)
  max_BC <- max(bcmvn_sample$MeanBC)
  max_BCmetric <- max(bcmvn_sample$BCmetric)
} else {
  ## pK Identification (no ground-truth)
  sweep.stats_sample <- summarizeSweep(sweep.res.seurat_data, GT = FALSE)
  
  pdf(file.path(save_dir, "images", "doublet_finder_peak.pdf"))
  
  bcmvn_sample <- find.pK(sweep.stats_sample)
  dev.off()
  max_BC <- max(bcmvn_sample$MeanBC)
  max_BCmetric <- max(bcmvn_sample$BCmetric)
}


# Pull out best pK
best_pk <- bcmvn_sample[bcmvn_sample$BCmetric == max_BCmetric,]$pK

pK <- as.double(as.character(best_pk))

## Homotypic Doublet Proportion Estimate ---------------------------------------
annotations <- seurat_data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 

# Find expected proportion -----------------------------------------------------

# Retrieve the closest number from the vector
expected_proportion <- 0.02 # Parse generally has low doublet rate

nExp_poi <- round(0.032*nrow(seurat_data@meta.data))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


## Find pK from bcmvn output, pN selection is less important--------------------

seurat_data <- doubletFinder_v3(seurat_data, PCs = 1:npcs, pN = 0.25,
                                pK = pK, nExp = nExp_poi.adj,
                                reuse.pANN = FALSE, sct = SCT)

# Rename to be a consistent name between all samples
df_column <- colnames(seurat_data[[]])[grep("DF.classifications",
                                            colnames(seurat_data[[]]))]


seurat_data$Doublet_finder <- seurat_data[[df_column]]

# Remove the original name
seurat_data[[df_column]] <- NULL

pdf(file.path(save_dir, "images", "doublet_finder.pdf"))

print(plotDimRed(seurat_data, col_by = "Doublet_finder",
                 plot_type = "rna.umap"))
print(plotDimRed(seurat_data, col_by = "JCHAIN", plot_type = "rna.umap"))

dev.off()

saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_doublet.rds"))