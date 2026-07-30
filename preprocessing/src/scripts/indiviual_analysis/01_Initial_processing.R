library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)
library(scuttle)
library(here)
library(openxlsx)
library(Matrix)

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

vars.to.regress <- NULL

# Set directories
save_dir <- file.path(results_dir, "R_analysis", sample)

# Make directories
ifelse(!dir.exists(save_dir), dir.create(save_dir, recursive = TRUE), FALSE)

ifelse(!dir.exists(file.path(save_dir, "images")),
       dir.create(file.path(save_dir, "images")), FALSE)

ifelse(!dir.exists(file.path(save_dir, "files")),
       dir.create(file.path(save_dir, "files")), FALSE)

ifelse(!dir.exists(file.path(save_dir, "rda_obj")),
       dir.create(file.path(save_dir, "rda_obj")), FALSE)


mito_features <- read.csv(here("files/mito_genes.csv"), header = FALSE)$V1

# Create seurat object
count_matrix <- readMM(file.path(results_dir, sample, "all-sample",
                                 "DGE_filtered", "count_matrix.mtx"))

gene_info <- read.table(file.path(results_dir, sample, "all-sample",
                                  "DGE_filtered", "all_genes.csv"),
                        header = TRUE, sep = ",")

colnames(count_matrix) <- gene_info$gene_id

cell_info <- read.table(file.path(results_dir, sample, "all-sample",
                                  "DGE_filtered", "cell_metadata.csv"),
                        header = TRUE, sep = ",", row.names = 1)

rownames(count_matrix) <- rownames(cell_info)

count_matrix <- t(count_matrix)

seurat_object <- CreateSeuratObject(counts = count_matrix,
                                    meta.data = cell_info)

seurat_object$orig.ident <- seurat_object$sample

Idents(seurat_object) <- seurat_object$orig.ident

# Add mitochondrial percent
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object,
                                                      features = mito_features)

# Quality plots to determin cutoffs
rna_qual <- VlnPlot(seurat_object,
                    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                    ncol = 3)

# Use scuttle for cutoffs ------------------------------------------------------
se <- as.SingleCellExperiment(seurat_object)
is.mito <- which(rownames(se) %in% mito_features)
# per.cell <- perCellQCMetrics(se, subsets=list(Mito=is.mito))
# 
# qc.stats <- perCellQCFilters(per.cell,
#                              sum.field = "sum",
#                              detected.field = "detected",
#                              sub.fields = "subsets_Mito_percent")



per.cell <- perCellQCMetrics(se)

qc.stats <- perCellQCFilters(per.cell,
                             sum.field = "sum",
                             detected.field = "detected")


all_info <- cbind(per.cell, qc.stats)
colnames(all_info) <- paste0("cell_qc_", colnames(all_info))
all_info <- data.frame(all_info)
seurat_object <- AddMetaData(seurat_object, metadata = all_info)

rna_qual <- featDistPlot(seurat_object,
                         geneset = c("nFeature_RNA", "nCount_RNA",
                                     "percent.mt"),
                         plot_type = "violin", col_by = "cell_qc_discard",
                         combine = FALSE)

seurat_object <- seurat_object %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()

# Save before moving on
saveRDS(seurat_object, file = file.path(save_dir, "rda_obj",
                                        "seurat_unfilt.rds"))

seurat_object <- subset(seurat_object, subset = cell_qc_discard == FALSE)

seurat_object <- seurat_object %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()

saveRDS(seurat_object, file = file.path(save_dir, "rda_obj",
                                        "seurat_start.rds"))