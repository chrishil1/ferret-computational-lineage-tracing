# Next figure out how to add this to the snakemake pipeline
# Should run all other scripts first
# Maybe add in a "merged_scripts" section
# In getting scripts, if sample is merged, return all other scripts?

library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)
library(here)
library(openxlsx)


# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

args <- commandArgs(trailingOnly = TRUE)

sample <- args[[1]]
sample <- gsub("__.*", "", sample)
#sample <- "merged"

sample_info <- args[[3]]
#sample_info <- here("files/sample_info.tsv")

results_dir <- args[[2]]
#results_dir <- here("results")

sample_info <- read.table(sample_info, fill = TRUE, header = TRUE)

samples_use <- sample_info[sample_info$sample != sample, ]$sample

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
ifelse(!dir.exists(save_dir), dir.create(save_dir), FALSE)

ifelse(!dir.exists(file.path(save_dir, "images")),
       dir.create(file.path(save_dir, "images")), FALSE)

ifelse(!dir.exists(file.path(save_dir, "files")),
       dir.create(file.path(save_dir, "files")), FALSE)

ifelse(!dir.exists(file.path(save_dir, "rda_obj")),
       dir.create(file.path(save_dir, "rda_obj")), FALSE)

all_objs <- lapply(samples_use, function(x){
  so <- readRDS(file.path(results_dir, "R_analysis", x, 
                          "rda_obj", "seurat_processed.rds"))
  
  return(so)
})

# Merge Seurat objects
seurat_data <- merge(all_objs[[1]], all_objs[2:length(all_objs)])

rm(all_objs)

seurat_data$batch <- ifelse(grepl("_b2$", seurat_data$sample), "b2", "b1")
seurat_data$day <- gsub("_b.*", "", seurat_data$sample)

DefaultAssay(seurat_data) <- "RNA"

# Renormalize
seurat_data <- seurat_data %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()


DefaultAssay(seurat_data) <- "RNA"

saveRDS(seurat_data, file = file.path(save_dir, "rda_obj",
                                      "seurat_start.rds"))
