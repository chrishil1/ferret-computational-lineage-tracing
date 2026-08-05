library(Seurat)
library(tidyverse)
library(cowplot)
library(harmony)
library(here)
library(scAnalysisR)
library(viridis)
library(clustifyr)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

source(here("src/scripts/functions.R"))

remove_ambience <- FALSE

normalization_method <- "log" # can be SCT or log

cor_cutoff <- 0.5

all_ref_dir <-
  "/beevol/home/wellskri/Analysis/references/single_cell_references"
mapping_file <- read.csv(here("files/species_mapping_file.csv"))

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


ADT_pca <- sample_info$adt_pca
run_adt_umap <- sample_info$adt_umap

RNA_pcs <- sample_info$PCs
resolution <- sample_info$resolution
batch_correction <- sample_info$batch_correction

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
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))

pca_list <- list("rna" = c(assay = seurat_assay, reduction_name = "rna"))

if(remove_ambience){
  pca_list <- c(pca_list,
                list("ambience" = c(assay = "AMBRNA", 
                                    reduction_name = "ambpca")))
}

#-------------------------------------------------------------------------------
for(assay_type in names(pca_list)){
  seurat_assay <- pca_list[[assay_type]][["assay"]]
  reduction_name <- pca_list[[assay_type]][["reduction_name"]]

  if(batch_correction == "none"){
    reduction_save <- reduction_name
  } else {
    reduction_save <- paste(assay_type, batch_correction, sep = "_")
  }
  
  plot_type <- paste0(reduction_save, ".umap")
  
  
  # NOTE update the section below to include any references pertanent to your
  # sample.
  
  # Mouse ------------------------------------------------------------------------
  
  ########################
  # Pancreas development #
  ########################
  
  # Information for cell mapping
  ref_dir <- file.path(all_ref_dir, "pancreas/byrnes_2018_mouse")
  
  ref_mat <- read.csv(file.path(ref_dir, "E14_epithelial_average.csv"),
                      header = TRUE, row.names = 1)
  
  
  cluster_res <- clustifyr_orthologs(seurat_data, ref_mat,
                                     save_dir = save_dir,
                                     save_name = "celltype_byrnes",
                                     mapping_file = mapping_file,
                                     mapping_gene_col = "gene_id",
                                     mapping_ortholog_col = "Mouse.gene.name",
                                     assay = "RNA",
                                     nfeatures = 2000, clusters = "RNA_cluster",
                                     plot_type = plot_type,
                                     cor_cutoff = cor_cutoff)
  
  seurat_data <- cluster_res$object
  
  seurat_res_byrnes <- cluster_res$RNA
  
  plotDimRed(seurat_data, col_by = "RNA_celltype_byrnes",
             plot_type = plot_type)
  
  write.csv(seurat_res_byrnes, file.path(save_dir,
                                         "files/celltype_mapping_byrnes_2018.csv"))
  
  
  #-------------------------------------------------------------------------------
  
  ##################
  # Pancreas atlas #
  ##################
  
  # Information for cell mapping
  ref_dir <- file.path(all_ref_dir, "pancreas/Baron_2016")
  
  ref_mat <- read.csv(file.path(ref_dir, "clustifyr_mouse_reference.csv"),
                      header = TRUE, row.names = 1)
  
  cluster_res <- clustifyr_orthologs(seurat_data, ref_mat,
                                     save_dir = save_dir,
                                     save_name = "baron_celltype",
                                     mapping_file = mapping_file,
                                     mapping_gene_col = "gene_id",
                                     mapping_ortholog_col = "Mouse.gene.name",
                                     assay = "RNA",
                                     nfeatures = 2000, clusters = "RNA_cluster",
                                     plot_type = plot_type,
                                     cor_cutoff = cor_cutoff)
  
  
  seurat_data <- cluster_res$object
  
  seurat_res_baron <- cluster_res$RNA
  
  plotDimRed(seurat_data, col_by = "RNA_baron_celltype", plot_type = plot_type)
  
  write.csv(seurat_res_baron,
            file.path(save_dir, "files/celltype_mapping_baron_2016.csv"))
  
  #-------------------------------------------------------------------------------
  
  ################
  # Tabula muris #
  ################
  
  # Information for cell mapping
  ref_dir <- file.path(all_ref_dir, "pancreas/tabula_muris")
  
  ref_mat <- read.csv(file.path(ref_dir, "clustifyr_reference.csv"),
                      header = TRUE, row.names = 1)
  
  cluster_res <- clustifyr_orthologs(seurat_data, ref_mat,
                                     save_dir = save_dir,
                                     save_name = "tabula_muris_celltype",
                                     mapping_file = mapping_file,
                                     mapping_gene_col = "gene_id",
                                     mapping_ortholog_col = "Mouse.gene.name",
                                     assay = "RNA",
                                     nfeatures = 2000, clusters = "RNA_cluster",
                                     plot_type = plot_type,
                                     cor_cutoff = cor_cutoff)
  
  
  seurat_data <- cluster_res$object
  
  seurat_res_tabula_muris <- cluster_res$RNA
  
  plotDimRed(seurat_data, col_by = "RNA_tabula_muris_celltype",
             plot_type = plot_type)
  
  write.csv(seurat_res_tabula_muris,
            file.path(save_dir, "files/celltype_mapping_tabula_muris.csv"))
  
  
  # Human ------------------------------------------------------------------------
  
  
  ##################
  # Pancreas atlas #
  ##################
  
  # Information for cell mapping
  ref_dir <- file.path(all_ref_dir, "pancreas/Baron_2016")
  
  ref_mat <- read.csv(file.path(ref_dir, "clustifyr_human_reference.csv"),
                      header = TRUE, row.names = 1)
  
  cluster_res <- clustifyr_orthologs(seurat_data, ref_mat,
                                     save_dir = save_dir,
                                     save_name = "baron_human_celltype",
                                     mapping_file = mapping_file,
                                     mapping_gene_col = "gene_id",
                                     mapping_ortholog_col = "Human.gene.name",
                                     assay = "RNA",
                                     nfeatures = 2000, clusters = "RNA_cluster",
                                     plot_type = plot_type,
                                     cor_cutoff = cor_cutoff)
  
  
  seurat_data <- cluster_res$object
  
  seurat_res_baron_human <- cluster_res$RNA
  
  plotDimRed(seurat_data, col_by = "RNA_baron_human_celltype",
             plot_type = plot_type)
  
  write.csv(seurat_res_baron_human,
            file.path(save_dir, "files/celltype_mapping_baron_human_2016.csv"))
  
  
  #-------------------------------------------------------------------------------
  
  ####################
  # Ductal reference #
  ####################
  
  # Information for cell mapping
  ref_dir <- file.path(all_ref_dir, "pancreas/Qadir_2020_human")
  
  ref_mat <- read.csv(file.path(ref_dir, "clustifyr_reference.csv"),
                      header = TRUE, row.names = 1)
  
  cluster_res <- clustifyr_orthologs(seurat_data, ref_mat,
                                     save_dir = save_dir,
                                     save_name = "qadir_celltype",
                                     mapping_file = mapping_file,
                                     mapping_gene_col = "gene_id",
                                     mapping_ortholog_col = "Human.gene.name",
                                     assay = "RNA",
                                     nfeatures = 2000, clusters = "RNA_cluster",
                                     plot_type = plot_type,
                                     cor_cutoff = cor_cutoff)
  
  
  seurat_data <- cluster_res$object
  
  seurat_res_qadir <- cluster_res$RNA
  
  plotDimRed(seurat_data, col_by = "RNA_qadir_celltype",
             plot_type = plot_type)
  
  write.csv(seurat_res_qadir,
            file.path(save_dir, "files/celltype_mapping_qadir_2020.csv"))
  
  
  
  #-------------------------------------------------------------------------------
  
  ##################
  # Pancreas atlas #
  ##################
  
  # Information for cell mapping
  ref_dir <- file.path(all_ref_dir, "pancreas/Muraro_2016_human")
  
  ref_mat <- read.csv(file.path(ref_dir, "clustifyr_reference.csv"),
                      header = TRUE, row.names = 1)
  
  cluster_res <- clustifyr_orthologs(seurat_data, ref_mat,
                                     save_dir = save_dir,
                                     save_name = "muraro_celltype",
                                     mapping_file = mapping_file,
                                     mapping_gene_col = "gene_id",
                                     mapping_ortholog_col = "Human.gene.name",
                                     assay = "RNA",
                                     nfeatures = 2000, clusters = "RNA_cluster",
                                     plot_type = plot_type,
                                     cor_cutoff = cor_cutoff)
  
  
  seurat_data <- cluster_res$object
  
  seurat_res_muraro <- cluster_res$RNA
  
  plotDimRed(seurat_data, col_by = "RNA_muraro_celltype",
             plot_type = plot_type)
  
  write.csv(seurat_res_muraro,
            file.path(save_dir, "files/celltype_mapping_muraro_2016.csv"))
  
  
  #-------------------------------------------------------------------------------
  
  ########################
  # Pancreas development #
  ########################
  
  # Information for cell mapping
  ref_dir <- file.path(all_ref_dir, "pancreas/krentz_2018_human_mouse")
  
  ref_mat <- read.csv(file.path(ref_dir,
                                "S6D1_GFP_clustifyr_reference_celltype.csv"),
                      header = TRUE, row.names = 1)
  
  cluster_res <- clustifyr_orthologs(seurat_data, ref_mat,
                                     save_dir = save_dir,
                                     save_name = "krentz_celltype",
                                     mapping_file = mapping_file,
                                     mapping_gene_col = "gene_id",
                                     mapping_ortholog_col = "Human.gene.name",
                                     assay = "RNA",
                                     nfeatures = 2500, clusters = "RNA_cluster",
                                     plot_type = plot_type,
                                     cor_cutoff = cor_cutoff)
  
  
  seurat_data <- cluster_res$object
  
  seurat_res_krentz <- cluster_res$RNA
  
  plotDimRed(seurat_data, col_by = "RNA_krentz_celltype",
             plot_type = plot_type)
  
  write.csv(seurat_res_krentz,
            file.path(save_dir, "files/celltype_mapping_krentz_2018_2020.csv"))
  
  
  
  #-------------------------------------------------------------------------------
  
  ############
  # Combined #
  ############
  
  colnames(seurat_res_baron) <- paste0(colnames(seurat_res_baron), "_m")
  
  colnames(seurat_res_baron_human) <- paste0(colnames(seurat_res_baron_human), "_h")
  
  full_res <- cbind(cbind(seurat_res_baron, seurat_res_byrnes),
                    seurat_res_tabula_muris, seurat_res_baron_human,
                    seurat_res_qadir, seurat_res_muraro,
                    seurat_res_krentz)
  
  seurat_cluster <- cor_to_call(full_res) %>% 
    mutate(type = ifelse(r < cor_cutoff, "undetermined", type))
  
  new_clusters <- seurat_cluster$type
  
  names(new_clusters) <- seurat_cluster$cluster
  
  seurat_data$RNA_combined_celltype <- new_clusters[seurat_data$RNA_cluster]
  
  plotDimRed(seurat_data, col_by = "RNA_combined_celltype",
             plot_type = plot_type)
  
  seurat_data$RNA_combined_celltype <- gsub("_h$", "",
                                            seurat_data$RNA_combined_celltype)
  seurat_data$RNA_combined_celltype <- gsub("_m$", "",
                                            seurat_data$RNA_combined_celltype)
  seurat_data$RNA_combined_celltype <- gsub("Acinar", "acinar",
                                            seurat_data$RNA_combined_celltype)
  
  seurat_data$RNA_combined_celltype <- gsub("Ductal", "ductal",
                                            seurat_data$RNA_combined_celltype)
  
  
  plotDimRed(seurat_data, col_by = "RNA_combined_celltype",
             plot_type = plot_type)
  
  #-------------------------------------------------------------------------------
  ############################
  # Combined with ref winner #
  ############################
  
  colnames(seurat_res_baron) <- paste0(colnames(seurat_res_baron), "_baron")
  
  colnames(seurat_res_baron_human) <- paste0(colnames(seurat_res_baron_human),
                                             "_baron")
  colnames(seurat_res_byrnes) <- paste0(colnames(seurat_res_byrnes),
                                        "_m_byrnes")
  
  colnames(seurat_res_tabula_muris) <- paste0(colnames(seurat_res_tabula_muris),
                                              "_m_tabula_muris")
  
  colnames(seurat_res_qadir) <- paste0(colnames(seurat_res_qadir),
                                       "_h_qadir")
  
  colnames(seurat_res_muraro) <- paste0(colnames(seurat_res_muraro),
                                        "_h_muraro")
  
  colnames(seurat_res_krentz) <- paste0(colnames(seurat_res_krentz),
                                        "_h_krentz")
  
  full_res <- cbind(cbind(seurat_res_baron, seurat_res_byrnes),
                    seurat_res_tabula_muris, seurat_res_baron_human,
                    seurat_res_qadir, seurat_res_muraro,
                    seurat_res_krentz)
  
  seurat_cluster <- cor_to_call(full_res) %>% 
    mutate(type = ifelse(r < cor_cutoff, "undetermined", type))
  
  new_clusters <- seurat_cluster$type
  
  names(new_clusters) <- seurat_cluster$cluster
  
  seurat_data$RNA_combined_celltype_ref_name <- new_clusters[seurat_data$RNA_cluster]
  
  plotDimRed(seurat_data, col_by = "RNA_combined_celltype_ref_name",
             plot_type = plot_type)
  
  
  write.csv(full_res, file = file.path(save_dir, "files", "all_clustifyr_cors.csv"))
  
}

saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_processed.rds"))