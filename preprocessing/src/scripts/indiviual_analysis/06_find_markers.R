library(Seurat)
library(tidyverse)
library(cowplot)
library(openxlsx)
library(here)
library(scAnalysisR)

source("src/scripts/functions.R")

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

args <- commandArgs(trailingOnly = TRUE)

sample <- args[[1]]
sample <- gsub("__.*", "", sample)
#sample <- "combined_group_2"
#sample <- "combined_group_1"

sample_info <- args[[3]]
#sample_info <- here("files/sample_info.tsv")

results_dir <- args[[2]]
#results_dir <- here("results")

sample_info <- read.table(sample_info, fill = TRUE, header = TRUE)

sample_info <- sample_info[sample_info$sample == sample,]

HTO <- sample_info$HTO
ADT <- sample_info$ADT
hash_ident <- sample_info$hash_ident

ADT_pca <- sample_info$adt_pca
run_adt_umap <- sample_info$adt_umap

RNA_pcs <- sample_info$PCs
resolution <- sample_info$resolution


cell_types <- "RNA_combined_celltype"
clusters <- "RNA_cluster"
pval <- 0.05
logfc <- 0.5

HTO <- FALSE
ADT <- FALSE

if(normalization_method == "SCT"){
  SCT <- TRUE
  seurat_assay <- "SCT"
} else {
  SCT <- FALSE
  seurat_assay <- "RNA"
}

# Set directories
save_dir <- file.path(results_dir, "R_analysis", sample)

# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))

mapping_file <- read.csv(here("files/species_mapping_file.csv"))

gene_lists <- read.table(file.path("/beevol", "home", "wellskri", "Analysis", 
                                   "references", "single_cell_references",
                                   "gene_lists", 
                                   "PanglaoDB_markers_27_Mar_2020.tsv.gz"),
                         sep = "\t", header = TRUE)

# Cell type DE -----------------------------------------------------------------


# RNA cluster DE ---------------------------------------------------------------

seurat_data$cluster_celltype <- paste0(seurat_data[[clusters]][[1]], "_",
                                       seurat_data[[cell_types]][[1]])

marker_list <- find_write_markers_orthologs(seurat_object = seurat_data,
                                            meta_col = "cluster_celltype",
                                            pval = pval,
                                            logfc = logfc,
                                            assay = "RNA",
                                            save_dir = save_dir,
                                            mapping_file = mapping_file,
                                            mapping_gene_col = "gene_id",
                                            mapping_ortholog_col = c("Mouse.gene.name",
                                                                     "Human.gene.name",
                                                                     "Dog.gene.name",
                                                                     "Pig.gene.name"),
                                            gene_lists = NULL)

write.csv(marker_list, file.path(save_dir, "files", "all_markers.csv"))


gene_list <- c("INS", "SST", "GCG", "MKI67", "CFTR")

plots <- featDistPlot(seurat_data, geneset = gene_list, sep_by = "RNA_cluster",
                      combine = FALSE)

pdf(file.path(save_dir, "images", "marker_gene_violin.pdf"),
    height = 6, width = 15)
print(plots)
dev.off()


library(enrichR)
websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes   
}

if (websiteLive) dbs <- listEnrichrDbs()

dbs_use <- c("Human_Gene_Atlas", "Mouse_Gene_Atlas",
             "PanglaoDB_Augmented_2021", "Tabula_Muris",
             "Tabula_Sapiens")

wb <- openxlsx::createWorkbook()

human <- lapply(unique(marker_list$cluster), function(x){
  markers_use <- marker_list %>%
    dplyr::filter(cluster == x, p_val_adj < 0.05)
  
  enriched <- enrichr(unique(markers_use$Human.gene.name), dbs_use)
  
  top_enriched <- lapply(enriched, function(enriched_df){
    enriched_df <- enriched_df %>%
      dplyr::arrange(desc(Combined.Score))
    
    top_term <- enriched_df[1,]
  })
  top_enriched <- do.call(rbind, top_enriched)
  
  top_enriched$database <- dbs_use

  x <- as.character(x)
  # Save as excel file
  openxlsx::addWorksheet(wb = wb, sheetName = x)
  openxlsx::writeData(wb = wb, sheet = x, x = top_enriched)
  
  return(top_enriched)
})

openxlsx::saveWorkbook(wb = wb, file = file.path(save_dir, "files",
                                                 "enrichr_human_gene.xlsx"),
                       overwrite = TRUE)

# Try with module score --------------------------------------------------------
celltypes_use <- c("Acinar cells", "Alpha cells", "B cells",
                   "Beta cells", "Chromaffin cells", "Delta cells",
                   "Ductal cells", "Endothelial cells", "Gamma (PP) cells",
                   "Macrophages", "Monocytes")

all_plots <- lapply(unique(gene_lists$cell.type), function(x){
  celltype_genes <- gene_lists %>%
    dplyr::filter(cell.type == x)
  
  all_genes <- celltype_genes$official.gene.symbol
  
  mapping_genes <- mapping_file %>%
    dplyr::filter(Human.gene.name %in% all_genes) %>%
    dplyr::filter(!is.na(Human.gene.name)) %>%
    dplyr::filter(Human.gene.name != "")
  
  human_genes <- unique(mapping_genes$gene_id)
  
  ferret_genes <- setdiff(all_genes, unique(mapping_genes$Human.gene.name))
  
  x <- make.names(x)
  
  so <- AddModuleScore(seurat_data,
                       features = list(c(human_genes, ferret_genes)),
                                       name = x)
  
  return_plot <- featDistPlot(so, geneset = paste0(x, "1"),
                              sep_by = "cluster_celltype",
                              combine = FALSE, col_by = "cluster_celltype")
  
  return(return_plot)
})

pdf(file.path(save_dir, "images", "panglao_db_module_violin.pdf"),
    height = 6, width = 15)
print(all_plots)
dev.off()

marker_genes <- c("INS", # Beta
                  "GCG", # Alpha
                  "SST", # Delta
                  "NKX2-2", "NKX6-1", "PAX6", # Endocrine
                  "CFTR", # ductal
                  "PDGFRA", # fibroblast
                  "CD3D", # T cells
                  "ITGAM", # Myeloid
                  "SOX2", "NOTCH1", # Neurons
                  "ALB", # Hepatocyte
                  "CD19", # B cell
                  "IL2RA", # T cell
                  "NES", # Stellate cell
                  "CDH5", "PECAM1", # Endothelail
                  "OLFM4" # epithelail
                  )

marker_plot <- featDistPlot(seurat_data, geneset = marker_genes,
                            sep_by = "cluster_celltype",
                            combine = FALSE, col_by = "cluster_celltype")

pdf(file.path(save_dir, "images", "RNA_markers_violin.pdf"),
    height = 6, width = 15)
print(marker_plot)
dev.off()



if(sample == "combined_group_1"){
  rename_celltype <- c("0" = "fibroblast", # Updated
                       "1" = "fibroblast", # Updated
                       "2" = "fibroblast", # Updated
                       "3" = "fibroblast", # Updated
                       "4" = "ductal", # Updated
                       "5" = "acinar", # Updated
                       "6" = "beta", # Updated
                       "7" = "fibroblast", # Updated
                       "8" = "stellate", # Updated
                       "9" = "endothelial", # Updated
                       "10" = "ductal", # Updated
                       "11" = "neurons", # Updated
                       "12" = "fibroblast", # Updated
                       "13" = "fibroblast",# Updated
                       "14" = "immune_myeloid", # Updated
                       "15" = "fibroblast", # Updated
                       "16" = "immune_myeloid", # Updated
                       "17" = "fev_endocrine_progenitor", # Updated
                       "18" = "alpha", # Updated
                       "19" = "delta", # Updated
                       "20" = "ductal", # Updated 
                       "21" = "immune_t", # Updated
                       "22" = "endothelial", # Updated
                       "23" = "immune_b", # Updated
                       "24" = "hepatocyte") # Updated

} else if(sample == "combined_group_2"){
  rename_celltype <- c("0" = "fibroblast", # Updated
                       "1" = "ductal", # Updated 
                       "2" = "fibroblast", # Updated
                       "3" = "fibroblast", # Updated
                       "4" = "fibroblast", # Updated
                       "5" = "acinar", # Updated
                       "6" = "endothelial", # Updated
                       "7" = "ductal", # Updated
                       "8" = "stellate", # Updated
                       "9" = "fibroblast", # Updated
                       "10" = "fibroblast", # Updated
                       "11" = "beta", # Updated
                       "12" = "fibroblast", # Updated
                       "13" = "kidney", # Updated 
                       "14" = "fibroblast", # Updated
                       "15" = "neurons", # Updated
                       "16" = "immune_myeloid", # Updated
                       "17" = "acinar", # Updated
                       "18" = "fibroblast", # Updated
                       "19" = "fibroblast", # Updated
                       "20" = "epithelial", # Updated
                       "21" = "fev_endocrine_progenitor", # Updated
                       "22" = "kidney", # Updated
                       "23" = "fibroblast", # Updated
                       "24" = "immune_myeloid", # Updated
                       "25" = "stellate", # Updated
                       "26" = "ductal", # Updated
                       "27" = "alpha", # Updated
                       "28" = "delta", # Updated
                       "29" = "neurons", # Updated
                       "30" = "ductal", # Updated 
                       "31" = "immune_t", # Updated
                       "32" = "immune_t", # Updated
                       "33" = "fibroblast", # Updated
                       "34" = "fibroblast", # Updated
                       "35" = "immune_b", # Updated
                       "36" = "endothelial", # Updated
                       "37" = "hepatocyte", # Updated
                       "38" = "ductal") # Updated  
}

seurat_data$final_celltype <- rename_celltype[seurat_data$RNA_cluster]
# 
# seurat_data$final_celltype <- ifelse(is.na(seurat_data$final_celltype),
#                                      seurat_data$RNA_combined_celltype,
#                            seurat_data$final_celltype)
# 
# seurat_data$final_celltype <- ifelse(seurat_data$final_celltype == "pancreatic_D_cell",
#                                      "fev_endocrine_progenitor",
#                                      seurat_data$final_celltype)

seurat_data$cluster_celltype <- paste(seurat_data$RNA_cluster,
                                      seurat_data$final_celltype,
                                      sep = "_")

pdf(file.path(save_dir, "images", "final_celltype_umap.pdf"))
print(plotDimRed(seurat_data, col_by = "final_celltype", plot_type = "rna.umap"))
dev.off()

all_plots <- lapply(unique(seurat_data$final_celltype), function(x){
  plotDimRed(seurat_data, col_by = "RNA_cluster", plot_type = "rna.umap",
             highlight_group = TRUE, group = x, meta_data_col = "final_celltype",
             ggrastr = TRUE)
})


pdf(file.path(save_dir, "images", "UMAP_celltype_single.pdf"),
    height = 6, width = 8)
print(all_plots)
dev.off()


marker_plot <- featDistPlot(seurat_data, geneset = marker_genes,
                            sep_by = "cluster_celltype",
                            combine = FALSE, col_by = "cluster_celltype")

pdf(file.path(save_dir, "images", "RNA_markers_violin_final.pdf"),
    height = 6, width = 15)
print(marker_plot)
dev.off()


saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_processed.rds"))

# Some notes from enrichR combined_1
# Cluster 8 - mesothelial cells in liver (descartes), pancreas fibroblast (tabula sapiens), stem cell?
# Cluster 0 - stromal cells in adrenal (descartes), pancreas fibroblast (tabula sapiens)
# Cluster 2 - stellate cells (descartes), gamma pp cells (panglao db)
# Cluster 13 - muscle almost all say smooth muscle
# Cluster 15 - no solid hits, maybe fibroblast
# Cluster 14 - macrophage
# Cluster 11 - neuron
# Cluster 21 - Nothing strong. Ductal by cftr expression
# Cluster 24 - Liver
# Cluster 9 - Endothelial
# Cluster 10 - Acinar
# Cluster 12 - Pancreas fibroblast
# Cluster 5 - Beta cells
# Cluster 3 - Ductal cell
# Cluster 20 - Stem cell (many different tissues, but all stem cell)
# Cluster 4 - Fibroblast (keep as stellate?)
# Cluster 18 - Alpha cell
# Cluster 16 - Immune cell (T or NK)
# Cluster 7 - Fibroblast (keep as stelalte?)
# Cluster 22 - Endothelial (keep)
# Cluster 1 - Fibroblast (keep as stelate?)
# Cluster 19 - Delta cell (keep)
# Cluster 18 - Acinar cells
# Cluster 6 - Acinar cells
# Cluster 23 - Immune

# Based on CFTR expression, CLusters 3 and 21 are ductal - this agrees with the cluster naming
# Based on PRSS1 (LOC101683959) and KLK1 (LOC101690899) no obvious acinar
# Based on GCG, cluster 18 is alpha
# Based on INS, there is no obvious Beta
# Based on SST, cluster 19 is delta
# Clusters 17, 18, 19, 5 have some expression of Nkx2-2
# No OCT4 (POU5F1) detected