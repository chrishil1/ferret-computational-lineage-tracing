library(Seurat)
library(here)
library(ggplot2)
library(scAnalysisR)

source(here("src/common_setup.R"))

#' Generate Gene Expression Dotplots by Cell Type and Developmental Stage
#'
#' Creates and saves Seurat dotplots visualizing the expression of specified
#' genes across cell types and developmental stages. Grouping variables are
#' determined by species, with ferret data grouped by cell type and day,
#' and human data grouped by cell type and post-conception week (PCW).
#'
#' Parameters
#' ----------
#' object : Seurat
#'     Seurat object containing gene expression data and metadata.
#' species : character
#'     Species associated with the Seurat object. Use "ferret" for ferret
#'     data; all other values are treated as human data.
#' genes : character vector
#'     Gene symbols to include in the dotplots.
#' sample : character
#'     Sample identifier used to label output files (e.g., "FS" or "WA").
#' width : numeric, optional
#'     Width of the output PDF files in inches. Default is 10.
plot_dotplots <- function(object,
                          species,
                          genes,
                          sample,
                          width = 10) {
  if (species == "ferret") {
    pdf(file.path(outdir, paste0(species, "_", sample, "_dotplot_by_celltype.pdf")),
      width = width, height = 6
    )
    print(DotPlot(object,
      features = genes, group.by = "final_renamed_celltype",
      cols = c("blue", "red")
    ) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_y_discrete(limits = rev(levels(factor(object$final_renamed_celltype)))) +
      ggtitle("Gene Expression by Celltype"))
    dev.off()

    pdf(file.path(outdir, paste0(species, "_", sample, "_dotplot_by_day.pdf")),
      width = width, height = 6
    )
    print(DotPlot(object,
      features = genes, group.by = "day",
      cols = c("blue", "red")
    ) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_y_discrete(limits = rev(levels(factor(object$day)))) +
      ggtitle("Gene Expression by Day"))
    dev.off()
  } else {
    pdf(file.path(outdir, paste0(species, "_", sample, "_dotplot_by_celltype.pdf")),
      width = width, height = 6
    )
    print(DotPlot(object,
      features = genes, group.by = "Celltype",
      cols = c("blue", "red")
    ) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_y_discrete(limits = rev(levels(factor(object$Celltype)))) +
      ggtitle("Gene Expression by Celltype"))
    dev.off()

    pdf(file.path(outdir, paste0(species, "_", sample, "_dotplot_by_PCW.pdf")),
      width = width, height = 6
    )
    print(DotPlot(object,
      features = genes, group.by = "PCW",
      cols = c("blue", "red")
    ) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_y_discrete(limits = rev(levels(factor(object$PCW)))) +
      ggtitle("Gene Expression by PCW"))
    dev.off()
  }
}


#' Generate Gene Expression UMAPs and Violin Plots
#'
#' Creates and saves Seurat feature plots and violin plots visualizing
#' the expression of specified genes across cells, cell types, and
#' developmental stages. For each gene, a UMAP feature plot and two
#' violin plots are generated. Grouping variables are determined by
#' species, with ferret data grouped by cell type and day, and human
#' data grouped by cell type and post-conception week (PCW).
#'
#' Parameters
#' ----------
#' object : Seurat
#'     Seurat object containing gene expression data, UMAP embeddings,
#'     and associated metadata.
#' species : character
#'     Species associated with the Seurat object. Use "ferret" for ferret
#'     data; all other values are treated as human data.
#' genes : character vector
#'     Gene symbols to visualize in UMAPs and violin plots.
#' sample : character
#'     Sample identifier used to label output files (e.g., "FS" or "WA").
plot_gene_set <- function(object,
                          species,
                          genes,
                          sample) {
  for (gene in genes) {
    pdf(file.path(outdir, paste0(species, "_", sample, "_", gene, "_umap.pdf")),
      width = 8, height = 8
    )
    print(FeaturePlot(object, features = gene, reduction = "umap"))
    dev.off()

    if (species == "ferret") {
      pdf(file.path(outdir, paste0(species, "_", sample, "_", gene, "_vln_day.pdf")),
        width = 8, height = 6
      )
      print(VlnPlot(object,
        cols = c("#4285F4", "#34A853", "#FBBC05", "#EA4335"),
        features = gene, group.by = "day", pt.size = 0
      ))
      dev.off()

      pdf(file.path(outdir, paste0(species, "_", sample, "_", gene, "_vln_celltype.pdf")),
        width = 10, height = 6
      )
      print(VlnPlot(object,
        cols = human_palette, features = gene,
        group.by = "final_renamed_celltype", pt.size = 0
      )
      + NoLegend())
      dev.off()
    } else {
      pdf(file.path(outdir, paste0(species, "_", sample, "_", gene, "_vln_PCW.pdf")),
        width = 8, height = 6
      )
      print(VlnPlot(object,
        cols = c("#4285F4", "#34A853", "#FBBC05"),
        features = gene, group.by = "PCW", pt.size = 0
      ))
      dev.off()

      pdf(file.path(outdir, paste0(species, "_", sample, "_", gene, "_vln_celltype.pdf")),
        width = 10, height = 6
      )
      print(VlnPlot(object,
        cols = human_palette, features = gene,
        group.by = "Celltype", pt.size = 0
      ) + NoLegend())
      dev.off()
    }
  }
}

# Create output directory
outdir <- "figures/15_GOI_plots/"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Load in objects of interest
ferret_object <- readRDS(file.path(data_dir, "05_final_ferret.rds"))
human_object <- readRDS(file.path(data_dir, "05_final_human.rds"))

objects <- list(ferret = ferret_object, human = human_object)

# Define gene lists for the four categories
gene_sets <- list(
  FS = list(
    acinar = c(
      "BARX2", "BHLHA15", "CELA1", "CPA4", "EGFR", "ELF3", "ERN1",
      "IL6ST", "KLF15", "MYC", "PDIA2", "PRSS3", "RET", "ZBTB16"
    ),
    ductal = c(
      "CA2", "EGFR", "ELF3", "HHEX", "ITGB4", "KCNK5", "NFIB",
      "NFKB1", "VEGFA", "ZBTB4"
    )
  ),
  WA = list(
    acinar = c(
      "AQP8", "CCKAR", "CEL", "CELA3A", "CPA1", "CPA2", "CUZD1",
      "ERP27", "GATA4", "GP2", "HNF4G", "NFKBIZ", "PLA2G1B",
      "PNLIPRP1", "RBPJL", "XBP1"
    ),
    ductal = c(
      "ANXA4", "BICC1", "CFTR", "DPP4", "EHF", "FGFR3", "GLIS3",
      "HNF1B", "LAMA1", "LAMB1", "SCTR", "SLC4A4", "TNFRSF19",
      "WNT5B"
    )
  )
)

# Iterate through all combinations and generate plots
for (species in names(objects)) {
  for (sample in names(gene_sets)) {
    for (celltype in names(gene_sets[[sample]])) {
      plot_gene_set(
        objects[[species]],
        species,
        gene_sets[[sample]][[celltype]],
        sample
      )
      plot_dotplots(
        objects[[species]],
        species,
        gene_sets[[sample]][[celltype]],
        sample
      )
    }
  }
}
