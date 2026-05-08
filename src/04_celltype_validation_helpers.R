#' Generate a Validation Heatmap of Gene Expression Across Cell Types
#'
#' Computes average gene expression per cell type from a Seurat object and
#' renders a scaled heatmap via pheatmap. Markers can either be discovered
#' automatically with FindAllMarkers or supplied as a pre-defined gene list.
#' Results are optionally reordered by a custom cell-type ordering.
#'
#' @param seurat_object A Seurat object containing the single-cell data.
#' @param n Integer. Number of top differentially expressed genes to select
#'   per cluster when markers are computed automatically. Default: 10.
#' @param custom_order Character vector. Optional custom ordering of cell-type
#'   column names for the heatmap. If NULL, columns appear in their default
#'   order. Default: NULL.
#' @param upregulated Logical. If TRUE, selects only upregulated genes
#'   (avg_log2FC > 0.5); if FALSE, selects genes by absolute log2FC.
#'   Only applies when genes_of_interest is NULL. Default: FALSE.
#' @param genes_of_interest Character vector. Optional pre-defined list of
#'   genes to display. If provided, skips marker discovery entirely.
#'   Default: NULL.
#' @param cluster_rows Logical. Whether to hierarchically cluster heatmap
#'   rows. Passed directly to pheatmap. Default: TRUE.
#' @param idents Character. Name of the Seurat metadata column to use as
#'   cell identity/grouping variable. Default: "RNA_baron_celltype".
#' @param assay Character. Seurat assay to pull expression data from.
#'   Either "RNA" or "SCT". Default: "RNA".
#' @param show_rownames Logical. Whether to display gene names on the heatmap
#'   y-axis. Default: TRUE.
#'
#' @return A named list with two elements:
#'     heatmap: A pheatmap object of the scaled expression heatmap.
#'     expression_data: A numeric matrix of scaled, reordered average
#'       expression values (genes x cell types).
#'
#' @examples
#' result <- get_validation_heatmap(
#'   seurat_object   = my_seurat,
#'   n               = 10,
#'   custom_order    = c("Alpha", "Beta", "Delta"),
#'   upregulated     = TRUE,
#'   idents          = "RNA_baron_celltype",
#'   assay           = "RNA"
#' )
#' result$heatmap
#' head(result$expression_data)
get_validation_heatmap <- function(seurat_object,
                                   n = 10,
                                   custom_order = NULL,
                                   upregulated = FALSE,
                                   genes_of_interest = NULL,
                                   cluster_rows = TRUE,
                                   idents = "RNA_baron_celltype",
                                   assay = "RNA",
                                   show_rownames = TRUE) {
  # Set cluster identities for the Seurat object
  Idents(seurat_object) <- idents

  if (is.null(genes_of_interest)) {
    # Find all markers if no specific gene list is provided
    all_markers <- FindAllMarkers(seurat_object, assay = assay)

    write.csv(
      all_markers,
      file = paste0("all_markers_", idents, ".csv"),
      row.names = TRUE
    )

    if (upregulated) {
      top_genes_csv <- all_markers %>% filter(
        avg_log2FC > 0.5,
        p_val_adj < 0.05
      )
      write.csv(top_genes_csv,
        file = here(
          image_dir,
          "all_markers_sig_up.csv"
        ),
        row.names = FALSE
      )

      # Filter top upregulated genes based on log2FC and padj
      top_genes_df <- all_markers %>%
        filter(abs(avg_log2FC) > 0.5, p_val_adj < 0.05) %>%
        group_by(cluster) %>%
        arrange(cluster, desc(avg_log2FC)) %>%
        slice_head(n = n)

      title <- paste0("Top ", n, " Upregulated DE Genes per cluster")
    } else {
      top_genes_csv <- all_markers %>% filter(
        abs(avg_log2FC) > 0.5,
        p_val_adj < 0.05
      )
      write.csv(top_genes_csv,
        file = here("all_markers_sig.csv"),
        row.names = FALSE
      )

      # Filter top genes based on log2FC and padj
      top_genes_df <- all_markers %>%
        filter(abs(avg_log2FC) > 0.5, p_val_adj < 0.05) %>%
        group_by(cluster) %>%
        arrange(cluster, desc(abs(avg_log2FC))) %>%
        slice_head(n = n)

      title <- paste0("Top ", n, " DE Genes per cluster")
    }

    selected_genes <- unique(top_genes_df$gene)
  } else {
    # Otherwise, just use that gene list
    selected_genes <- genes_of_interest
    title <- "Normalized Expression by Celltype - RNA"
  }

  # Calculate the average expression for each gene per cell type
  seurat_avg <- AverageExpression(seurat_object, group.by = idents)

  # Subset expression data for selected genes
  if (assay == "RNA") {
    expression_data <- seurat_avg$RNA[selected_genes, , drop = FALSE]
  } else {
    expression_data <- seurat_avg$SCT[selected_genes, , drop = FALSE]
  }

  # Prevent overflow by setting to some maximum value
  expression_data[expression_data > 1e150] <- 1e150

  # Scale and reorder expression data for heatmap
  scaled_expression_data <- t(scale(t(expression_data)))

  # If a custom order exists, reorder
  if (!is.null(custom_order)) {
    reordered_expression_data <-
      scaled_expression_data[, custom_order, drop = FALSE]
  } else {
    reordered_expression_data <- scaled_expression_data
  }

  # Construct heatmap of expression data
  heatmap <- pheatmap(reordered_expression_data,
    cluster_rows = cluster_rows,
    cluster_cols = FALSE,
    show_rownames = show_rownames,
    show_colnames = TRUE,
    main = title,
    fontsize_row = 5,
    color = blue_yellow
  )

  return(list(heatmap = heatmap, expression_data = reordered_expression_data))
}
