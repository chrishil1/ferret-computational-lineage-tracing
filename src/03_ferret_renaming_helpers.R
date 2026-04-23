# Define color palette for use with umap/barplots
human_palette <- c(
  "Acinar Immature" = "#e64463",
  "Acinar Maturing" = "#62be4a",
  "Acinar Progenitor" = "#9255cb",
  "Acinar Proliferating" = "#b2b336",
  "Acinar Secretory" = "#5069d7",
  "Alpha" = "#dc8f2f",
  "B lymphocyte" = "#d579e1",
  "Beta" = "#4c8a35",
  "Delta" = "#ca42a9",
  "Dendritic Cell" = "#54c081",
  "Ductal ECM" = "#d83f86",
  "Ductal Immature" = "#52bfad",
  "Ductal Maturing" = "#cc4724",
  "Ductal Proliferating" = "#4ab4dd",
  "Ducto-Acinar" = "#dc6c4d",
  "Endothelial" = "#5982c5",
  "Endothelial Lymphatic" = "#cca15d",
  "Erythroblast" = "#9384d9",
  "Mesenchyme Proliferating" = "#727229",
  "Mesenchyme Subtype 1" = "#ae57a1",
  "Mesenchyme Subtype 2" = "#9cb26a",
  "Mesenchyme Subtype 3" = "#755298",
  "Mesenchyme Subtype 4" = "#37835d",
  "Mono-Macrophage" = "#b83849",
  "Neuron" = "#ca96d2",
  "Pancreatic Progenitor" = "#9a5e2c",
  "Phagocytic Macrophage" = "#e47da1",
  "Schwann" = "#e5907a",
  "Stromal" = "#9f4970",
  "T lymphocyte" = "#aa5659"
)

blue_yellow <- c(
  "#352A86", "#343DAE", "#0262E0", "#1389D2",
  "#2DB7A3", "#A5BE6A", "#F8BA43", "#F6DA23", "#F8FA0D"
)

#' Create a Stacked Percentage Bar Plot
#'
#' Generates a stacked bar plot showing the percentage breakdown of a grouping
#' variable across categories of a given attribute.
#'
#' @param column Character. The name of the metadata column in seurat_obj to use
#'   as the stacking/fill variable (e.g. sample, day, batch). Will be coerced
#'   to a sorted factor.
#' @param seurat_obj A Seurat object containing cell metadata accessible via
#'   meta.data.
#' @param attribute Character. The name of the metadata column to place on the
#'   x-axis (e.g. cell type, cluster). Each unique value becomes one bar.
#' @param palette Named character vector. Colours to use for the fill variable,
#'   passed to scale_fill_manual(). Names should match the levels of
#'   column.
#'
#' @return A ggplot object showing a stacked bar chart where each bar
#'   sums to 100%, with attribute on the x-axis and column as
#'   the fill variable.
#'
#' @examples
#' create_barplot("sample", seurat_object, "celltype", sample_palette)
create_barplot <- function(column, seurat_obj, attribute, palette) {
  # Order the column variable alphabetically or numerically
  unique_samples <- sort(unique(seurat_obj[[column, drop = TRUE]]))
  seurat_obj[[column]] <- factor(seurat_obj[[column, drop = TRUE]],
    levels = unique_samples
  )

  # Create a table of counts by day and the column variable
  data <- seurat_obj@meta.data
  summary_table <- table(data[[attribute]], data[[column]])
  summary_df <- as.data.frame(summary_table)
  summary_df <- summary_df %>%
    group_by(Var1) %>%
    mutate(Percent = Freq / sum(Freq) * 100) %>%
    ungroup()

  # Create the stacked bar plot, with each bar summing to 100%
  figure <- ggplot(summary_df, aes(x = Var1, y = Percent, fill = Var2)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = palette, name = column) +
    xlab("PCW") +
    ylab("Percentage") +
    theme_minimal()

  return(figure)
}

#' Generate a Differential Expression Heatmap
#'
#' Runs FindAllMarkers() on a Seurat object, selects the top differentially
#' expressed genes per cluster, and visualises average scaled expression as a
#' heatmap using pheatmap. Marker results are also written to a CSV file.
#'
#' @param seurat_object A Seurat object to perform differential expression on.
#' @param group_by Character. Metadata column used to group cells for average
#'   expression calculation and for naming the output CSV. If set to
#'   "final_combined_celltype", the label "fev_endocrine_progenitor" is
#'   shortened to "fev" before processing.
#' @param n Integer. Number of top DE genes to select per cluster. Default 10.
#' @param custom_order Character vector. Column order for the heatmap, should
#'   match the group levels in group_by. If NULL, columns are not reordered.
#' @param upregulated Logical. If TRUE, selects the top n genes ranked by
#'   avg_log2FC (upregulated only). If FALSE, selects by abs(avg_log2FC)
#'   (up and downregulated). Default FALSE.
#' @param palette Character vector. Colour palette passed to pheatmap.
#'   If NULL, the pheatmap default is used.
#'
#' @return A pheatmap object showing scaled average expression of the top DE
#'   genes, with rows clustered and columns ordered by custom_order. As a side
#'   effect, writes all_markers_<group_by>.csv to the working directory.
#'
#' @examples
#' get_de_heatmap(
#'   seurat_object = seurat_obj,
#'   group_by = "final_combined_celltype",
#'   n = 10,
#'   custom_order = c("Alpha", "Beta", "Delta", "Ductal"),
#'   upregulated = TRUE,
#'   palette = viridis(100),
#'   ident = "final_combined_celltype"
#' )
get_de_heatmap <- function(seurat_object,
                           group_by,
                           n = 10,
                           custom_order = NULL,
                           upregulated = FALSE,
                           palette = NULL) {
  # Rename for better readability
  if (group_by == "final_combined_celltype") {
    seurat_object$final_combined_celltype <- gsub(
      "fev_endocrine_progenitor",
      "fev",
      seurat_object$final_combined_celltype
    )
  }

  # Set cluster identities for the Seurat object
  markers_file <- paste0("all_markers_", group_by, ".csv")
  if (file.exists(markers_file)) {
    all_markers <- read.csv(markers_file)
  } else {
    all_markers <- FindAllMarkers(seurat_object)
    write.csv(all_markers, file = markers_file, row.names = FALSE)
  }

  if (upregulated) {
    # Filter top upregulated genes based on log2FC and padj
    top_genes_df <- all_markers %>%
      filter(abs(avg_log2FC) > 0.5, p_val_adj < 0.05) %>%
      group_by(cluster) %>%
      arrange(cluster, desc(avg_log2FC)) %>%
      slice_head(n = n)

    title <- paste0("Top ", n, " Upregulated DE Genes per cluster")
  } else {
    # Filter top genes based on log2FC and padj
    top_genes_df <- all_markers %>%
      filter(abs(avg_log2FC) > 0.5, p_val_adj < 0.05) %>%
      group_by(cluster) %>%
      arrange(cluster, desc(abs(avg_log2FC))) %>%
      slice_head(n = n)

    title <- paste0("Top ", n, " DE Genes per cluster")
  }

  # Calculate the average expression for each gene per cell type
  seurat_avg <- AverageExpression(seurat_object,
    group.by = group_by
  )

  # Subset out by just looking at the top genes extracted from n
  expression_data <- seurat_avg$RNA[unique(top_genes_df$gene), ]

  # Scale and reorder expression data for heatmap
  scaled_expression_data <- t(scale(t(expression_data)))
  reordered_expression_data <- scaled_expression_data[,
    custom_order,
    drop = FALSE
  ]

  # Create an output heatmap summarizing this data
  heatmap <- pheatmap(reordered_expression_data,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_rownames = FALSE,
    show_colnames = TRUE,
    main = title,
    color = palette
  )

  return(heatmap)
}
