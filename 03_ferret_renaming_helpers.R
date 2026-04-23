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
