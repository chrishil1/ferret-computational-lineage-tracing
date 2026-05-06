library(clusterProfiler)
library(gprofiler2)
library(here)
library(org.Hs.eg.db)
library(pathview)
library(Seurat)
library(tidyverse)

#' Run KEGG pathway enrichment and generate bubble plots.
#'
#' This function performs gene set enrichment analysis using `gost`
#' with a custom background, extracts KEGG pathway results, and
#' generates two bubble plots summarizing the top enriched pathways.
#'
#' @param query_names Character vector of gene identifiers to test for
#'        enrichment.
#' @param background Character vector of background gene identifiers.
#' @param name Character string used for plot titles and output filenames.
#' @param organism Character string specifying the organism
#'        (default: "hsapiens").
#'
#' @return NULL
#'   The function saves two PDF plots to disk and does not return a value.
#'
#' @examples
#' query_genes <- c("INS", "GCG", "SST", "PDX1", "MAFA")
#' background_genes <- c("INS", "GCG", "SST", "PDX1", "MAFA", "NKX6-1")
#' run_enrichment(
#'   query_names = query_genes,
#'   background = background_genes,
#'   name = "example_celltype"
#' )
run_enrichment <- function(query_names,
                           background,
                           name,
                           organism = "hsapiens") {
  # Run gost with the query genes and custom background, ordered by default
  enrichment_results <- gost(
    query = query_names,
    custom_bg = background,
    organism = organism,
    ordered_query = TRUE,
    significant = FALSE
  )

  # Extract the KEGG pathways to run with pathview
  kegg_terms <-
    enrichment_results$result[grep(
      "^KEGG:",
      enrichment_results$result$term_id
    ), ]
  kegg_ids <- sub("^KEGG:", "", kegg_terms$term_id)
  kegg_data <- data.frame(
    term_id = kegg_ids,
    term_name = kegg_terms$term_name,
    p_value = kegg_terms$p_value,
    recall = kegg_terms$recall,
    intersection_size = kegg_terms$intersection_size
  )

  # Sort the dataframe by p-value
  kegg_data_sorted <- kegg_data[order(kegg_data$p_value), ]

  # Select the top 25 entries
  kegg_data_top25 <- head(kegg_data_sorted, 25)

  # Create a bubble plot with recall as the x axis
  bubble_plot <- ggplot(kegg_data_top25, aes(
    x = recall,
    y = reorder(term_name, recall)
  )) +
    geom_point(aes(
      size = intersection_size,
      color = -log10(p_value)
    ), alpha = 0.7) +
    scale_size_continuous(range = c(3, 10)) +
    scale_color_gradient(low = "blue", high = "red") +
    labs(
      x = "Recall",
      y = "KEGG Term",
      size = "intersection_size",
      color = "-log10(p_value)",
      title = paste0("Top 25 KEGG Enrichment Pathways, ", name)
    ) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10))

  # Create a bubble plot with -log10(p-value) as the x axis
  bubble_plot_p_value <- ggplot(
    kegg_data_top25,
    aes(
      x = -log10(p_value),
      y = reorder(term_name, -log10(p_value))
    )
  ) +
    geom_point(aes(size = intersection_size, color = recall), alpha = 0.7) +
    scale_size_continuous(range = c(3, 10)) +
    scale_color_gradient(low = "blue", high = "red") +
    labs(
      x = "-log10(p-value)",
      y = "KEGG Term",
      size = "intersection_size",
      color = "Recall",
      title = paste0("Top 25 KEGG Enrichment Pathways, ", name)
    ) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10))

  # Save the bubble plots to the output directory
  ggsave(paste0("output_images/kegg_bubble_recall_", name, ".pdf"),
    plot = bubble_plot, width = 10, height = 8, bg = "white"
  )
  ggsave(paste0("output_images/kegg_bubble_pvalue_", name, ".pdf"),
    plot = bubble_plot_p_value, width = 10, height = 8, bg = "white"
  )
}

#' Filter gene scores and run enrichment analysis across up/down/all gene sets.
#'
#' This helper function reads a gene score CSV file, filters genes based on
#' correlation magnitude and a specified q-value threshold, and then runs
#' KEGG enrichment analysis for three gene subsets: all significant genes,
#' positively correlated genes, and negatively correlated genes.
#'
#' @param file Character string specifying the path to the input CSV file.
#' @param qval_col Character string giving the column name for the q-value
#'   (adjusted p-value) used for filtering significance.
#' @param name_prefix Character string used as a prefix for output file names
#'   in enrichment plots.
#' @param threshold Numeric value defining the absolute correlation cutoff
#'   for selecting genes (default = 0.7).
#' @param qval_cutoff Numeric value defining the maximum q-value threshold
#'   for significance filtering (default = 0.05).
#'
#' @return NULL
#'   This function does not return a value. It triggers enrichment analysis
#'   and saves output plots to disk via `run_enrichment()`.
#'
#' @examples
#' process_and_run(
#'   file = alpha_driver_genes_human.csv,
#'   qval_col = "Alpha_qval",
#'   name_prefix = "human_alpha"
#' )
process_and_run <- function(file,
                            qval_col,
                            name_prefix,
                            threshold = 0.7,
                            qval_cutoff = 0.05) {
  # Read in the investigation specific dataframe
  df <- read.csv(here(file)) |> na.omit()

  # Calculate all significant genes, as well as up and downregulated genes only
  filtered_all <- df[abs(df$scaled_correlation) > threshold &
    df[[qval_col]] < qval_cutoff, ]$X
  filtered_up <- df[df$scaled_correlation > threshold &
    df[[qval_col]] < qval_cutoff, ]$X
  filtered_down <- df[df$scaled_correlation < -threshold &
    df[[qval_col]] < qval_cutoff, ]$X

  # Run enrichment for each of the three gene sets of interest
  run_enrichment(filtered_all,
    background_genes,
    name = paste0(name_prefix, "_all")
  )
  run_enrichment(filtered_up,
    background_genes,
    name = paste0(name_prefix, "_up")
  )
  run_enrichment(filtered_down,
    background_genes,
    name = paste0(name_prefix, "_down")
  )
}

# Load in the background genes common in both objects
background_genes <- read.csv(here("csv_files/background_genes.csv"))$Genes

# Run the processing for human acinar driving genes
process_and_run(
  "go_files/acinar_maturing_driver_genes_human_scaled.csv",
  "Acinar.Maturing_qval",
  "human_acinar"
)

# Run the processing for human ductal driving genes
process_and_run(
  "go_files/ductal_immature_1_driver_genes_human_scaled.csv",
  "Ductal.Immature_1_qval",
  "human_ductal"
)

# Run the processing for ferret acinar driving genes
process_and_run(
  "go_files/acinar_maturing_driver_genes_ferret_scaled.csv",
  "Acinar.Maturing_qval",
  "ferret_acinar"
)

# Run the processing for ferret ductal driving genes
process_and_run(
  "go_files/ductal_proliferating_driver_genes_ferret_scaled.csv",
  "Ductal.Proliferating_qval",
  "ferret_ductal"
)
