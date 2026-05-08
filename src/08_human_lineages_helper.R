#' Get Query Gene Names for a Specified Cluster
#'
#' Retrieves the significant differentially expressed genes for a given cluster
#' from a Seurat object, mapping them to their human orthologs using a species
#' mapping file.
#'
#' @param seurat_object A Seurat object containing single-cell RNA-seq data.
#' @param cluster_number An integer specifying the cluster number for which to
#'   retrieve differentially expressed genes.
#'
#' @return A character vector of gene names associated with the specified
#'   cluster, including both positively and negatively expressed genes.
#'
#' @examples
#' query_genes <- get_query_names(
#'   seurat_object = my_seurat,
#'   cluster_number = 3
#' )
get_query_names <- function(seurat_object,
                            cluster_number) {
  # Set seed for reproducibility
  set.seed(0)

  # Set parameters for helper functions
  find_conversion <- FALSE
  all_seurat_genes <- FALSE
  only_positive <- FALSE

  # Load Seurat Object
  mapping_file <- read.csv(here("species_mapping_file.csv"))

  # Get all genes mapped to their ortholog
  all_genes <- get_all_genes(
    seurat_object = seurat_object,
    mapping_file = mapping_file,
    mapping_ortholog_col = "Human.gene.stable.ID",
    find_conversions = find_conversion,
    all_seurat_genes = all_seurat_genes
  )

  # Get the query genes from the significant cluster expression
  query_names <- get_cluster_genes(
    all_genes = all_genes,
    cluster_number = cluster_number,
    only_positive = only_positive,
    conversions = (find_conversion || all_seurat_genes),
    seurat_object = seurat_object
  )

  return(query_names)
}

#' Fit and Embed Slingshot Trajectory Curves
#'
#' Runs Slingshot trajectory inference on a SingleCellExperiment object and
#' projects the resulting principal curves into a specified dimensionality
#' reduction embedding (e.g. UMAP) from a Seurat object.
#'
#' @param seurat_object A Seurat object containing the target dimensionality
#'   reduction embedding.
#' @param sce A SingleCellExperiment object with cluster labels stored
#'   in sce$seurat_clusters and a reduced dimension matrix accessible
#'   via reducedDim().
#' @param reduced_dim A string specifying the name of the reduced dimension in
#'   sce to use for trajectory inference (e.g. "PCA").
#' @param reduction A string specifying the name of the dimensionality
#'   reduction in seurat_object to embed the curves into (e.g. "umap").
#' @param start_cluster An optional integer or string specifying the cluster
#'   to use as the root of the trajectory. Default is NULL (inferred
#'   automatically by Slingshot).
#' @param end_cluster An optional integer or string specifying the cluster(s)
#'   to use as terminal states of the trajectory. Default is NULL
#'   (inferred automatically by Slingshot).
#'
#' @return A named list with two elements:
#'   - all_curves: A list of principal curve objects embedded in the
#'     target reduction space, as returned by slingshot slingCurves.
#'   - slingshot_sce: The SingleCellExperiment object returned by slingshot,
#'     containing pseudotime nd lineage information.
#' }
#'
#' @examples
#' curves <- get_curves(
#'   seurat_object = my_seurat,
#'   sce           = my_sce,
#'   reduced_dim   = "PCA",
#'   reduction     = "umap",
#'   start_cluster = "0",
#'   end_cluster   = "5"
#' )
get_curves <- function(seurat_object,
                       sce,
                       reduced_dim,
                       reduction,
                       start_cluster = NULL,
                       end_cluster = NULL) {
  slingshot_sce <- slingshot(
    sce,
    reducedDim = reduced_dim,
    clusterLabels = sce$Leiden_cluster,
    start.clus = start_cluster,
    end.clus = end_cluster
  )

  pseudotime <- slingPseudotime(slingshot_sce, na = FALSE)
  all_umap <- Embeddings(seurat_object, reduction = reduction)
  sce_umap <- all_umap[rownames(all_umap) %in% rownames(pseudotime), ]
  sce_umap_slingshot <- embedCurves(slingshot_sce, newDimRed = sce_umap)
  all_curves <- slingCurves(sce_umap_slingshot)

  return(list(all_curves = all_curves, slingshot_sce = slingshot_sce))
}

#' Generate and display a set of UMAP and violin plots for a specific lineage
#'
#' @param sce SingleCellExperiment object for the dataset
#' @param slingshot_sce Slingshot object for the dataset
#' @param all_curves List ofcurves corresponding to different lineages
#' @param seurat_object Seurat object containing UMAP embeddings
#' @param lineage The specific lineage to focus on
#' @param celltype If TRUE, only display plots related to celltype
getplots <- function(sce,
                     slingshot_sce,
                     all_curves,
                     seurat_object,
                     slingshot_seurat,
                     lineage,
                     dimension = "rna.umap",
                     celltype = FALSE,
                     pseudotime_limit = 150,
                     trim_index = NULL) {
  # Specify the seurat object columns to be plotted
  output_plot_columns <- c(
    "Leiden_cluster",
    "PCW",
    "Celltype"
  )

  # Extract the current curve for the specified lineage
  current_curve <- all_curves[[lineage]]
  options(warn = -1)

  # Generate UMAP plots for the given lineage
  lineage_plots <- umap_all_plots(
    curve = current_curve,
    seurat_object = seurat_object,
    palette = custom_palette,
    column_list = output_plot_columns,
    dimension = dimension,
    lineage = lineage
  )

  lineage_pseudotime <- umap_pseudotime(
    curve = current_curve,
    seurat_object = slingshot_seurat,
    original_object = seurat_object,
    slingshot_sce = slingshot_sce,
    dimension = dimension,
    lineage = lineage,
    pseudotime_limit = pseudotime_limit,
    trim_index = trim_index
  )

  # Generate violin plots for the given lineage
  violin_plots <- violin_all_plots(
    sce = sce,
    slingshot_sce = slingshot_sce,
    palette = custom_palette,
    column_list = output_plot_columns,
    lineage = lineage
  )

  # Extract the four pairs of UMAP and violin plots
  p1 <- lineage_plots[[1]]
  p2 <- violin_plots[[1]]
  p3 <- lineage_plots[[2]]
  p4 <- violin_plots[[2]]
  p5 <- lineage_plots[[3]]
  p6 <- violin_plots[[3]]

  # Create grids of UMAP and violin plots side by side
  grid_1 <- plot_grid(p1, p2,
    labels = c("A", "B"),
    ncol = 2, rel_widths = c(0.9, 1)
  )
  grid_2 <- plot_grid(p3, p4,
    labels = c("C", "D"),
    ncol = 2, rel_widths = c(0.85, 1)
  )
  grid_3 <- plot_grid(p5, p6,
    labels = c("E", "F"),
    ncol = 2, rel_widths = c(1.2, 1)
  )

  print(grid_1)
  print(grid_2)
  print(grid_3)

  options(warn = 0)
  return(lineage_pseudotime)
}

#' Generate pseudotime-ordered heatmaps for lineage-associated genes
#'
#' This function filters lineage-associated differential expression results,
#' identifies significant genes based on p-value and Wald statistic thresholds,
#' and generates a heatmap of gene expression ordered along pseudotime for a
#' specified lineage.
#'
#' @param slingshot_sce A SingleCellExperiment object containing Slingshot
#'        results.
#' @param object A Seurat object containing expression data.
#' @param gene_list A data frame or list containing genes of interest
#'        (expects a `genes` column).
#' @param lineage Integer or character specifying the lineage to analyze.
#' @param lineage_asso_res Data frame of lineage association results
#'        (e.g., from tradeSeq).
#' @param waldcutoff Numeric threshold for filtering genes by Wald statistic
#'        (default = 200).
#' @param pval_cutoff Numeric threshold for filtering genes by p-value
#'        (default = 0.01).
#' @param save_dir Directory path where output CSV files will be saved.
#' @param fraction Numeric (0–1) indicating the fraction of pseudotime to retain
#'        (default = 1.0).
#' @param day Logical; if TRUE, include "day" metadata in the heatmap, otherwise
#'        use cluster labels.
#'
#' @return sig_genes_plot: Heatmap data or object for significant genes ordered
#'         by pseudotime
#'
#' @examples
#' results <- make_lineage_heatmaps(
#'   slingshot_sce = sce_obj,
#'   object = seurat_obj,
#'   gene_list = marker_df,
#'   lineage = 1,
#'   lineage_asso_res = assoc_res,
#'   waldcutoff = 200,
#'   pval_cutoff = 0.01,
#'   save_dir = "results/heatmaps",
#'   fraction = 0.8,
#'   day = TRUE
#' )
make_lineage_heatmaps <- function(slingshot_sce,
                                  object,
                                  gene_list,
                                  lineage,
                                  lineage_asso_res,
                                  waldcutoff = 200,
                                  pval_cutoff = 0.01,
                                  save_dir = save_dir,
                                  fraction = 1.0,
                                  day = TRUE) {
  # Which columns to extract from the dataframe
  keep_cols <- c("waldStat", "pvalue", "df")
  keep_cols <- paste0(keep_cols, "_", lineage)

  # Filter by pvalue
  filtered_res <- lineage_asso_res %>%
    dplyr::select(dplyr::all_of(keep_cols)) %>%
    setNames(gsub("_[0-9]*", "", names(.))) %>%
    dplyr::filter(pvalue < pval_cutoff) %>%
    dplyr::arrange(-waldStat)

  write.csv(filtered_res, file.path(save_dir, "fitgam_sig_genes.csv"))

  rownames(filtered_res)

  opts <- list(species = 9606)
  tfs <- getMatrixSet(JASPAR2020, opts)
  tf_gene_symbols <- sapply(tfs, function(tf) tf@name)
  filtered_res$tf_mapped <- rownames(filtered_res) %in% tf_gene_symbols
  deg_tfs <- filtered_res[filtered_res$tf_mapped, ]

  # Filter by gene list
  filtered_genes_res <-
    filtered_res[rownames(filtered_res) %in% gene_list$genes, ]
  write.csv(filtered_genes_res, file.path(
    save_dir,
    "fitgam_findmarkers_sig_genes.csv"
  ))

  # Filter by wald cutoff (fewer genes)
  filtered_res_high <- filtered_res %>%
    dplyr::filter(waldStat > waldcutoff)

  write.csv(filtered_res_high, file.path(
    save_dir,
    "fitgam_waldcutoff_genes.csv"
  ))

  pseudotime <- slingPseudotime(slingshot_sce, na = TRUE)[, lineage]

  # Ensure pseudotime is a named vector
  pseudotime <- as.vector(pseudotime)
  names(pseudotime) <- colnames(slingshot_sce)

  # Remove NA values
  valid_indices <- !is.na(pseudotime)
  pseudotime <- pseudotime[valid_indices]

  # Get the cell names for the valid pseudotime values
  cells <- names(pseudotime)

  # Keep only cells in the lineage
  seurat_sub <- subset(object, cells = cells)

  # Add pseudotime to the Seurat metadata
  seurat_sub <- AddMetaData(seurat_sub,
    metadata = pseudotime,
    col.name = "pseudotime"
  )

  # Order cells by pseudotime
  ordered_cells <- seurat_sub@meta.data %>%
    dplyr::arrange(pseudotime) %>%
    rownames()

  pseudotime_cutoff <- fraction

  ordered_cells_subset <-
    ordered_cells[pseudotime[ordered_cells] <= pseudotime_cutoff]

  # Subset the Seurat object based on the selected cells
  seurat_sub <- seurat_sub[, ordered_cells_subset]


  top_genes <- rownames(head(filtered_res, 200))
  expr_mat <- GetAssayData(seurat_sub, layer = "data", assay = "RNA")
  expr_mat <- expr_mat[top_genes, ordered_cells_subset, drop = FALSE]

  pseudotime_vals <- pseudotime[ordered_cells_subset]

  window_size <- 10 # adjust based on dataset size
  n <- length(pseudotime_vals)

  gene_peak <- apply(expr_mat, 1, function(x) {
    if (sum(x) == 0) {
      return(Inf)
    } # skip empty genes

    # compute rolling sums over pseudotime
    window_sums <- sapply(seq_len(n - window_size + 1), function(i) {
      sum(x[i:(i + window_size - 1)])
    })

    # pick the pseudotime midpoint of the max window
    best_idx <- which.max(window_sums)
    mean(pseudotime_vals[best_idx:(best_idx + window_size - 1)])
  })

  # order genes by windowed peak pseudotime
  gene_order <- names(sort(gene_peak))

  seurat_sub$pseudotime <- pseudotime_vals

  if (day) {
    meta_col <- c("day", "seurat_clusters")
  } else {
    meta_col <- c("pseudotime")
  }

  # all sig genes plot
  sig_genes_plot <- plot_heatmap(seurat_sub,
    gene_list = gene_order,
    meta_col = meta_col,
    cluster_rows = FALSE,
    plot_rownames = FALSE,
    cell_order = ordered_cells_subset,
    return_data = TRUE
  )

  png(file.path(save_dir, "fitgam_sig_genes.png"),
    width = 1600, height = 1000, res = 125
  )
  print(sig_genes_plot$heatmap)
  dev.off()

  return(list(
    sig_genes_plot = sig_genes_plot,
    cluster_sig_plot = NULL
  ))
}

#' Get all genes with orthologs
#'
#' Extract all genes from a Seurat object and attempt to map them to orthologs
#' using the inputted mapping file. Optionally, it can use the gconvert function
#' to find orthologs for genes that are not present in the mapping file.
#'
#' @param seurat_object A Seurat object containing the gene expression data.
#' @param mapping_file A data frame containing the gene mapping information. It
#' should have columns for gene IDs and their corresponding orthologs.
#' @param mapping_ortholog_col A string specifying the column name in the
#' mapping file that contains the orthologs.
#' @param find_conversions A boolean indicating whether to use gconvert to find
#' orthologs for genes not present in the mapping file. Defaults to FALSE.
#' @param all_seurat_genes A boolean indicating whether to attempt to convert
#' all genes from the Seurat object using gconvert, regardless of their presence
#' in the mapping file. This parameter is only considered if find_conversions
#' is set to FALSE. Defaults to FALSE.
#'
#' @return A data frame with two columns: "genes" and "orthologs".
#' Genes without orthologsare excluded from the result.
#'
#' @examples
#' mapping_df <- data.frame(
#'   gene_id = c("GeneA", "GeneB"),
#'   human_ortholog = c("HUMAN1", "HUMAN2")
#' )
get_all_genes <- function(seurat_object,
                          mapping_file,
                          mapping_ortholog_col,
                          find_conversions = FALSE,
                          all_seurat_genes = FALSE) {
  # Extract the gene_id column from the mapping file
  mapping_gene_col <- "gene_id"

  # Get orthologs for all genes first
  all_genes <- data.frame("genes" = rownames(seurat_object)) %>%
    dplyr::filter(genes %in% mapping_file[[mapping_gene_col]])

  mapping_list <- mapping_file[[mapping_ortholog_col]]
  names(mapping_list) <- mapping_file[[mapping_gene_col]]

  # Create column of orthologs next to gene names
  all_genes$orthologs <- mapping_list[all_genes$genes]

  if (find_conversions || all_seurat_genes) {
    # Extract genes that didn't have a match
    no_orthologs <- all_genes %>%
      filter(is.na(orthologs) | orthologs == "")

    genes_to_convert <- if (find_conversions) no_orthologs$genes else rownames(seurat_object)
    gconvert_result <- gconvert(query = genes_to_convert)

    # Reformat for merging with original all_genes object
    gconvert_result <- gconvert_result[, c("input", "target")]
    colnames(gconvert_result) <- c("genes", "new_orthologs")

    if (find_conversions) {
      all_genes <- merge(all_genes, gconvert_result, by = "genes", all.x = TRUE)

      all_genes$orthologs <-
        ifelse(is.na(all_genes$orthologs) | all_genes$orthologs == "",
          all_genes$new_orthologs, all_genes$orthologs
        )
    } else {
      all_genes <- merge(gconvert_result, all_genes, by = "genes", all.x = TRUE)

      all_genes$orthologs <-
        ifelse(is.na(all_genes$new_orthologs) | all_genes$new_orthologs == "",
          all_genes$orthologs, all_genes$new_orthologs
        )
    }

    # Remove temporary column
    all_genes$new_orthologs <- NULL
  }

  # Remove all rows that are NA or blank
  filtered_all_genes <- all_genes %>%
    filter(!is.na(orthologs) & orthologs != "")

  return(filtered_all_genes)
}

#' Get significant genes with orthologs for a given cluster
#'
#' This function identifies significant genes for a specified cluster in a
#' Seurat object and attempts to map them to orthologs using the provided
#' mapping file. Optionally, it can use the gconvert function to find orthologs
#' for genes that are not mapped.
#'
#' @param all_genes A data frame containing all genes and their orthologs. It
#' must contain columns for gene names and their corresponding orthologs.
#' @param cluster_number An integer specifying the cluster number for which to
#' find significant genes.
#' @param only_positive A boolean indicating whether to include only positive
#' markers. Defaults to FALSE.
#' @param conversions A boolean indicating whether to use gconvert to find
#' orthologs for genes that are not present in the mapping file. Defaults to
#' FALSE.
#'
#' @return A vector of orthologs for the significant genes in the input cluster.
#'
#' @examples
#' gene_map <- data.frame(
#'   genes = c("GeneA", "GeneB", "GeneC"),
#'   orthologs = c("HUMAN1", "HUMAN2", NA)
#' )
get_cluster_genes <- function(all_genes,
                              cluster_number,
                              only_positive = FALSE,
                              conversions = FALSE,
                              seurat_object) {
  # Get gene names for cluster
  cluster_markers <- FindMarkers(seurat_object,
    ident.1 = cluster_number,
    only.pos = only_positive
  )

  # Extract the markers that are significant only
  significant_markers <- subset(cluster_markers, p_val_adj < 0.05)

  gene_names <- rownames(significant_markers)
  ordered_genes <- data.frame(genes = gene_names)

  # If we want to run gconvert to get unmapped genes, do that here
  if (conversions) {
    significant_genes <- ordered_genes %>%
      left_join(all_genes, by = "genes")

    # Gene the unmapped genes for gconvert
    no_orthologs <- significant_genes %>%
      filter(is.na(orthologs) | orthologs == "")

    gconvert_result <- gconvert(query = no_orthologs$genes)

    # Check to see if we have success with any of the mappings
    if (!is.null(gconvert_result)) {
      # Reformat for merging with original all_genes object
      gconvert_result <- gconvert_result[, c("input", "target")]
      colnames(gconvert_result) <- c("genes", "new_orthologs")

      # Join the columns correctly
      significant_genes <- significant_genes %>%
        left_join(gconvert_result, by = "genes")

      significant_genes$orthologs <-
        ifelse(is.na(significant_genes$orthologs) |
          significant_genes$orthologs == "",
        significant_genes$new_orthologs, significant_genes$orthologs
        )

      # Remove temporary column
      significant_genes$new_orthologs <- NULL

      # Filter out any NA's
      significant_genes <- significant_genes %>%
        filter(!is.na(orthologs))
    }
  } else {
    # Filter orthologs for significant cluster
    significant_genes <- ordered_genes %>%
      left_join(all_genes, by = "genes") %>%
      filter(!is.na(orthologs))
  }

  return(significant_genes)
}


#' Generate multiple UMAP plots with lineage curve overlay
#'
#' This function creates a list of UMAP plots from a Seurat object, coloring
#' cells by specified metadata columns and overlaying a lineage curve (e.g.,
#' from Slingshot). Each plot corresponds to a different metadata column.
#'
#' @param curve A Slingshot curve object containing lineage coordinates.
#' @param seurat_object Seurat object containing dimensional reduction results.
#' @param palette A vector of colors used for non-"day" metadata columns.
#' @param column_list A character vector of metadata columns to visualize.
#' @param dimension Character specifying which dimensional reduction to use
#'   (default = "rna.umap").
#' @param lineage Integer or identifier for the lineage (currently unused but
#'   included for consistency with upstream workflows).
#'
#' @return A list of ggplot objects, for each metadata column in `column_list`.
#'
#' @examples
#' plots <- umap_all_plots(
#'   curve = sling_curves[[1]],
#'   seurat_object = seurat_obj,
#'   palette = c("red", "blue", "green", "purple"),
#'   column_list = c("day", "Leiden_cluster"),
#'   dimension = "rna.umap",
#'   lineage = 1
#' )
umap_all_plots <- function(curve,
                           seurat_object,
                           palette,
                           column_list,
                           dimension = "rna.umap",
                           lineage) {
  # Create the umaps using the specified columns in the parameter list
  all_plots <- lapply(1:length(column_list), function(x) {
    if (column_list[x] == "day") {
      cur_palette <- c("#4285F4", "#34A853", "#FBBC05", "#EA4335")
    } else {
      cur_palette <- palette
    }

    current_column <- column_list[[x]]

    curve1_coord <- data.frame(curve$s[curve$ord, c(1, 2)])
    curve1_coord$stage <- "line"

    base_plot <- plotDimRed(seurat_object,
      col_by = current_column, color = cur_palette,
      plot_type = dimension,
      ggrastr = TRUE
    )[[1]]

    if (dimension == "rna.umap") {
      base_plot <- base_plot + ggplot2::geom_path(
        data = curve1_coord,
        ggplot2::aes(rnaUMAP_1, rnaUMAP_2),
        color = "black", size = 1
      )
    } else {
      base_plot <- base_plot + ggplot2::geom_path(
        data = curve1_coord,
        ggplot2::aes(UMAP_1, UMAP_2),
        color = "black", size = 1
      )
    }

    base_plot <- base_plot + theme(
      legend.text = element_text(size = 6), # Legend text size
      legend.key.size = unit(0.5, "cm") # Legend key size
    )

    return(base_plot)
  })

  return(all_plots)
}

#' Plot UMAP colored by pseudotime with lineage overlay
#'
#' This function generates a UMAP plot where cells are colored by pseudotime
#' values derived from Slingshot. It maps pseudotime from a
#' SingleCellExperiment object onto a Seurat object, applies an optional
#' pseudotime cutoff, and overlays the lineage trajectory curve.
#'
#' @param curve A Slingshot curve object containing lineage coordinates.
#' @param slingshot_sce A SingleCellExperiment object containing Slingshot
#'        results.
#' @param seurat_object A Seurat object used to align pseudotime values.
#' @param original_object A Seurat object to visualize (receives pseudotime
#'        metadata).
#' @param dimension Character specifying which dimensional reduction to use
#'   (default = "rna.umap").
#' @param lineage Integer or character specifying the lineage to analyze.
#' @param pseudotime_limit Numeric value to cap pseudotime; values above this
#'   threshold are set to NA (default = 150).
#'
#' @return A ggplot (cowplot grid) object showing UMAP colored by pseudotime
#'   with the lineage trajectory overlaid.
#'
#' @examples
#' pt_plot <- umap_pseudotime(
#'   curve = sling_curves[[1]],
#'   slingshot_sce = sce_obj,
#'   seurat_object = seurat_obj,
#'   original_object = full_seurat_obj,
#'   dimension = "rna.umap",
#'   lineage = 1,
#'   pseudotime_limit = 120
#' )
umap_pseudotime <- function(curve,
                            slingshot_sce,
                            seurat_object,
                            original_object,
                            dimension = "rna.umap",
                            lineage,
                            pseudotime_limit = 150,
                            trim_index = NULL) {
  pseudotime_column <- paste0("slingPseudotime_", lineage)
  pseudotime_values <- colData(slingshot_sce)[[pseudotime_column]]
  names(pseudotime_values) <- colnames(seurat_object)

  original_object$pseudotime <- rep(NA, ncol(original_object))
  matching_cells <- intersect(
    colnames(seurat_object),
    colnames(original_object)
  )
  original_object$pseudotime[matching_cells] <-
    pseudotime_values[matching_cells]

  original_object$pseudotime[original_object$pseudotime > pseudotime_limit] <- NA

  # Create UMAP plot colored by pseudotime
  curve1_coord <- data.frame(curve$s[curve$ord, c(1, 2)])
  curve1_coord$stage <- "line"

  if (!is.null(trim_index)) {
    curve1_coord <- curve1_coord[1:trim_index, ]
  }

  umap_coords <- Embeddings(original_object, reduction = dimension)
  umap_df <- data.frame(
    UMAP_1 = umap_coords[, 1],
    UMAP_2 = umap_coords[, 2],
    pseudotime = original_object$pseudotime
  )

  # Create a base plot using ggplot2 and scale_color_gradient
  base_plot <- ggplot(umap_df, aes(
    x = UMAP_1,
    y = UMAP_2,
    color = pseudotime
  )) +
    geom_point_rast(size = 1) + # Plot the points with the pseudotime coloring
    scale_color_gradientn(
      colors = c("blue", "green", "yellow", "red"),
      na.value = "#D3D3D3", # Handle missing values
      name = "Pseudotime"
    ) + # Legend name
    labs(
      title = "UMAP colored by Pseudotime",
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    theme(legend.position = "right")

  # Add the trajectory line on top
  if (dimension == "rna.umap") {
    base_plot <- base_plot + ggplot2::geom_path(
      data = curve1_coord,
      ggplot2::aes(rnaUMAP_1, rnaUMAP_2),
      color = "black", size = 1
    )
  } else {
    base_plot <- base_plot + ggplot2::geom_path(
      data = curve1_coord,
      ggplot2::aes(UMAP_1, UMAP_2),
      color = "black", size = 1
    )
  }
  grid_plot <- plot_grid(base_plot, NULL,
    labels = c("I", ""),
    ncol = 2, rel_widths = c(2, 1),
    rel_heights = c(3, 1)
  )

  return(grid_plot)
}


#' Generate violin-style pseudotime distribution plots for multiple metadata
#' columns
#'
#' This function creates a series of quasi-random (violin-like) plots showing
#' the distribution of cells along Slingshot pseudotime, grouped and colored by
#' specified metadata columns. Each plot corresponds to a different column.
#'
#' @param sce A SingleCellExperiment object containing cell metadata.
#' @param slingshot_sce A SingleCellExperiment object containing Slingshot
#'        results.
#' @param palette A vector of colors used for non-"day" metadata columns.
#' @param column_list A character vector of metadata columns to visualize.
#' @param lineage Integer or character specifying the lineage to analyze.
#'
#' @return A list of ggplot objects, each representing pseudotime distributions
#'   for a metadata column.
#'
#' @examples
#' violin_plots <- violin_all_plots(
#'   sce = sce_obj,
#'   slingshot_sce = sce_obj,
#'   palette = c("red", "blue", "green"),
#'   column_list = c("day", "Leiden_cluster"),
#'   lineage = 1
#' )
violin_all_plots <- function(sce,
                             slingshot_sce,
                             palette,
                             column_list,
                             lineage) {
  # Specify the pseudotime column to utilize
  pseudotime_col <- paste0("slingPseudotime_", lineage)
  slingshot_sce$current_sling <- slingshot_sce[[pseudotime_col]]

  # Iterate over the pseudotime columns
  all_plots <- lapply(1:length(column_list), function(i) {
    if (column_list[i] == "day") {
      cur_palette <- c("#4285F4", "#34A853", "#FBBC05", "#EA4335")
    } else {
      cur_palette <- palette
    }

    column <- column_list[[i]]

    unique_samples <- sort(unique(sce[[column]]))
    unique_samples_rev <- rev(unique_samples)

    sce[[column]] <- factor(sce[[column]], levels = unique_samples)
    sce[["reverse"]] <- factor(sce[[column]], levels = unique_samples_rev)

    # Create the plot for the current pseudotime
    plot <- ggplot(as.data.frame(colData(sce)), aes(
      x = slingshot_sce$current_sling,
      y = sce[["reverse"]],
      colour = sce[[column]]
    )) +
      geom_quasirandom(groupOnX = FALSE, orientation = "y") +
      scale_color_manual(values = cur_palette, name = column, drop = FALSE) +
      theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none"
      ) +
      xlab("Slingshot pseudotime")

    return(plot)
  })

  return(all_plots)
}

#' NOTE: This is the same function as was used in scAnalysis, but tweaked so the
#' labels/colors could be more customizable when used with make_lineage_heatmaps
#'
#' @param seurat_object A Seurat object containing the single-cell RNA-seq data.
#' @param gene_list A character vector of gene names to be included in the
#'        heatmap.
#' @param meta_col A character string specifying the metadata column to use for
#'        grouping cells.
#' @param colors Optional. A named vector of colors to use for the metadata
#'        groups.
#' @param meta_df Optional. A data frame containing metadata information.
#' @param color_list Optional. A list of named color vectors for metadata
#'        columns.
#' @param max_val Numeric. The maximum value for scaling the heatmap.
#'        Default is 2.5.
#' @param min_val Numeric. The minimum value for scaling the heatmap.
#'        Default is -2.5.
#' @param cluster_rows Logical. Whether to cluster rows. Default is FALSE.
#' @param cluster_cols Logical. Whether to cluster columns. Default is FALSE.
#' @param average_expression Logical. Whether to average expression by the
#'        metadata groups. Default is FALSE.
#' @param plot_meta_col Logical. Whether to plot the metadata column as an
#'        annotation. Default is TRUE.
#' @param plot_rownames Logical. Whether to display row names (gene names) in
#'        the heatmap. Default is TRUE.
#' @param cell_order Optional. A character vector specifying the order of cells.
#' @param return_data Logical. Whether to return the heatmap object and the
#'        underlying data. Default is FALSE.
#' @param assay Character. The assay to use from the Seurat object.
#'        Default is "RNA".
#' @param ... Additional parameters to pass to the pheatmap function.
plot_heatmap <- function(seurat_object,
                         gene_list,
                         meta_col,
                         colors = NULL,
                         fontsize_row = 5,
                         meta_df = NULL,
                         color_list = NULL,
                         max_val = 2.5,
                         min_val = -2.5,
                         cluster_rows = FALSE,
                         cluster_cols = FALSE,
                         average_expression = FALSE,
                         plot_meta_col = TRUE,
                         plot_rownames = TRUE,
                         cell_order = NULL,
                         return_data = FALSE,
                         assay = "RNA", ...) {
  if (average_expression) {
    Idents(seurat_object) <- meta_col
    heatmap_df <- AverageExpression(seurat_object,
      seurat = FALSE,
      group.by = meta_col, assays = assay
    )
    heatmap_df <- heatmap_df[[assay]]
    character_vals <- suppressWarnings(all(!is.na(as.numeric(
      as.character(colnames(heatmap_df))
    ))))
    if (is.null(meta_df)) {
      sample_info <- seurat_object[[meta_col]]
      if (is.null(levels(sample_info[[meta_col]]))) {
        sample_info[[meta_col]] <- factor(sample_info[[meta_col]])
      }
      meta_df <- data.frame(levels(sample_info[[meta_col]]))
      colnames(meta_df) <- meta_col
      if (character_vals) {
        rownames(meta_df) <- paste0("X", meta_df[[meta_col]])
      } else {
        rownames(meta_df) <- meta_df[[meta_col]]
      }
      if (is.null(colors)) {
        colorcount <- length(levels(sample_info[[meta_col]]))
        colors <- (grDevices::colorRampPalette(RColorBrewer::brewer.pal(
          9,
          "Set1"
        )))(colorcount)
        names(colors) <- levels(sample_info[[meta_col]])
      }
      color_list <- list(colors)
      names(color_list) <- meta_col
    }
  } else {
    heatmap_df <- GetAssayData(seurat_object,
      layer = "data",
      assay = assay
    )
  }
  heatmap_df <- heatmap_df[rownames(heatmap_df) %in% gene_list, ]
  heatmap_df <- data.frame(heatmap_df)
  heatmap_df <- heatmap_df[order(match(
    rownames(heatmap_df),
    gene_list
  )), ]
  heatmap_df <- heatmap_df[rowSums(heatmap_df) > 0, ]
  if (is.null(meta_df)) {
    sample_info <- seurat_object[[meta_col]]
  } else {
    sample_info <- meta_df
    coloring <- color_list
  }
  heatmap_scale <- t(scale(t(heatmap_df), scale = TRUE))
  blueYellow <- c(
    "#352A86", "#343DAE", "#0262E0", "#1389D2",
    "#2DB7A3", "#A5BE6A", "#F8BA43", "#F6DA23", "#F8FA0D"
  )
  if (!is.null(cell_order)) {
    sample_info <- sample_info[order(match(
      rownames(sample_info),
      cell_order
    )), , drop = FALSE]
    rownames(sample_info) <- make.names(rownames(sample_info))
    if (!identical(colnames(heatmap_scale), rownames(sample_info))) {
      heatmap_scale <- heatmap_scale[, rownames(sample_info)]
    }
  } else if (!cluster_cols) {
    cluster_order <- levels(sample_info[[meta_col]])
    sample_info <- sample_info[order(match(
      sample_info[[meta_col]],
      cluster_order
    )), , drop = FALSE]
    rownames(sample_info) <- make.names(rownames(sample_info))
    if (!identical(colnames(heatmap_scale), rownames(sample_info))) {
      heatmap_scale <- heatmap_scale[, rownames(sample_info)]
    }
  } else {
    rownames(sample_info) <- make.names(rownames(sample_info))
  }
  if (!plot_meta_col) {
    sample_info[[meta_col]] <- NULL
  }
  heatmap_scale <- ifelse(heatmap_scale > max_val, max_val,
    heatmap_scale
  )
  heatmap_scale <- ifelse(heatmap_scale < min_val, min_val,
    heatmap_scale
  )

  annotation_colors <- list(
    pseudotime = colorRampPalette(c("blue", "green", "yellow", "red"))(100)
  )

  ensg_ids <- rownames(heatmap_scale)
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

  gene_mapping <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name"),
    filters = "ensembl_gene_id",
    values = ensg_ids,
    mart = mart
  )

  ens_to_gene <- setNames(
    gene_mapping$external_gene_name,
    gene_mapping$ensembl_gene_id
  )
  output_mapping <- ifelse(rownames(heatmap_scale) %in% names(ens_to_gene) &
    ens_to_gene[rownames(heatmap_scale)] != "",
  ens_to_gene[rownames(heatmap_scale)],
  rownames(heatmap_scale)
  )

  rownames(heatmap_scale) <- as.character(output_mapping)

  heatmap <- pheatmap::pheatmap(heatmap_scale,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols, show_rownames = plot_rownames,
    show_colnames = FALSE, annotation_col = sample_info,
    annotation_colors = annotation_colors,
    color = blueYellow, border_color = NA, legend = FALSE,
    annotation_legend = FALSE, fontsize_row = fontsize_row,
    clustering_method = "complete", silent = TRUE, ...
  )
  if (return_data) {
    return(list(
      heatmap = heatmap, z_score = heatmap_scale,
      counts = heatmap_df
    ))
  } else {
    return(heatmap)
  }
}
