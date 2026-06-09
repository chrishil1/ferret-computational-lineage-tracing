# -----------------------------------------------------------------------------
# Dataset Configuration
# -----------------------------------------------------------------------------

#' Get dataset-specific plot column names
#'
#' Returns the metadata columns used for UMAP and violin plots depending on
#' the dataset. Human uses seurat_clusters/day/original_labels/majority_celltype;
#' ferret uses Leiden_cluster/PCW/Celltype.
#'
#' @param dataset Character; either "human" or "ferret".
#' @return A character vector of metadata column names.
get_plot_columns <- function(dataset = "human") {
  if (dataset == "human") {
    c("Leiden_cluster", "PCW", "Celltype")
  } else {
    c("seurat_clusters", "day", "original_labels", "final_renamed_celltype")
  }
}

#' Get the day/time column name for a dataset
#'
#' @param dataset Character; either "human" or "ferret".
#' @return A single column name string.
get_day_column <- function(dataset = "human") {
  if (dataset == "human") "day" else "PCW"
}

#' Get the day palette for a dataset
#'
#' @return A character vector of hex color codes.
get_day_palette <- function() {
  c("#4285F4", "#34A853", "#FBBC05", "#EA4335")
}

# -----------------------------------------------------------------------------
# Palette Utilities
# -----------------------------------------------------------------------------

#' Build a named color palette for a metadata column
#'
#' Creates a named vector mapping each unique level in a Seurat object metadata
#' column to a color from the supplied palette. Levels are derived from the
#' full Seurat object so that color assignments remain consistent even when a
#' downstream subset (e.g. a lineage SCE) contains only a fraction of cell
#' types.
#'
#' @param seurat_object A Seurat object whose metadata defines all possible
#'   levels for column.
#' @param column Character; the metadata column name.
#' @param palette A character vector of colors. Must be at least as long as the
#'   number of unique levels.
#'
#' @return A named character vector of colors, one per unique level.
#'
#' @examples
#' pal <- make_named_palette(seurat_obj, "majority_celltype", custom_palette)
make_named_palette <- function(seurat_object, column, palette) {
  all_levels <- sort(unique(seurat_object[[column]][, 1]))
  setNames(palette[seq_along(all_levels)], all_levels)
}

#' Build named palettes for all plot columns
#'
#' Calls make_named_palette for every column returned by get_plot_columns, using
#' the day palette for the time column and custom_palette for all others.
#'
#' @param seurat_object A Seurat object containing all cells.
#' @param dataset Character; either "human" or "ferret".
#' @param custom_palette A character vector of colors for non-day columns.
#'
#' @return A named list of color vectors, one per plot column.
build_named_palettes <- function(seurat_object, dataset, custom_palette) {
  columns <- get_plot_columns(dataset)
  day_col <- get_day_column(dataset)
  
  palettes <- lapply(columns, function(col) {
    pal <- if (col == day_col) get_day_palette() else custom_palette
    make_named_palette(seurat_object, col, pal)
  })
  names(palettes) <- columns
  palettes
}

# -----------------------------------------------------------------------------
# Query Gene Helpers
# -----------------------------------------------------------------------------

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
#' query_genes <- get_query_names(seurat_object = my_seurat, cluster_number = 3)
get_query_names <- function(seurat_object, cluster_number) {
  set.seed(0)
  
  find_conversion  <- FALSE
  all_seurat_genes <- FALSE
  only_positive    <- FALSE
  
  mapping_file <- read.csv(file.path(csv_dir, "species_mapping_file.csv"))
  
  all_genes <- get_all_genes(
    seurat_object        = seurat_object,
    mapping_file         = mapping_file,
    mapping_ortholog_col = "Human.gene.stable.ID",
    find_conversions     = find_conversion,
    all_seurat_genes     = all_seurat_genes
  )
  
  query_names <- get_cluster_genes(
    all_genes      = all_genes,
    cluster_number = cluster_number,
    only_positive  = only_positive,
    conversions    = (find_conversion || all_seurat_genes),
    seurat_object  = seurat_object
  )
  
  return(query_names)
}

#' Get all genes with orthologs
#'
#' Extract all genes from a Seurat object and attempt to map them to orthologs
#' using the inputted mapping file. Optionally uses gconvert for unmapped genes.
#'
#' @param seurat_object A Seurat object containing the gene expression data.
#' @param mapping_file A data frame with gene mapping information containing
#'   columns for gene IDs and their corresponding orthologs.
#' @param mapping_ortholog_col Column name in mapping_file for orthologs.
#' @param find_conversions Logical; use gconvert for genes absent from the
#'   mapping file. Defaults to FALSE.
#' @param all_seurat_genes Logical; attempt to convert all Seurat genes via
#'   gconvert regardless of mapping file presence. Only considered when
#'   find_conversions is FALSE. Defaults to FALSE.
#'
#' @return A data frame with columns "genes" and "orthologs". Genes without
#'   orthologs are excluded.
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
  mapping_gene_col <- "gene_id"
  
  all_genes <- data.frame("genes" = rownames(seurat_object)) %>%
    dplyr::filter(genes %in% mapping_file[[mapping_gene_col]])
  
  mapping_list <- mapping_file[[mapping_ortholog_col]]
  names(mapping_list) <- mapping_file[[mapping_gene_col]]
  all_genes$orthologs <- mapping_list[all_genes$genes]
  
  if (find_conversions || all_seurat_genes) {
    no_orthologs <- all_genes %>%
      filter(is.na(orthologs) | orthologs == "")
    
    genes_to_convert <- if (find_conversions) no_orthologs$genes else rownames(seurat_object)
    gconvert_result  <- gconvert(query = genes_to_convert)
    
    gconvert_result <- gconvert_result[, c("input", "target")]
    colnames(gconvert_result) <- c("genes", "new_orthologs")
    
    if (find_conversions) {
      all_genes <- merge(all_genes, gconvert_result, by = "genes", all.x = TRUE)
      all_genes$orthologs <- ifelse(
        is.na(all_genes$orthologs) | all_genes$orthologs == "",
        all_genes$new_orthologs, all_genes$orthologs
      )
    } else {
      all_genes <- merge(gconvert_result, all_genes, by = "genes", all.x = TRUE)
      all_genes$orthologs <- ifelse(
        is.na(all_genes$new_orthologs) | all_genes$new_orthologs == "",
        all_genes$orthologs, all_genes$new_orthologs
      )
    }
    all_genes$new_orthologs <- NULL
  }
  
  all_genes %>% filter(!is.na(orthologs) & orthologs != "")
}

#' Get significant genes with orthologs for a given cluster
#'
#' Identifies significant markers for a cluster and maps them to orthologs.
#' Optionally uses gconvert for unmapped genes.
#'
#' @param all_genes A data frame with "genes" and "orthologs" columns.
#' @param cluster_number Integer cluster to query.
#' @param only_positive Logical; include only positive markers. Default FALSE.
#' @param conversions Logical; use gconvert for unmapped genes. Default FALSE.
#' @param seurat_object A Seurat object.
#'
#' @return A data frame of significant genes and their orthologs.
#'
#' @examples
#' gene_map <- data.frame(
#'   genes = c("GeneA", "GeneB", "GeneC"),
#'   orthologs = c("HUMAN1", "HUMAN2", NA)
#' )
get_cluster_genes <- function(all_genes,
                              cluster_number,
                              only_positive = FALSE,
                              conversions   = FALSE,
                              seurat_object) {
  cluster_markers <- FindMarkers(seurat_object,
                                 ident.1  = cluster_number,
                                 only.pos = only_positive
  )
  
  significant_markers <- subset(cluster_markers, p_val_adj < 0.05)
  ordered_genes <- data.frame(genes = rownames(significant_markers))
  
  if (conversions) {
    significant_genes <- ordered_genes %>%
      left_join(all_genes, by = "genes")
    
    no_orthologs    <- significant_genes %>% filter(is.na(orthologs) | orthologs == "")
    gconvert_result <- gconvert(query = no_orthologs$genes)
    
    if (!is.null(gconvert_result)) {
      gconvert_result <- gconvert_result[, c("input", "target")]
      colnames(gconvert_result) <- c("genes", "new_orthologs")
      
      significant_genes <- significant_genes %>%
        left_join(gconvert_result, by = "genes")
      
      significant_genes$orthologs <- ifelse(
        is.na(significant_genes$orthologs) | significant_genes$orthologs == "",
        significant_genes$new_orthologs, significant_genes$orthologs
      )
      significant_genes$new_orthologs <- NULL
      significant_genes <- significant_genes %>% filter(!is.na(orthologs))
    }
  } else {
    significant_genes <- ordered_genes %>%
      left_join(all_genes, by = "genes") %>%
      filter(!is.na(orthologs))
  }
  
  return(significant_genes)
}

# -----------------------------------------------------------------------------
# Slingshot Trajectory
# -----------------------------------------------------------------------------

#' Fit and Embed Slingshot Trajectory Curves
#'
#' Runs Slingshot on a SingleCellExperiment object and projects the resulting
#' principal curves into a UMAP embedding from a Seurat object.
#'
#' @param seurat_object A Seurat object containing the target embedding.
#' @param sce A SingleCellExperiment with cluster labels in sce$seurat_clusters
#'   (human) or sce$Leiden_cluster (ferret).
#' @param reduced_dim Name of the reduced dimension in sce for inference.
#' @param reduction Name of the Seurat reduction to embed curves into.
#' @param start_cluster Optional root cluster. Default NULL (auto-inferred).
#' @param end_cluster Optional terminal cluster(s). Default NULL (auto-inferred).
#' @param dataset Character; "human" or "ferret". Determines which cluster
#'   label column is used. Default "human".
#'
#' @return A named list: all_curves (embedded principal curves) and
#'   slingshot_sce (SCE with pseudotime and lineage information).
#'
#' @examples
#' curves <- get_curves(
#'   seurat_object = my_seurat, sce = my_sce,
#'   reduced_dim = "PCA", reduction = "umap",
#'   start_cluster = "0", end_cluster = "5",
#'   dataset = "human"
#' )
get_curves <- function(seurat_object,
                       sce,
                       reduced_dim,
                       reduction,
                       start_cluster = NULL,
                       end_cluster   = NULL,
                       dataset       = "human") {
  cluster_col <- if (dataset == "human") "Leiden_cluster" else "seurat_clusters"
  
  slingshot_sce <- slingshot(
    sce,
    reducedDim  = reduced_dim,
    clusterLabels = sce[[cluster_col]],
    start.clus  = start_cluster,
    end.clus    = end_cluster
  )
  
  pseudotime        <- slingPseudotime(slingshot_sce, na = FALSE)
  all_umap          <- Embeddings(seurat_object, reduction = reduction)
  sce_umap          <- all_umap[rownames(all_umap) %in% rownames(pseudotime), ]
  sce_umap_slingshot <- embedCurves(slingshot_sce, newDimRed = sce_umap)
  all_curves        <- slingCurves(sce_umap_slingshot)
  
  return(list(all_curves = all_curves, slingshot_sce = slingshot_sce))
}

# -----------------------------------------------------------------------------
# Plot Orchestration
# -----------------------------------------------------------------------------

#' Generate and display trajectory analysis plots for a specified lineage
#'
#' Produces paired UMAP and violin plots for a Slingshot lineage, organised
#' into side-by-side grids. The number of plot pairs and metadata columns used
#' depend on dataset: human produces 4 pairs, ferret produces 3.
#' Named palettes are built from the full Seurat object so colors remain
#' consistent between UMAPs and violins even when a lineage contains only a
#' subset of cell types.
#'
#' @param sce A SingleCellExperiment object for violin plot generation.
#' @param slingshot_sce A SingleCellExperiment with Slingshot pseudotime.
#' @param all_curves A named list of Slingshot curve objects.
#' @param seurat_object A Seurat object with cell metadata and reductions.
#' @param slingshot_seurat A Seurat object with pseudotime in metadata.
#' @param lineage String or index selecting the target lineage from all_curves.
#' @param dataset Character; "human" or "ferret". Controls column names,
#'   palette construction, and number of plot pairs. Default "human".
#' @param dimension Dimensional reduction for UMAP plots. Default "rna.umap".
#' @param pseudotime_limit Numeric upper cap for pseudotime display. Cells
#'   above this are greyed out. Default 150.
#' @param trim_index Optional integer; trims the trajectory curve to this many
#'   points (human only). Default NULL.
#'
#' @return A ggplot object showing the pseudotime UMAP for the lineage. Plot
#'   grids are printed to the graphics device as a side effect.
#'
#' @examples
#' # Human
#' pseudotime_plot <- getplots(
#'   sce = my_sce, slingshot_sce = my_slingshot_sce,
#'   all_curves = my_curves, seurat_object = my_seurat,
#'   slingshot_seurat = my_slingshot_seurat,
#'   lineage = "Lineage1", dataset = "human"
#' )
#'
#' # Ferret
#' pseudotime_plot <- getplots(
#'   sce = my_sce, slingshot_sce = my_slingshot_sce,
#'   all_curves = my_curves, seurat_object = my_seurat,
#'   slingshot_seurat = my_slingshot_seurat,
#'   lineage = "Lineage1", dataset = "ferret", trim_index = 80
#' )
getplots <- function(sce,
                     slingshot_sce,
                     all_curves,
                     seurat_object,
                     slingshot_seurat,
                     lineage,
                     dataset          = "human",
                     dimension        = "rna.umap",
                     pseudotime_limit = 150,
                     trim_index       = NULL) {
  
  output_plot_columns <- get_plot_columns(dataset)
  
  # Build named palettes from the full Seurat object so that colors are
  # consistent between UMAP and violin even when the lineage SCE is a subset
  named_palettes <- build_named_palettes(seurat_object, dataset, custom_palette)
  
  current_curve <- all_curves[[lineage]]
  options(warn = -1)
  
  lineage_plots <- umap_all_plots(
    curve          = current_curve,
    seurat_object  = seurat_object,
    named_palettes = named_palettes,
    column_list    = output_plot_columns,
    dimension      = dimension,
    lineage        = lineage,
    dataset        = dataset
  )
  
  violin_plots <- violin_all_plots(
    sce            = sce,
    slingshot_sce  = slingshot_sce,
    named_palettes = named_palettes,
    column_list    = output_plot_columns,
    lineage        = lineage
  )
  
  # Build and print grids — human gets 4 pairs, ferret gets 3
  grid_labels <- list(
    list("A", "B", 0.9,  1),
    list("C", "D", 0.85, 1),
    list("E", "F", 1.2,  1),
    list("G", "H", 1.2,  1)
  )
  n_pairs <- length(output_plot_columns)
  
  for (i in seq_len(n_pairs)) {
    cfg <- grid_labels[[i]]
    g <- plot_grid(
      lineage_plots[[i]], violin_plots[[i]],
      labels     = c(cfg[[1]], cfg[[2]]),
      ncol       = 2,
      rel_widths = c(cfg[[3]], cfg[[4]])
    )
    print(g)
  }
  
  options(warn = 0)
  
  lineage_pseudotime <- umap_pseudotime(
    curve            = current_curve,
    seurat_object    = slingshot_seurat,
    original_object  = seurat_object,
    slingshot_sce    = slingshot_sce,
    dimension        = dimension,
    lineage          = lineage,
    pseudotime_limit = pseudotime_limit,
    trim_index       = trim_index,
    dataset          = dataset
  )
  
  return(lineage_pseudotime)
}

# -----------------------------------------------------------------------------
# Individual Plot Functions
# -----------------------------------------------------------------------------

#' Generate multiple UMAP plots with lineage curve overlay
#'
#' Creates a list of UMAP plots from a Seurat object, coloring cells by
#' specified metadata columns and overlaying a Slingshot lineage curve.
#'
#' @param curve A Slingshot curve object containing lineage coordinates.
#' @param seurat_object Seurat object containing dimensional reduction results.
#' @param named_palettes A named list of color vectors, one per column in
#'   \code{column_list}. Each vector must be named by the full set of levels
#'   from the complete Seurat object so colors map consistently by value.
#'   Build with \code{build_named_palettes} or \code{make_named_palette}.
#' @param column_list A character vector of metadata columns to visualize.
#' @param dimension Dimensional reduction to use. Default "rna.umap".
#' @param lineage Lineage identifier (unused internally; kept for consistency).
#'
#' @return A list of ggplot objects, one per column in \code{column_list}.
#'
#' @examples
#' named_palettes <- build_named_palettes(seurat_obj, "human", custom_palette)
#' plots <- umap_all_plots(
#'   curve = sling_curves[[1]], seurat_object = seurat_obj,
#'   named_palettes = named_palettes,
#'   column_list = c("day", "majority_celltype"),
#'   dimension = "rna.umap", lineage = 1
#' )
umap_all_plots <- function(curve,
                           seurat_object,
                           named_palettes,
                           column_list,
                           dimension = "rna.umap",
                           lineage,
                           dataset = "human") {
  all_plots <- lapply(seq_along(column_list), function(x) {
    current_column <- column_list[[x]]
    cur_palette    <- named_palettes[[current_column]]
    
    curve1_coord <- data.frame(curve$s[curve$ord, c(1, 2)])
    curve1_coord$stage <- "line"
    
    base_plot <- plotDimRed(seurat_object,
                            col_by    = current_column,
                            color     = cur_palette,
                            plot_type = dimension,
                            ggrastr   = TRUE
    )[[1]]
    
    if (dimension == "rna.umap") {
      base_plot <- base_plot + ggplot2::geom_path(
        data = curve1_coord,
        ggplot2::aes(rnaUMAP_1, rnaUMAP_2),
        color = "black", size = 1
      )
    } else {
      if (dataset == "human") {
        base_plot <- base_plot + ggplot2::geom_path(
          data = curve1_coord,
          ggplot2::aes(UMAP_1, UMAP_2),
          color = "black", size = 1
        )
      } else {
        base_plot <- base_plot + ggplot2::geom_path(
          data = curve1_coord,
          ggplot2::aes(umap_1, umap_2),
          color = "black", size = 1
        )
      }
    }
    
    base_plot + theme(
      legend.text     = element_text(size = 6),
      legend.key.size = unit(0.5, "cm")
    )
  })
  
  return(all_plots)
}

#' Plot UMAP colored by pseudotime with lineage overlay
#'
#' Generates a UMAP plot where cells are colored by Slingshot pseudotime.
#' Cells above \code{pseudotime_limit} are shown in grey. The trajectory
#' curve is overlaid and can optionally be trimmed via \code{trim_index}
#' (useful for ferret datasets where the curve extends beyond the region of
#' interest).
#'
#' @param curve A Slingshot curve object containing lineage coordinates.
#' @param slingshot_sce A SingleCellExperiment with Slingshot results.
#' @param seurat_object A Seurat object used to align pseudotime values.
#' @param original_object A Seurat object to visualize (receives pseudotime).
#' @param dimension Dimensional reduction to use. Default "rna.umap".
#' @param lineage Lineage identifier.
#' @param pseudotime_limit Numeric cap for pseudotime; values above are set to
#'   NA and shown in grey. Default 150.
#' @param trim_index Optional integer; if set, the curve is trimmed to this
#'   many points before plotting. Default NULL.
#'
#' @return A cowplot grid object with the pseudotime UMAP.
#'
#' @examples
#' pt_plot <- umap_pseudotime(
#'   curve = sling_curves[[1]], slingshot_sce = sce_obj,
#'   seurat_object = seurat_obj, original_object = full_seurat_obj,
#'   dimension = "rna.umap", lineage = 1, pseudotime_limit = 120
#' )
umap_pseudotime <- function(curve,
                            slingshot_sce,
                            seurat_object,
                            original_object,
                            dimension        = "rna.umap",
                            lineage,
                            pseudotime_limit = 150,
                            trim_index       = NULL,
                            dataset = "human") {

  pseudotime_column <- paste0("slingPseudotime_", lineage)
  pseudotime_values <- colData(slingshot_sce)[[pseudotime_column]]
  names(pseudotime_values) <- colnames(seurat_object)
  
  original_object$pseudotime <- rep(NA, ncol(original_object))
  matching_cells <- intersect(colnames(seurat_object), colnames(original_object))
  original_object$pseudotime[matching_cells] <- pseudotime_values[matching_cells]
  original_object$pseudotime[original_object$pseudotime > pseudotime_limit] <- NA
  
  curve1_coord <- data.frame(curve$s[curve$ord, c(1, 2)])
  curve1_coord$stage <- "line"
  
  if (!is.null(trim_index)) {
    curve1_coord <- curve1_coord[seq_len(trim_index), ]
  }
  
  umap_coords <- Embeddings(original_object, reduction = dimension)
  umap_df <- data.frame(
    UMAP_1    = umap_coords[, 1],
    UMAP_2    = umap_coords[, 2],
    pseudotime = original_object$pseudotime
  )
  
  base_plot <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
    geom_point_rast(size = 1) +
    scale_color_gradientn(
      colors   = c("blue", "green", "yellow", "red"),
      na.value = "#D3D3D3",
      name     = "Pseudotime"
    ) +
    labs(title = "UMAP colored by Pseudotime", x = "UMAP 1", y = "UMAP 2") +
    theme(legend.position = "right")
  
  if (dimension == "rna.umap") {
    base_plot <- base_plot + ggplot2::geom_path(
      data = curve1_coord,
      ggplot2::aes(rnaUMAP_1, rnaUMAP_2),
      color = "black", size = 1
    )
  } else {
    if (dataset == "human") {
      base_plot <- base_plot + ggplot2::geom_path(
        data = curve1_coord,
        ggplot2::aes(UMAP_1, UMAP_2),
        color = "black", size = 1
      )
    } else {
      base_plot <- base_plot + ggplot2::geom_path(
        data = curve1_coord,
        ggplot2::aes(umap_1, umap_2),
        color = "black", size = 1
      )
    }
  }
  
  plot_grid(base_plot, NULL,
            labels      = c("I", ""),
            ncol        = 2,
            rel_heights = c(3, 1)
  )
}

#' Generate violin-style pseudotime distribution plots for multiple metadata
#' columns
#'
#' Creates quasi-random (violin-like) plots of pseudotime distributions grouped
#' and colored by metadata columns. Named palettes built from the full Seurat
#' object ensure that colors match the corresponding UMAP plots even when the
#' lineage SCE contains only a subset of cell types.
#'
#' @param sce A SingleCellExperiment object containing cell metadata.
#' @param slingshot_sce A SingleCellExperiment with Slingshot results.
#' @param named_palettes A named list of color vectors, one per column in
#'   \code{column_list}. Each vector must be named by the full set of levels
#'   from the complete Seurat object. Build with \code{build_named_palettes}.
#' @param column_list A character vector of metadata columns to visualize.
#' @param lineage Lineage identifier.
#'
#' @return A list of ggplot objects, one per column in \code{column_list}.
#'
#' @examples
#' named_palettes <- build_named_palettes(seurat_obj, "human", custom_palette)
#' violin_plots <- violin_all_plots(
#'   sce = sce_obj, slingshot_sce = sce_obj,
#'   named_palettes = named_palettes,
#'   column_list = c("day", "majority_celltype"),
#'   lineage = 1
#' )
violin_all_plots <- function(sce,
                             slingshot_sce,
                             named_palettes,
                             column_list,
                             lineage) {
  
  pseudotime_col <- paste0("slingPseudotime_", lineage)
  slingshot_sce$current_sling <- slingshot_sce[[pseudotime_col]]
  
  all_plots <- lapply(seq_along(column_list), function(i) {
    column      <- column_list[[i]]
    cur_palette <- named_palettes[[column]]
    
    unique_samples     <- sort(unique(sce[[column]]))
    unique_samples_rev <- rev(unique_samples)
    
    sce[[column]]    <- factor(sce[[column]], levels = unique_samples)
    sce[["reverse"]] <- factor(sce[[column]], levels = unique_samples_rev)
    
    ggplot(as.data.frame(colData(sce)), aes(
      x      = slingshot_sce$current_sling,
      y      = sce[["reverse"]],
      colour = sce[[column]]
    )) +
      geom_quasirandom(groupOnX = FALSE, orientation = "y") +
      scale_color_manual(values = cur_palette, name = column, drop = FALSE) +
      theme(
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        legend.position = "none"
      ) +
      xlab("Slingshot pseudotime")
  })
  
  return(all_plots)
}

# -----------------------------------------------------------------------------
# Heatmap
# -----------------------------------------------------------------------------

#' Generate pseudotime-ordered heatmaps for lineage-associated genes
#'
#' Filters lineage-associated differential expression results, identifies
#' significant genes by p-value and Wald statistic, and generates a heatmap
#' ordered along pseudotime. The \code{dataset} parameter controls which
#' metadata columns appear as heatmap annotations.
#'
#' @param slingshot_sce A SingleCellExperiment with Slingshot results.
#' @param object A Seurat object containing expression data.
#' @param gene_list A data frame with a \code{genes} column.
#' @param lineage Integer or character specifying the lineage.
#' @param lineage_asso_res Data frame of lineage association results from
#'   tradeSeq.
#' @param waldcutoff Wald statistic threshold. Default 200.
#' @param pval_cutoff P-value threshold. Default 0.01.
#' @param save_dir Directory for output CSV files.
#' @param fraction Numeric (0-1); fraction of pseudotime to retain. Default 1.
#' @param dataset Character; "human" or "ferret". Human includes "day" as a
#'   heatmap annotation; ferret uses pseudotime only. Default "human".
#'
#' @return A list with \code{heatmap}, \code{z_score}, and \code{counts} for
#'   significant genes ordered by pseudotime.
#'
#' @examples
#' results <- make_lineage_heatmaps(
#'   slingshot_sce = sce_obj, object = seurat_obj,
#'   gene_list = marker_df, lineage = 1,
#'   lineage_asso_res = assoc_res, waldcutoff = 200,
#'   pval_cutoff = 0.01, save_dir = "results/heatmaps",
#'   fraction = 0.8, dataset = "human"
#' )
make_lineage_heatmaps <- function(slingshot_sce,
                                  object,
                                  gene_list,
                                  lineage,
                                  lineage_asso_res,
                                  waldcutoff = 200,
                                  pval_cutoff = 0.01,
                                  save_dir    = save_dir,
                                  fraction    = 1.0,
                                  dataset     = "human") {
  keep_cols <- paste0(c("waldStat", "pvalue", "df"), "_", lineage)
  
  filtered_res <- lineage_asso_res %>%
    dplyr::select(dplyr::all_of(keep_cols)) %>%
    setNames(gsub("_[0-9]*", "", names(.))) %>%
    dplyr::filter(pvalue < pval_cutoff) %>%
    dplyr::arrange(-waldStat)
  
  # write.csv(filtered_res, file.path(save_dir, "fitgam_sig_genes.csv"))
  
  opts <- list(species = 9606)
  tfs <- getMatrixSet(JASPAR2020, opts)
  tf_gene_symbols <- sapply(tfs, function(tf) tf@name)
  filtered_res$tf_mapped <- rownames(filtered_res) %in% tf_gene_symbols
  
  filtered_genes_res <- filtered_res[rownames(filtered_res) %in% gene_list$genes, ]
  # write.csv(filtered_genes_res, file.path(save_dir, "fitgam_findmarkers_sig_genes.csv"))
  
  filtered_res_high <- filtered_res %>% dplyr::filter(waldStat > waldcutoff)
  # write.csv(filtered_res_high, file.path(save_dir, "fitgam_waldcutoff_genes.csv"))
  
  pseudotime        <- slingPseudotime(slingshot_sce, na = TRUE)[, lineage]
  pseudotime        <- as.vector(pseudotime)
  names(pseudotime) <- colnames(slingshot_sce)
  pseudotime        <- pseudotime[!is.na(pseudotime)]
  
  cells      <- names(pseudotime)
  seurat_sub <- subset(object, cells = cells)
  seurat_sub <- AddMetaData(seurat_sub, metadata = pseudotime, col.name = "pseudotime")
  
  ordered_cells        <- seurat_sub@meta.data %>% dplyr::arrange(pseudotime) %>% rownames()
  ordered_cells_subset <- ordered_cells[pseudotime[ordered_cells] <= fraction]
  seurat_sub           <- seurat_sub[, ordered_cells_subset]
  
  top_genes      <- rownames(head(filtered_res, 200))
  #top_genes      <- top_genes[top_genes %in% rownames(expr_mat)]
  expr_mat       <- GetAssayData(seurat_sub, layer = "data", assay = "RNA")
  expr_mat       <- expr_mat[top_genes, ordered_cells_subset, drop = FALSE]
  pseudotime_vals <- pseudotime[ordered_cells_subset]
  
  window_size <- 10
  n           <- length(pseudotime_vals)
  
  gene_peak <- apply(expr_mat, 1, function(x) {
    if (sum(x) == 0) return(Inf)
    window_sums <- sapply(seq_len(n - window_size + 1), function(i) sum(x[i:(i + window_size - 1)]))
    best_idx <- which.max(window_sums)
    mean(pseudotime_vals[best_idx:(best_idx + window_size - 1)])
  })
  
  gene_order         <- names(sort(gene_peak))
  seurat_sub$pseudotime <- pseudotime_vals
  
  # Human includes day annotation; ferret uses pseudotime only
  meta_col <- if (dataset == "human") c("pseudotime") else c("day", "pseudotime")
  
  sig_genes_plot <- plot_heatmap(seurat_sub,
                                 gene_list    = gene_order,
                                 meta_col     = meta_col,
                                 cluster_rows = FALSE,
                                 plot_rownames = FALSE,
                                 cell_order   = ordered_cells_subset,
                                 return_data  = TRUE
  )
  
  return(sig_genes_plot)
}

# -----------------------------------------------------------------------------
# Heatmap Helper
# -----------------------------------------------------------------------------

#' Plot a scaled gene expression heatmap
#'
#' NOTE: Adapted from scAnalysis with customizable labels and colors for use
#' with \code{make_lineage_heatmaps}.
#'
#' @param seurat_object A Seurat object.
#' @param gene_list Character vector of genes to include.
#' @param meta_col Character vector of metadata columns for annotation.
#' @param colors Optional named color vector for metadata groups.
#' @param meta_df Optional metadata data frame.
#' @param color_list Optional named list of color vectors for metadata columns.
#' @param max_val Numeric heatmap maximum. Default 2.5.
#' @param min_val Numeric heatmap minimum. Default -2.5.
#' @param cluster_rows Logical. Default FALSE.
#' @param cluster_cols Logical. Default FALSE.
#' @param average_expression Logical; average by metadata group. Default FALSE.
#' @param plot_meta_col Logical; show metadata annotation. Default TRUE.
#' @param plot_rownames Logical; show gene names. Default TRUE.
#' @param cell_order Optional character vector of cell ordering.
#' @param return_data Logical; return heatmap object and underlying data.
#'   Default FALSE.
#' @param assay Assay to use. Default "RNA".
#' @param ... Additional arguments passed to \code{pheatmap}.
plot_heatmap <- function(seurat_object,
                         gene_list,
                         meta_col,
                         colors             = NULL,
                         fontsize_row       = 5,
                         meta_df            = NULL,
                         color_list         = NULL,
                         max_val            = 2.5,
                         min_val            = -2.5,
                         cluster_rows       = FALSE,
                         cluster_cols       = FALSE,
                         average_expression = FALSE,
                         plot_meta_col      = TRUE,
                         plot_rownames      = TRUE,
                         cell_order         = NULL,
                         return_data        = FALSE,
                         assay              = "RNA", ...) {
  if (average_expression) {
    Idents(seurat_object) <- meta_col
    heatmap_df <- AverageExpression(seurat_object,
                                    seurat   = FALSE,
                                    group.by = meta_col, assays = assay
    )[[assay]]
    
    character_vals <- suppressWarnings(all(!is.na(as.numeric(as.character(colnames(heatmap_df))))))
    
    if (is.null(meta_df)) {
      sample_info <- seurat_object[[meta_col]]
      if (is.null(levels(sample_info[[meta_col]]))) {
        sample_info[[meta_col]] <- factor(sample_info[[meta_col]])
      }
      meta_df <- data.frame(levels(sample_info[[meta_col]]))
      colnames(meta_df) <- meta_col
      rownames(meta_df) <- if (character_vals) paste0("X", meta_df[[meta_col]]) else meta_df[[meta_col]]
      
      if (is.null(colors)) {
        colorcount <- length(levels(sample_info[[meta_col]]))
        colors     <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(colorcount)
        names(colors) <- levels(sample_info[[meta_col]])
      }
      color_list <- setNames(list(colors), meta_col)
    }
  } else {
    heatmap_df <- GetAssayData(seurat_object, layer = "data", assay = assay)
  }
  
  heatmap_df <- heatmap_df[rownames(heatmap_df) %in% gene_list, ]
  heatmap_df <- data.frame(heatmap_df)
  heatmap_df <- heatmap_df[order(match(rownames(heatmap_df), gene_list)), ]
  heatmap_df <- heatmap_df[rowSums(heatmap_df) > 0, ]
  
  sample_info <- if (is.null(meta_df)) seurat_object[[meta_col]] else { coloring <- color_list; meta_df }
  
  heatmap_scale <- t(scale(t(heatmap_df), scale = TRUE))
  blueYellow <- c("#352A86", "#343DAE", "#0262E0", "#1389D2",
                  "#2DB7A3", "#A5BE6A", "#F8BA43", "#F6DA23", "#F8FA0D")
  
  if (!is.null(cell_order)) {
    sample_info <- sample_info[order(match(rownames(sample_info), cell_order)), , drop = FALSE]
    rownames(sample_info) <- make.names(rownames(sample_info))
    if (!identical(colnames(heatmap_scale), rownames(sample_info))) {
      heatmap_scale <- heatmap_scale[, rownames(sample_info)]
    }
  } else if (!cluster_cols) {
    cluster_order <- levels(sample_info[[meta_col]])
    sample_info   <- sample_info[order(match(sample_info[[meta_col]], cluster_order)), , drop = FALSE]
    rownames(sample_info) <- make.names(rownames(sample_info))
    if (!identical(colnames(heatmap_scale), rownames(sample_info))) {
      heatmap_scale <- heatmap_scale[, rownames(sample_info)]
    }
  } else {
    rownames(sample_info) <- make.names(rownames(sample_info))
  }
  
  if (!plot_meta_col) sample_info[[meta_col]] <- NULL
  
  heatmap_scale <- ifelse(heatmap_scale > max_val, max_val, heatmap_scale)
  heatmap_scale <- ifelse(heatmap_scale < min_val, min_val, heatmap_scale)
  
  # Ferret datasets include a day annotation; both include pseudotime
  annotation_colors <- list(
    day        = c("22day" = "#4285F4", "24day" = "#34A853",
                   "28day" = "#FBBC05", "35day" = "#EA4335"),
    pseudotime = colorRampPalette(c("blue", "green", "yellow", "red"))(100)
  )
  
  # Map Ensembl IDs to gene symbols for row labels
  ensg_ids     <- rownames(heatmap_scale)
  mart         <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  gene_mapping <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name"),
    filters    = "ensembl_gene_id",
    values     = ensg_ids,
    mart       = mart
  )
  ens_to_gene    <- setNames(gene_mapping$external_gene_name, gene_mapping$ensembl_gene_id)
  output_mapping <- ifelse(
    rownames(heatmap_scale) %in% names(ens_to_gene) & ens_to_gene[rownames(heatmap_scale)] != "",
    ens_to_gene[rownames(heatmap_scale)],
    rownames(heatmap_scale)
  )
  rownames(heatmap_scale) <- as.character(output_mapping)
  
  heatmap <- pheatmap::pheatmap(heatmap_scale,
                                cluster_rows       = cluster_rows,
                                cluster_cols       = cluster_cols,
                                show_rownames      = plot_rownames,
                                show_colnames      = FALSE,
                                annotation_col     = sample_info,
                                annotation_colors  = annotation_colors,
                                color              = blueYellow,
                                border_color       = NA,
                                legend             = FALSE,
                                annotation_legend  = FALSE,
                                annotation_names_row = FALSE,
                                fontsize_row       = fontsize_row,
                                clustering_method  = "complete",
                                silent             = TRUE, ...
  )
  
  if (return_data) {
    return(list(heatmap = heatmap, z_score = heatmap_scale, counts = heatmap_df))
  } else {
    return(heatmap)
  }
}