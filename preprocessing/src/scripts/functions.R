clustifyr_orthologs <- function (seurat_object, ref_mat, save_dir,
                                 save_name, mapping_file,
                                 mapping_gene_col, mapping_ortholog_col,
                                 assay = "RNA", nfeatures = 1000,
                                 clusters = "RNA_cluster", 
                                 plot_type = "rna.umap", cor_cutoff = 0.5,
                                 features = NULL) {
  ifelse(!dir.exists(file.path(save_dir, "images")), 
         dir.create(file.path(save_dir, "images")), FALSE)
  ifelse(!dir.exists(file.path(save_dir, "files")), 
         dir.create(file.path(save_dir,  "files")), FALSE)
  
  if(is.null(features)){
    seurat_var <- FindVariableFeatures(seurat_object, assay = assay, 
                                       selection.method = "vst",
                                       nfeatures = nfeatures)    
  } else {
    seurat_var <- seurat_object
    VariableFeatures(seurat_var) <- features
  }

  
  # Keep only genes in both
  seurat_mat <- GetAssayData(object = seurat_var, slot = "data")
  seurat_mat <- seurat_mat[rownames(seurat_mat) %in%
                             mapping_file[[mapping_gene_col]],]
  
  # keep only variable genes
  seurat_mat <- seurat_mat[rownames(seurat_mat) %in% 
                             VariableFeatures(seurat_var),]
  
  if(nrow(seurat_mat) < 500){
    warning(paste0("After filtering for genes in your ortholog set only ",
                   nrow(seurat_mat), " gene remain!"))
  }
  
  # Map genes
  all_genes <- data.frame("genes" = rownames(seurat_mat)) %>%
    dplyr::filter(genes %in% mapping_file[[mapping_gene_col]])
  
  mapping_list <- mapping_file[[mapping_ortholog_col]]
  names(mapping_list) <- mapping_file[[mapping_gene_col]]
  
  all_genes$orthologs <- mapping_list[all_genes$genes]
  
  if(!identical(all_genes$genes, rownames(seurat_mat))){
    stop("Not sure how to handle cases where the row numbers don't add up")
  }
  
  # Change gene names based on the mapping
  rownames(seurat_mat) <- all_genes$orthologs
  
  ref_mat <- ref_mat[rownames(ref_mat) %in% rownames(seurat_mat), ]
  return_list <- list()
  DefaultAssay(seurat_object) <- seurat_assay
  seurat_metadata <- seurat_object[[clusters]]
  seurat_metadata[[clusters]] <- as.character(seurat_metadata[[clusters]])
  seurat_object[[clusters]] <- seurat_metadata[[clusters]]
  
  seurat_genes <- rownames(seurat_mat)
  
  seurat_res <- clustify(input = seurat_mat, metadata = seurat_metadata, 
                         ref_mat = ref_mat, query_genes = seurat_genes,
                         cluster_col = clusters)
  
  return_list$RNA <- seurat_res
  pheatmap::pheatmap(seurat_res, color = viridisLite::viridis(10))
  seurat_cluster <- cor_to_call(seurat_res) %>%
    mutate(type = ifelse(r < cor_cutoff, "undetermined", type))
  new_clusters <- seurat_cluster$type
  names(new_clusters) <- seurat_cluster$cluster
  colname <- paste0(assay, "_", save_name)
  seurat_object[[colname]] <- new_clusters[seurat_object[[clusters]][[1]]]
  plot1 <- plotDimRed(seurat_object, col_by = colname, plot_type = plot_type)
  
  pdf(file.path(save_dir, "images", paste0("RNA_mapping_", 
                                           save_name, ".pdf")))
  print(plot1)
  dev.off()
  return_list$object <- seurat_object
  return(return_list)
}

find_write_markers_orthologs <- function(seurat_object, save_dir,
                                         mapping_file = NULL,
                                         mapping_gene_col = NULL,
                                         mapping_ortholog_col = NULL,
                                         meta_col = "RNA_cluster",
                                         assay = "RNA", pval = 0.05,
                                         logfc = 0.5, gene_lists = NULL,
                                         pairwise = FALSE) {
  
  ifelse(!dir.exists(file.path(save_dir, "files", "DE")),
         dir.create(file.path(save_dir, "files", "DE")), FALSE)
  
  Idents(seurat_object) <- meta_col
  
  if (pairwise) {
    marker_genes <- find_write_markers_pairwise_ortj(seurat_object = seurat_object, 
                                                save_dir = save_dir,
                                                meta_col = meta_col,
                                                assay = assay, 
                                                pval = pval, logfc = logfc,
                                                gene_lists = gene_lists,
                                                mapping_file = mapping_file,
                                                mapping_gene_col = mapping_gene_col,
                                                mapping_ortholog_col = mapping_ortholog_col)
  } else {
    marker_genes <- find_write_markers_orig_orth(seurat_object = seurat_object, 
                                            save_dir = save_dir,
                                            meta_col = meta_col, assay = assay, 
                                            pval = pval, logfc = logfc,
                                            gene_lists = gene_lists,
                                            mapping_file = mapping_file,
                                            mapping_gene_col = mapping_gene_col,
                                            mapping_ortholog_col = mapping_ortholog_col)
  }
  return(marker_genes)
}

find_write_markers_orig_orth <- function(seurat_object, save_dir,
                                         meta_col = "RNA_cluster", 
                                         mapping_file = NULL,
                                         mapping_gene_col = NULL,
                                         mapping_ortholog_col = NULL,
                                         assay = "RNA", pval = 0.05,
                                         logfc = 0.5, gene_lists = NULL) 
{
  marker_genes <- FindAllMarkers(seurat_object, assay = assay, 
                                 only.pos = TRUE)
  if(!is.null(mapping_file)){
    # Add in ortholog gene
    new_mapping <- mapping_file %>%
      dplyr::select(all_of(c(mapping_gene_col, mapping_ortholog_col)))
    
    marker_genes <- merge(marker_genes, new_mapping, by.x = "gene",
                          by.y = mapping_gene_col, all.x = TRUE,
                          all.y = FALSE) %>%
      distinct() %>%
      group_by(cluster) %>%
      arrange(p_val_adj, .by_group = TRUE)
  }
  
  write.csv(marker_genes, file = file.path(save_dir, "files", 
                                           "DE", paste0(assay, "_markers_",
                                                        meta_col, ".csv")))
  gene_wb <- createWorkbook()
  
  full_list <- lapply(unique(marker_genes$cluster), function(x) {
    x <- as.character(x)
    new_df <- marker_genes %>% dplyr::filter(cluster == x & 
                                               p_val_adj < pval &
                                               avg_log2FC > logfc)
    addWorksheet(gene_wb, x)
    writeData(gene_wb, x, new_df)
  })
  
  if (!is.null(gene_lists)) {
    hypergeometric <- hypergeometric_test(seurat_object = seurat_data, 
                                          gene_list = gene_lists,
                                          DE_table = marker_genes, 
                                          DE_p_cutoff = 0.05,
                                          DE_lfc_cutoff = 0.5,
                                          correction_method = "fdr")
    write.csv(hypergeometric, file = file.path(save_dir, 
                                               "files", "DE",
                                               paste0(assay, "_hypergeometric_", 
                                                      meta_col, ".csv")))
    
    full_list <- lapply(unique(hypergeometric$cluster), function(x) {
      x <- as.character(x)
      new_df <- hypergeometric %>% dplyr::filter(cluster == 
                                                   x)
      worksheet_name <- paste0(x, "_gse")
      addWorksheet(gene_wb, worksheet_name)
      writeData(gene_wb, worksheet_name, new_df)
    })
  }
  saveWorkbook(gene_wb, file = file.path(save_dir, "files", 
                                         "DE", paste0(assay, "_markers_",
                                                      meta_col, ".xlsx")), 
               overwrite = TRUE)
  return(marker_genes)
}

find_write_markers_pairwise_ortj <- function(seurat_object, save_dir,
                                             mapping_file = NULL,
                                             mapping_gene_col = NULL,
                                             mapping_ortholog_col = NULL,
                                             meta_col = "RNA_cluster",
                                             assay = "RNA", pval = 0.05,
                                             logfc = 0.5, gene_lists = NULL) {
  marker_genes <- pairwise_markers(seurat_object, assay = seurat_assay, 
                                   meta_col = meta_col)
  
  if(!is.null(mapping_file)){
    # Add in ortholog gene
    new_mapping <- mapping_file %>%
      dplyr::select(all_of(c(mapping_gene_col, mapping_ortholog_col)))
    
    marker_genes <- merge(marker_genes, new_mapping, by.x = "gene",
                          by.y = mapping_gene_col, all.x = TRUE,
                          all.y = FALSE) %>%
      distinct() %>%
      group_by(cluster) %>%
      arrange(p_val_adj, .by_group = TRUE)
  }
  
  write.csv(marker_genes, file = file.path(save_dir, "files", 
                                           "DE",
                                           paste0(assay,
                                                  "_pairwise_markers_",
                                                  meta_col, ".csv")))
  
  gene_wb <- createWorkbook()
  values <- unique(c(marker_genes$cluster_down, marker_genes$cluster_up))
  combinations <- combn(unique(c(marker_genes$cluster_down, 
                                 marker_genes$cluster_up)), m = 2)
  full_list <- lapply(1:ncol(combinations), function(x) {
    ident1 <- combinations[1, x]
    ident2 <- combinations[2, x]
    sheet_name <- paste0(ident1, "_vs_", ident2)
    new_df <- marker_genes %>% 
      dplyr::filter((cluster_up == ident1 & 
                       cluster_down == ident2) | 
                      (cluster_up == ident2 & cluster_down == ident1)) %>% 
      dplyr::filter(p_val_adj < pval & avg_log2FC > logfc)
    addWorksheet(gene_wb, sheet_name)
    writeData(gene_wb, sheet_name, new_df)
  })
  marker_genes$cluster <- paste0(marker_genes$cluster_up, "_vs_", 
                                 marker_genes$cluster_down)
  if (!is.null(gene_lists)) {
    hypergeometric <- hypergeometric_test(seurat_object = seurat_data, 
                                          gene_list = gene_lists,
                                          DE_table = marker_genes, 
                                          DE_p_cutoff = 0.05,
                                          DE_lfc_cutoff = 0.5, correction_method = "fdr")
    write.csv(hypergeometric, file = file.path(save_dir, 
                                               "files", "DE", paste0(assay,
                                                                     "_hypergeometric_", 
                                                                     meta_col,
                                                                     ".csv")))
    full_list <- lapply(unique(hypergeometric$cluster), function(x) {
      x <- as.character(x)
      new_df <- hypergeometric %>% dplyr::filter(cluster == 
                                                   x)
      worksheet_name <- paste0(x, "_gse")
      addWorksheet(gene_wb, worksheet_name)
      writeData(gene_wb, worksheet_name, new_df)
    })
  }
  saveWorkbook(gene_wb, file = file.path(save_dir, "files", 
                                         "DE", paste0(assay,
                                                      "_pairwise_markers_",
                                                      meta_col, ".xlsx")), 
               overwrite = TRUE)
  
  return(marker_genes)
}

make_plots <- function(seurat_object, cluster_name, celltype_name,
                       cluster_colors, celltype_colors, gene,
                       assay = "RNA", plot_type = "rna.umap",
                       save_name = NULL){
  
  violin1 <- featDistPlot(seurat_object, gene, sep_by = celltype_name,
                          col_by = celltype_name, color = celltype_colors,
                          assay = assay)
  
  
  umap1 <- plotDimRed(seurat_object, col_by = gene, plot_type = plot_type,
                      assay = assay)[[1]]
  
  umap3 <- plotDimRed(seurat_object, col_by = celltype_name,
                      plot_type = plot_type, color = celltype_colors)[[1]]
  
  umap2 <- NULL
  umap4 <- plotDimRed(seurat_object, col_by = cluster_name,
                      plot_type = plot_type, color = cluster_colors)[[1]]
  violin2 <- featDistPlot(seurat_object, gene, sep_by = cluster_name,
                          col_by = cluster_name, color = cluster_colors,
                          assay = assay)
  text_labels <- c("A", "", "B", "C", "D", "E")
  
  
  full_plot <- plot_grid(umap1, umap2,
                         umap3, violin1,
                         umap4, violin2,
                         labels = text_labels,
                         nrow = 3,
                         ncol = 2,
                         align = "hv",
                         axis = "l")
  
  if(!is.null(save_name)){
    pdf(save_name, height = 10, width = 10)
    print(full_plot)
    dev.off()
  }
  return(full_plot)
  
}

# Pseudotime -------------------------------------------------------------------
# Seurat object = seurat object
# lineages = list of lineages that are in the metadata of the seurat object
# Gene = gene to plot across pseudotime
# col_by = color for the plot (can be "lineage" to color by lineage,
# "pseudotime" to color by the pseudotime vals or
# a value from the metadata)
# color = custom colors to pass
# Max pseudotime = for aesthetics - should a few cells that are outliers 
# be removed? This should be improved
plotPseudotime <- function(seurat_object, lineages, gene_list,
                           col_by = "lineage", color = NULL,
                           save_plot = NULL, max_pseudotime = NULL, ...){
  
  plot_list <- lapply(gene_list, function(x) {
    print(x)
    plotPseudotimeSingle(seurat_object = seurat_object,
                         lineages = lineages, gene = x,
                         col_by = col_by, color = color,
                         max_pseudotime = max_pseudotime, ...)
  })
  if (!is.null(save_plot)){
    pdf(save_plot)
    print(plot_list)
    dev.off()
  }
  return(plot_list)
}

plotPseudotimeSingle <- function(seurat_object, lineages, gene,
                                 col_by = "lineage", color = NULL,
                                 max_pseudotime = NULL,
                                 line_color = "black", alpha = 0.75,
                                 ...){
  
  # Get data
  # If col_by is lineage or pseudotime, do this,
  # otherwise add col_by to the call
  if (col_by %in% c("pseudotime", "lineage")){
    data_info <- Seurat::FetchData(seurat_object, vars = c(lineages, gene))
  } else {
    data_info <- Seurat::FetchData(seurat_object,
                                   vars = c(lineages, gene, col_by))
  }
  
  # Here I'll need to pivot_longer by the lineages
  data_info <- data_info %>%
    tidyr::pivot_longer(cols = all_of(lineages), names_to = "lineage",
                        values_to = "pseudotime") %>%
    dplyr::rename(gene = dplyr::all_of(gene)) %>%
    dplyr::filter(!is.na(pseudotime))
  
  data_info$col_by <- data_info[[col_by]]
  
  if(!is.null(max_pseudotime)){
    data_info <- data_info %>%
      dplyr::filter(pseudotime < max_pseudotime)
  }
  
  # Make plot
  base_plot <- ggplot2::ggplot(data = data_info, ggplot2::aes(x = pseudotime,
                                                              y = gene,
                                                              color = col_by,
                                                              group = lineage)) +
    ggplot2::geom_point(alpha = alpha, size = 0.25) +
    ggplot2::ylab(gene)
  
  if(line_color == "lineage"){
    base_plot <- base_plot +
      ggplot2::geom_smooth(method = "gam")
  } else {
    base_plot <- base_plot + 
      ggplot2::geom_smooth(method = "gam", color = line_color)
  }
  
  # Next figure out colors, if col_by is continuous, default is 
  # viridis. If it's not continuous do RColorBrewer as I like
  if (is.numeric(data_info$col_by)){
    if (is.null(color)) {
      base_plot <- base_plot +
        ggplot2::scale_color_viridis_c(option = "magma") +
        labs(color = col_by)
    } else {
      low <- color[1]
      high <- color[2]
      base_plot <- base_plot + 
        ggplot2::scale_color_gradient(low = low, high = high, name = col_by)
    }
  } else {
    # Make colors if they don't already exist
    if (is.null(color)){
      # Factor in case it isn't already
      if(!is.factor(data_info$col_by)){
        data_info$col_by <- factor(data_info$col_by)
      }
      nColors <- length(levels(data_info$col_by))
      color <- grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(9, "Set1"))(nColors)
      names(color) <- levels(data_info$col_by)
    }
    base_plot <- base_plot + scale_color_manual(values = color,
                                                name = col_by)
  }
  
  return(base_plot)
  
}


make_lineage_heatmaps <- function(seurat_object, lineage_name,
                                  lineage_number, lineage_asso_res,
                                  image_dir, sample_colors,
                                  waldcutoff = 200, pval_cutoff = 0.05){
  
  keep_cols <- c("waldStat", "pvalue", "df")
  keep_cols <- paste0(keep_cols, "_", lineage_number)
  
  # Filter by pvalue
  filtered_res <- lineage_asso_res %>%
    dplyr::select(dplyr::all_of(keep_cols)) %>%
    setNames(gsub("_[0-9]*", "", names(.))) %>%
    dplyr::filter(pvalue < pval_cutoff)
  
  # Filter by wald cutoff (fewer genes)
  filtered_res_high <- filtered_res %>%
    dplyr::filter(waldStat > waldcutoff)
  
  # Subset to just the lineage
  cells <- seurat_object[[]] %>%
    dplyr::select(dplyr::all_of(lineage_name)) %>%
    tidyr::drop_na() %>%
    rownames()
  
  # Keep only cells in the lineage
  seurat_sub <- subset(seurat_object, cells = cells)
  
  # Order the cells by pseudotime
  order_cells <- seurat_sub[[]] %>%
    dplyr::select(all_of(lineage_name)) %>%
    dplyr::arrange_all()
  
  # Make plot of all sig genes
  pdf(file.path(image_dir, paste0(lineage_name, "_heatmap.pdf")))
  print(plot_heatmap(seurat_sub, gene_list = rownames(filtered_res),
                     meta_col = "sample", colors = sample_colors,
                     cell_order = rownames(order_cells),
                     cluster_rows = TRUE, plot_rownames = FALSE))
  dev.off()
  
  
  
  # Make plot of all genes with high wald
  pdf(file.path(image_dir,
                paste0(lineage_name, "_heatmap_wald", waldcutoff, ".pdf")))
  
  print(plot_heatmap(seurat_sub, gene_list = rownames(filtered_res_high),
                     meta_col = "sample", colors = sample_colors,
                     cell_order = rownames(order_cells),
                     cluster_rows = TRUE, plot_rownames = FALSE))
  dev.off()
  
  
  # Make plot of all genes with high wald as png
  png(file.path(image_dir,
                paste0(lineage_name, "_heatmap_wald", waldcutoff, ".png")))
  
  print(plot_heatmap(seurat_sub, gene_list = rownames(filtered_res_high),
                     meta_col = "sample", colors = sample_colors,
                     cell_order = rownames(order_cells),
                     cluster_rows = TRUE, plot_rownames = FALSE))
  dev.off()
  
  
  # Make a heatmap with rownames, return so that the size can be altered
  return_heatmap <- plot_heatmap(seurat_sub, 
                                 gene_list = rownames(filtered_res_high),
                                 meta_col = "sample", colors = sample_colors,
                                 cell_order = rownames(order_cells),
                                 cluster_rows = TRUE, plot_rownames = TRUE)
  
  return(list(sig_genes = filtered_res,
              sig_high_genes = filtered_res_high,
              heatmap = return_heatmap))
  
}

plot_overlapping_lineages <- function(seurat_object,
                                      cfko_lineage_genes,
                                      wt_lineage_genes,
                                      cfko_lineage_name,
                                      wt_lineage_name,
                                      save_dir = NULL,
                                      colors = NULL){
  
  if(is.null(colors)){
    lineage_colors <- met.brewer("Egypt", 2)
    names(lineage_colors) <- c(wt_lineage_name, cfko_lineage_name)
    
  }
  cfko_lineage_df <- data.frame("cfko_cluster" = cfko_lineage_genes,
                                "gene" = names(cfko_lineage_genes))
  
  wt_lineage_df <- data.frame("wt_cluster" = wt_lineage_genes,
                              "gene" = names(wt_lineage_genes))
  
  combined_lineage <- merge(cfko_lineage_df, wt_lineage_df,
                            by = "gene", all.x = FALSE, all.y = FALSE)
  
  
  if(!is.null(save_dir)){
    all_plots <- plotPseudotime(seurat_object,
                                lineages = c(cfko_lineage_name,
                                             wt_lineage_name),
                                gene_list = combined_lineage$gene,
                                col_by = "lineage", color = lineage_colors,
                                line_color = "lineage",
                                alpha = 0.25)
    
    
    pdf(file.path(save_dir,
                  paste0(wt_lineage_name, "_and_", cfko_lineage_name, ".pdf")),
        width = 4, height = 4)
    
    print(all_plots)
    
    dev.off()    
  }
  
  
  return(combined_lineage)
  
}

