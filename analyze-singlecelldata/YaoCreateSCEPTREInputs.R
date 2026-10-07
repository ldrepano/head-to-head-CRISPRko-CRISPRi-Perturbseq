# Author: Erika Li
# 10/2/2026

library("Seurat")
library("ggplot2")
library("dplyr")
library("patchwork")
library("Matrix")
library("stringr")
library("tidyr")

setwd("~/CRISPRko_CRISPRi_Benchmark_PerturbSeq_Manuscript")

# helper: 3-panel QC histogram (UMIs, genes, % mito)
make_qc_plot <- function(obj, title, show_cutoffs = FALSE) {
  plot_data <- FetchData(obj, vars = c("Total_unique_genes", "Total_RNA_count", "Percent_mitochondrial_reads")) %>%
    rename("# UMIS per Cell" = "Total_RNA_count",
           "# Genes per Cell" = "Total_unique_genes",
           "% mitochondrial transcripts" = "Percent_mitochondrial_reads")
  
  plot_list <- list()
  for (col_name in colnames(plot_data)) {
    col_mean <- mean(plot_data[[col_name]], na.rm = TRUE)
    col_sd   <- sd(plot_data[[col_name]], na.rm = TRUE)
    
    p <- ggplot(plot_data, aes(x = as.numeric(.data[[col_name]]))) +
      geom_histogram(bins = 20, color = "white", alpha = 0.7) +
      theme_minimal() + labs(x = col_name, y = "count") +
      annotate(x = Inf, y = Inf, geom = "text",
               label = paste0("mean = ", round(col_mean, digits = 2)),
               color = "red", size = 3, hjust = 1.1, vjust = 1.5)
    
    if (show_cutoffs) {
      xint <- switch(col_name,
                     "% mitochondrial transcripts" = 15,
                     "# Genes per Cell"            = col_mean - 2 * col_sd,
                     "# UMIS per Cell"             = col_mean - col_sd)
      p <- p + geom_vline(xintercept = xint, color = "darkred", linetype = "dashed", linewidth = 1)
    }
    plot_list[[col_name]] <- p
  }
  wrap_plots(plot_list, ncol = 3) + plot_annotation(title = title)
}

main <- function(samplename, rds_filename) {
  
  if (!dir.exists(samplename)) dir.create(samplename, recursive = TRUE)
  
  # load Data
  seurat_obj    <- readRDS(rds_filename)
  
  # mitochondrial transcript filter (< 15%)
  mito_filtered <- subset(seurat_obj, subset = Percent_mitochondrial_reads < 15)
  
  # MOI plot
  m <- GetAssayData(mito_filtered, assay = "perturbations", layer = "counts")
  nonzero_per_cell <- Matrix::colSums(m > 0)
  
  df_moi <- data.frame(num_guides = factor(nonzero_per_cell))
  p_moi <- ggplot(df_moi, aes(x = num_guides)) +
    geom_bar(fill = "steelblue") +
    geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
    labs(title = paste(samplename, "MOI - Yao et al."),
         x = "Number of Detected Guides per Cell",
         y = "Cell Count") +
    theme_minimal()
  ggsave(file.path(samplename, paste0(samplename, "-moi-plot.png")), plot = p_moi)
  
  # filter to single guides
  keep_single <- names(nonzero_per_cell)[nonzero_per_cell == 1]
  m_single    <- m[, keep_single]
  
  assigned_target <- rownames(m_single)[apply(m_single, 2, which.max)]
  names(assigned_target) <- keep_single
  
  full_guide_id <- mito_filtered$Guides[keep_single]
  names(full_guide_id) <- keep_single
  
  guide_gene_prefix <- sub("_.*", "", full_guide_id)
  guide_is_clean    <- !grepl("--", full_guide_id) & !is.na(full_guide_id) & full_guide_id != ""
  agrees            <- guide_gene_prefix == assigned_target
  
  keep_consistent <- keep_single[guide_is_clean & agrees]
  cat(length(keep_single) - length(keep_consistent),
      "cells dropped due to disagreement or missing/ambiguous Guides value\n")
  
  assigned_guide <- full_guide_id[keep_consistent]
  names(assigned_guide) <- keep_consistent
  
  guide_filtered <- subset(mito_filtered, cells = keep_consistent)
  
  # QC plot of guide-filtered cells (before quality filter) 
  qc_combined <- make_qc_plot(guide_filtered,
                              paste("Yao et al.,", samplename, "Screen"),
                              show_cutoffs = TRUE)
  ggsave(file.path(samplename, "combined_qc_plot.png"), plot = qc_combined,
         width = 20, height = 5, dpi = 300)
  
  # quality filter (UMIs / genes), thresholds from guide-filtered cells 
  mean_rna   <- mean(guide_filtered$Total_RNA_count)
  sd_rna     <- sd(guide_filtered$Total_RNA_count)
  mean_genes <- mean(guide_filtered$Total_unique_genes)
  sd_genes   <- sd(guide_filtered$Total_unique_genes)
  
  pass_qc <- guide_filtered$Total_RNA_count    > (mean_rna - sd_rna) &
    guide_filtered$Total_unique_genes > (mean_genes - 2 * sd_genes)
  qc_cells <- colnames(guide_filtered)[pass_qc]
  cat(sum(!pass_qc), "cells dropped by UMI/gene quality filter\n")
  
  qc_filtered <- subset(guide_filtered, cells = qc_cells)
  
  qc_combined_f <- make_qc_plot(qc_filtered,
                                paste("Yao et al.,", samplename, "Screen - QC Filtered"))
  ggsave(file.path(samplename, "combined_filtered_qc_plot.png"), plot = qc_combined_f,
         width = 20, height = 5, dpi = 300)
  
  # min cells per guide (>= 10)
  assigned_guide_qc <- assigned_guide[qc_cells]
  guide_tab   <- table(assigned_guide_qc)
  keep_guides <- names(guide_tab)[guide_tab >= 10]
  
  final_cells      <- names(assigned_guide_qc)[assigned_guide_qc %in% keep_guides]
  final_guide_call <- assigned_guide_qc[final_cells]
  
  filtered_obj <- subset(qc_filtered, cells = final_cells)
  filtered_obj$assigned_guide <- final_guide_call[colnames(filtered_obj)]
  stopifnot(!any(is.na(filtered_obj$assigned_guide)))
  
  # cell coverage plot (post-QC counts per guide)
  df_guides <- data.frame(guide = names(guide_tab), count = as.numeric(guide_tab))
  p_cov <- ggplot(df_guides, aes(x = count)) +
    geom_histogram(binwidth = 10, fill = "darkseagreen4", color = "white") +
    geom_vline(xintercept = 10, color = "firebrick", linetype = "dashed", linewidth = 1) +
    labs(title = paste(samplename, "Cell Counts per Guide - Yao et al."),
         subtitle = "Dashed line indicates 10-cell cutoff",
         x = "Cells per Guide", y = "Number of Guides") +
    coord_cartesian(xlim = c(0, 300)) +
    theme_minimal()
  ggsave(file.path(samplename, paste0(samplename, "-cell-coverage-plot.png")), plot = p_cov)
  
  # assemble SCEPTRE input
  guide_ids_f  <- filtered_obj$assigned_guide
  cells_f      <- colnames(filtered_obj)
  guide_factor <- factor(guide_ids_f)
  
  grna_matrix <- sparseMatrix(
    i = as.integer(guide_factor),
    j = seq_along(cells_f),
    x = 1,
    dimnames = list(levels(guide_factor), cells_f)
  )
  response_matrix <- filtered_obj[["RNA"]]$counts
  
  grna_ids <- rownames(grna_matrix)
  is_non_targeting <- grepl("^non-targeting", grna_ids, ignore.case = TRUE)
  
  # drop non-targeting rows and use safe-targeting as controls
  grna_matrix <- grna_matrix[!is_non_targeting, , drop = FALSE]
  grna_ids    <- rownames(grna_matrix)
  
  is_control   <- grepl("^safe-targeting", grna_ids, ignore.case = TRUE)
  gene_targets <- sub("_.+", "", grna_ids, ignore.case = TRUE)
  
  grna_target_data_frame <- data.frame(
    grna_id     = grna_ids,
    grna_target = ifelse(is_control, "non-targeting", gene_targets)
  )
  
  # drop cells whose only guide was non-targeting
  had_nt_only     <- Matrix::colSums(grna_matrix) == 0
  grna_matrix     <- grna_matrix[, !had_nt_only, drop = FALSE]
  response_matrix <- response_matrix[, !had_nt_only, drop = FALSE]
  
  # save SCEPTRE inputs
  saveRDS(response_matrix, file.path(samplename, "response_matrix.RDS"))
  saveRDS(grna_target_data_frame, file.path(samplename, "grna_target_data_frame.RDS"))
  saveRDS(grna_matrix, file.path(samplename, "grna_matrix.RDS"))
}

main("CRISPRko", "GSM6858447_KO_conventional.rds")
main("CRISPRi", "GSM6858449_KD_conventional.rds")