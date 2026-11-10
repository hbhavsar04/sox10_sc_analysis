library(Seurat); library(Signac); library(Matrix)
library(dplyr); library(tidyr); library(ggplot2)

aob <- readRDS("data/processed/5dpf_ATAC_OB_reclustered.rds")

imputed_heatmap <- function(obj, genes, group = "sub_chrom", scale_rows = TRUE) {
  dat <- GetAssayData(obj, assay = "ACTIVITY", layer = "data")
  found <- intersect(genes, rownames(dat))
  missing <- setdiff(genes, found)
  if (length(missing)) message("Not in imputed data: ", paste(missing, collapse = ", "))
  
  grp <- as.character(obj@meta.data[[group]])
  m <- sapply(split(seq_len(ncol(obj)), grp),
              function(i) Matrix::rowMeans(dat[found, i, drop = FALSE]))
  if (is.null(dim(m))) m <- matrix(m, nrow = 1, dimnames = list(found, names(m)))
  m <- m[, order(as.numeric(colnames(m))), drop = FALSE]
  
  vals <- if (scale_rows) t(scale(t(m))) else m
  df <- as.data.frame(vals) |>
    mutate(gene = rownames(vals)) |>
    pivot_longer(-gene, names_to = "subcluster", values_to = "value") |>
    mutate(raw = as.vector(t(m))[match(paste(gene, subcluster),
                                       paste(rep(rownames(m), each = ncol(m)),
                                             rep(colnames(m), nrow(m))))],
           gene = factor(gene, levels = rev(found)),
           subcluster = factor(subcluster, levels = colnames(m)))
  
  ggplot(df, aes(subcluster, gene, fill = value)) +
    geom_tile(colour = "grey90") +
    geom_text(aes(label = sprintf("%.2f", raw)), size = 2.8) +
    { if (scale_rows) scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", name = "z")
      else scale_fill_viridis_c(name = "mean") } +
    theme_bw() + labs(x = paste("Subcluster (", group, ")"), y = NULL,
                      title = "Mean imputed expression by subcluster",
                      subtitle = if (scale_rows) "color: z-score across subclusters; numbers: raw mean" else NULL)
}

genes <- c("ihha", "ptch1", "ptch2", "gli1",              # Hh ligand + targets
           "col10a1a", "ifitm5", "panx3", "slc8a4b",      # mature OB
           "sp7", "runx2a", "runx2b", "nfatc1", "abcc6a", # landmarks
           "scxa", "lum", "mafa", "matn4", "sox6")        # tendon/perichondrium-like

imputed_heatmap(aob, genes, group = "sub_chrom")   # chromatin subclusters
imputed_heatmap(aob, genes, group = "sub_imp")     # imputed-expression subclusters

# where does ihha sit on the OB chromatin UMAP?
DefaultAssay(aob) <- "ACTIVITY"
FeaturePlot(aob, "ihha", reduction = "umap_ob", order = TRUE) |
  DimPlot(aob, reduction = "umap_ob", group.by = "sub_chrom", label = TRUE) + NoLegend()