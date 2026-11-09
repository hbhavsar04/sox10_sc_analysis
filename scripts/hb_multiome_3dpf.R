library(Seurat)
library(Signac)
library(Matrix)
library(ggplot2)
library(patchwork)

setwd("~/Desktop/Hayden/git_repos/sox10_sc_analysis")

rna  <- UpdateSeuratObject(readRDS("data/raw/3dpf_scRNA.rds"))
atac <- UpdateSeuratObject(readRDS("data/raw/3dpf_snATAC_SnapATAC.rds"))

marker_genes <- c(
  # positional / arch patterning (paper Fig. 3)
  "pou3f3b", "dlx5a", "hand2", "meis2b", "pitx1", "nr5a2", "hoxb3a",
  # lineages emerging around 3 dpf (paper Figs. 2, 4)
  "hyal4", "cxcl12a", "fgf10a", "fgf10b", "gata3"
)

# check if imputed
#act_counts <- GetAssayData(atac, assay = "ACTIVITY", layer = "counts")
#vals <- if (inherits(act_counts, "dgCMatrix")) act_counts@x else as.vector(act_counts)
#vals <- sample(vals, min(length(vals), 1e6))
#cat("ACTIVITY:", nrow(act_counts), "genes x", ncol(act_counts), "cells\n")
#cat("Integer-valued:", all(vals == round(vals)), "\n")
#cat("Genes shared with RNA:", length(intersect(rownames(act_counts), rownames(rna))), "\n")

#act_data <- GetAssayData(atac, assay = "ACTIVITY", layer = "data")
#rna_data <- GetAssayData(rna, assay = "RNA", layer = "data")

#cat("fraction nonzero:", if (inherits(act_counts, "dgCMatrix"))
  #length(act_counts@x) / prod(dim(act_counts)) else mean(act_counts != 0), "\n")
#print(summary(vals))
#cat("counts == data:", identical(act_counts, act_data), "\n")
#cat("gene-mean correlation with RNA (Spearman):",
    #cor(Matrix::rowMeans(act_counts), Matrix::rowMeans(rna_data[rownames(act_counts), ]),
        #method = "spearman"), "\n")

DefaultAssay(rna)  <- "RNA"
DefaultAssay(atac) <- "ACTIVITY"   # imputed expression; already log

plot_genes <- intersect(marker_genes, intersect(rownames(rna), rownames(atac)))
missing <- setdiff(marker_genes, plot_genes)
if (length(missing)) message("Not found: ", paste(missing, collapse = ", "))
h <- 4 * ceiling(length(plot_genes) / 4)

# clusters in each space
p_clusters <- (DimPlot(rna, group.by = "RNA_snn_res.0.8", label = TRUE) +
                 ggtitle("scRNA-seq clusters (RNA UMAP)") + NoLegend()) |
  (DimPlot(atac, group.by = "ATAC_snn_res.0.8", label = TRUE) +
     ggtitle("snATAC clusters (ATAC UMAP)") + NoLegend())
ggsave("results/figures/3dpf_clusters_RNA_vs_ATAC.pdf", p_clusters, width = 12, height = 5.5)

## the key comparison: same genes, RNA-only vs ATAC-defined embedding
p_rna <- FeaturePlot(rna, plot_genes, reduction = "umap", order = TRUE, ncol = 4)
ggsave("results/figures/3dpf_markers_RNA_umap.pdf", p_rna, width = 16, height = h)
p_atac <- FeaturePlot(atac, plot_genes, reduction = "umap", order = TRUE, ncol = 4)
ggsave("results/figures/3dpf_markers_imputed_ATAC_umap.pdf", p_atac, width = 16, height = h)

# ATAC_markers
Idents(atac) <- "ATAC_snn_res.0.8"
atac_markers <- FindAllMarkers(atac, assay = "ACTIVITY", only.pos = TRUE,
                               min.pct = 0.25, logfc.threshold = 0.25)
write.csv(atac_markers, "results/3dpf_ATAC_cluster_markers_imputed.csv", row.names = FALSE)
# save results
saveRDS(atac, "data/processed/3dpf_pseudomultiome.rds")
writeLines(capture.output(sessionInfo()), "results/sessionInfo_01_pseudomultiome_3dpf.txt")