## 01_pseudomultiome_3dpf.R
## Explore the 3 dpf pseudo-multiome from Fabian et al. 2022 (Nat Commun 13:13).
##
## The authors' snATAC object already contains the pseudo-multiome: its ACTIVITY
## assay is scRNA-seq expression imputed onto each ATAC nucleus (log-normalized,
## counts == data, gene-mean Spearman r = 0.94 with RNA). So no anchoring or
## transfer is done here; we analyze imputed expression in the ATAC-defined space
## and compare it with the RNA-only embedding.

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

setwd("~/Desktop/Hayden/git_repos/sox10_sc_analysis")
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
set.seed(1)

## ---- parameters -------------------------------------------------------------
marker_genes <- c(
  # positional / arch patterning (paper Fig. 3)
  "pou3f3b", "dlx5a", "hand2", "meis2b", "pitx1", "nr5a2", "hoxb3a",
  # lineages emerging around 3 dpf (paper Figs. 2, 4)
  "hyal4", "cxcl12a", "fgf10a", "fgf10b", "gata3"
)

## ---- load -------------------------------------------------------------------
rna  <- UpdateSeuratObject(readRDS("data/raw/3dpf_scRNA.rds"))
atac <- UpdateSeuratObject(readRDS("data/raw/3dpf_snATAC_SnapATAC.rds"))

DefaultAssay(rna)  <- "RNA"
DefaultAssay(atac) <- "ACTIVITY"   # imputed expression; already log-normalized

plot_genes <- intersect(marker_genes, intersect(rownames(rna), rownames(atac)))
missing <- setdiff(marker_genes, plot_genes)
if (length(missing)) message("Not found: ", paste(missing, collapse = ", "))
h <- 4 * ceiling(length(plot_genes) / 4)

## ---- clusters in each space -------------------------------------------------
p_clusters <- (DimPlot(rna, group.by = "RNA_snn_res.0.8", label = TRUE) +
                 ggtitle("scRNA-seq clusters (RNA UMAP)") + NoLegend()) |
              (DimPlot(atac, group.by = "ATAC_snn_res.0.8", label = TRUE) +
                 ggtitle("snATAC clusters (ATAC UMAP)") + NoLegend())
ggsave("results/figures/3dpf_clusters_RNA_vs_ATAC.pdf", p_clusters, width = 12, height = 5.5)

## ---- the key comparison: same genes, RNA-only vs ATAC-defined embedding -----
p_rna <- FeaturePlot(rna, plot_genes, reduction = "umap", order = TRUE, ncol = 4)
ggsave("results/figures/3dpf_markers_RNA_umap.pdf", p_rna, width = 16, height = h)

p_atac <- FeaturePlot(atac, plot_genes, reduction = "umap", order = TRUE, ncol = 4)
ggsave("results/figures/3dpf_markers_imputed_ATAC_umap.pdf", p_atac, width = 16, height = h)

## ---- markers of each ATAC cluster, from imputed expression ------------------
Idents(atac) <- "ATAC_snn_res.0.8"
atac_markers <- FindAllMarkers(atac, assay = "ACTIVITY", only.pos = TRUE,
                               min.pct = 0.25, logfc.threshold = 0.25)
write.csv(atac_markers, "results/3dpf_ATAC_cluster_markers_imputed.csv", row.names = FALSE)

## ---- save -------------------------------------------------------------------
saveRDS(atac, "data/processed/3dpf_pseudomultiome.rds")
writeLines(capture.output(sessionInfo()), "results/sessionInfo_01_pseudomultiome_3dpf.txt")
