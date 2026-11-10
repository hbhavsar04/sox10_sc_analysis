## 01_pseudomultiome_3dpf.R
## Pseudo-multiome integration of 3 dpf Sox10>dsRed scRNA-seq and snATAC-seq,
## following Fabian et al. 2022 (Nat Commun 13:13), Methods:
##   ATAC gene-activity scores -> anchors against scRNA-seq -> RNA expression
##   transferred onto every ATAC nucleus, weighted by the ATAC LSI embedding.
## Result: each ATAC cell gets an imputed expression profile, analyzed in ATAC space.

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
lsi_dims    <- 2:30   # LSI 1 mostly tracks sequencing depth, so it's dropped
anchor_dims <- 1:30
min_cells   <- 10     # impute genes detected in >= this many RNA cells
marker_genes <- c(
  # positional / arch patterning (paper Fig. 3)
  "pou3f3b", "dlx5a", "hand2", "meis2b", "pitx1", "nr5a2", "hoxb3a",
  # lineages emerging around 3 dpf (paper Figs. 2, 4)
  "hyal4", "cxcl12a", "fgf10a", "fgf10b", "gata3"
)

## ---- load -------------------------------------------------------------------
rna  <- UpdateSeuratObject(readRDS("data/raw/3dpf_scRNA.rds"))
atac <- UpdateSeuratObject(readRDS("data/raw/3dpf_snATAC_SnapATAC.rds"))

## ---- what is in the ACTIVITY assay? -----------------------------------------
## If these are integer counts, it's the raw SnapATAC gene-body activity (gmat).
## If non-integer, the authors may have already stored imputed expression here.
act_counts <- GetAssayData(atac, assay = "ACTIVITY", layer = "counts")
vals <- if (inherits(act_counts, "dgCMatrix")) act_counts@x else as.vector(act_counts)
vals <- sample(vals, min(length(vals), 1e6))
cat("ACTIVITY:", nrow(act_counts), "genes x", ncol(act_counts), "cells\n")
cat("Integer-valued:", all(vals == round(vals)), "\n")
cat("Genes shared with RNA:", length(intersect(rownames(act_counts), rownames(rna))), "\n")

## ---- prepare RNA reference --------------------------------------------------
DefaultAssay(rna) <- "RNA"
if (length(VariableFeatures(rna)) == 0) {
  rna <- FindVariableFeatures(rna, nfeatures = 2000)
}

## ---- prepare ATAC gene activity as the query --------------------------------
DefaultAssay(atac) <- "ACTIVITY"
atac <- NormalizeData(atac, assay = "ACTIVITY",
                      normalization.method = "LogNormalize",
                      scale.factor = median(Matrix::colSums(act_counts)))

genes <- intersect(VariableFeatures(rna), rownames(act_counts))
genes <- genes[Matrix::rowSums(act_counts[genes, , drop = FALSE] > 0) > 0]
cat("Anchor features:", length(genes), "\n")
atac <- ScaleData(atac, assay = "ACTIVITY", features = genes)

## ---- anchors ----------------------------------------------------------------
anchors <- FindTransferAnchors(
  reference = rna, query = atac, features = genes,
  reference.assay = "RNA", query.assay = "ACTIVITY",
  reduction = "cca", dims = anchor_dims
)

## ---- sanity check: transfer RNA cluster labels onto ATAC --------------------
pred <- TransferData(anchorset = anchors, refdata = rna$RNA_snn_res.0.8,
                     weight.reduction = atac[["lsi"]], dims = lsi_dims)
atac <- AddMetaData(atac, pred[, c("predicted.id", "prediction.score.max")])

## ---- impute expression onto every ATAC nucleus (the pseudo-multiome) --------
rna_data <- GetAssayData(rna, assay = "RNA", layer = "data")
keep <- rownames(rna_data)[Matrix::rowSums(rna_data > 0) >= min_cells]
cat("Imputing", length(keep), "genes onto", ncol(atac), "ATAC cells\n")

imputed <- TransferData(anchorset = anchors, refdata = rna_data[keep, ],
                        weight.reduction = atac[["lsi"]], dims = lsi_dims)
atac[["imputedRNA"]] <- imputed

## ---- plots ------------------------------------------------------------------
plot_genes <- intersect(marker_genes, keep)
missing <- setdiff(marker_genes, plot_genes)
if (length(missing)) message("Not found / too sparse in RNA: ", paste(missing, collapse = ", "))
h <- 4 * ceiling(length(plot_genes) / 4)

# labels: RNA clusters vs. the same labels transferred onto the ATAC UMAP
p_labels <- (DimPlot(rna, group.by = "RNA_snn_res.0.8", label = TRUE) +
               ggtitle("scRNA-seq clusters (RNA UMAP)") + NoLegend()) |
            (DimPlot(atac, group.by = "predicted.id", label = TRUE) +
               ggtitle("Transferred labels (ATAC UMAP)") + NoLegend()) |
            (FeaturePlot(atac, "prediction.score.max") + ggtitle("Transfer confidence"))
ggsave("results/figures/3dpf_label_transfer.pdf", p_labels, width = 18, height = 5.5)

# the key comparison: same genes on the RNA-only UMAP vs. the ATAC-defined UMAP
DefaultAssay(rna) <- "RNA"
p_rna <- FeaturePlot(rna, plot_genes, reduction = "umap", order = TRUE, ncol = 4)
ggsave("results/figures/3dpf_markers_RNA_umap.pdf", p_rna, width = 16, height = h)

DefaultAssay(atac) <- "imputedRNA"
p_atac <- FeaturePlot(atac, plot_genes, reduction = "umap", order = TRUE, ncol = 4)
ggsave("results/figures/3dpf_markers_imputed_ATAC_umap.pdf", p_atac, width = 16, height = h)

## ---- save -------------------------------------------------------------------
DefaultAssay(atac) <- "ATAC"
saveRDS(atac, "data/processed/3dpf_pseudomultiome.rds")
writeLines(capture.output(sessionInfo()), "results/sessionInfo_01_pseudomultiome_3dpf.txt")
