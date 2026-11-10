## 06_recluster_atac_ob.R
## Recluster the 5 dpf ATAC OB nuclei two ways and compare:
##   A. chromatin: OB-only LSI (from script 05), depth components excluded
##   B. imputed expression (ACTIVITY): variable genes (ribosomal/mito removed), PCA,
##      optional cell-cycle regression
## For each partition: resolution sweep, agreement between partitions (ARI),
## subcluster annotation (maturation, landmarks, position, cell cycle, stress, depth),
## marker genes (imputed) and motifs (chromVAR), and an RNA check of each subcluster's
## marker genes (do they co-vary in raw RNA counts more than matched random genes?).
##
## Inputs: data/processed/5dpf_ATAC_OB_axes.rds (script 05), data/processed/5dpf_RNA_OB.rds

suppressPackageStartupMessages({
  library(Seurat); library(Signac); library(Matrix)
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(JASPAR2020); library(TFBSTools)
})

setwd("~/Desktop/Hayden/git_repos/sox10_sc_analysis")
stage   <- "5dpf"
fig_dir <- file.path("results/figures/ob_recluster", stage)
tab_dir <- file.path("results/ob_recluster", stage)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(1)

## ---- parameters -------------------------------------------------------------
res_grid      <- c(0.3, 0.5, 0.8, 1.0)
chrom_res_use <- 0.5      # resolution used for downstream analysis of each partition
imp_res_use   <- 0.5
n_axes        <- 8        # same as script 05
depth_cor_max <- 0.5
imp_pcs       <- 1:15
n_hvg         <- 2000
min_mean_imputed <- 0.05
regress_cell_cycle <- TRUE  # regress cell-cycle score out of imputed expression before PCA

module_n    <- 20      # top markers per subcluster for the RNA check
rna_min_pct <- 0.05
n_perm      <- 500

techn_pattern <- "^(rp[ls]|RP[LS]|mt-|MT-)"   # excluded from HVGs and markers

landmarks <- c("runx2a", "sp7", "nfatc1", "abcc6a")
annot_sets <- list(
  osteo_early = c("runx2a", "runx2b"),
  osteo_late  = c("sp7", "col10a1a", "ifitm5", "spp1", "bglap"),
  hox2        = c("hoxa2b", "hoxb2a"),
  hox3        = c("hoxa3a", "hoxb3a", "hoxd3a"),
  dorsal      = c("pou3f3b", "jag1b", "hey1"),
  ventral     = c("hand2", "dlx3b", "dlx5a"),
  cell_cycle  = c("mki67", "pcna", "top2a", "ccnb1", "cdk1", "mcm2"),
  stress      = c("fosab", "junba", "jun", "egr1", "atf3")
)

## ---- helpers ----------------------------------------------------------------
set_score <- function(mat, genes) {
  g <- intersect(genes, rownames(mat))
  if (length(g) == 0) return(rep(NA_real_, ncol(mat)))
  m <- mat[g, , drop = FALSE]
  m <- m[apply(m, 1, sd) > 0, , drop = FALSE]
  if (nrow(m) == 0) return(rep(NA_real_, ncol(mat)))
  colMeans(t(scale(t(m))))
}

ari <- function(a, b) {
  tab <- table(a, b); n <- sum(tab); s2 <- function(x) sum(choose(x, 2))
  sij <- s2(tab); sa <- s2(rowSums(tab)); sb <- s2(colSums(tab))
  e <- sa * sb / choose(n, 2)
  (sij - e) / ((sa + sb) / 2 - e)
}

jaspar_names <- local({
  pfm <- getMatrixSet(JASPAR2020, opts = list(collection = "CORE", all_versions = TRUE))
  setNames(sapply(pfm, name), names(pfm))
})
motif_label <- function(ids) {
  ma <- sub("^.*?(MA[0-9]+[.-][0-9]+).*$", "\\1", ids, perl = TRUE)
  ma <- sub("-", ".", ma, fixed = TRUE)
  lab <- unname(jaspar_names[ma])
  ifelse(is.na(lab), ids, lab)
}

coherence <- function(mat, genes) {
  m <- mat[intersect(genes, rownames(mat)), , drop = FALSE]
  m <- m[apply(m, 1, sd) > 0, , drop = FALSE]
  if (nrow(m) < 3) return(NA_real_)
  cc <- cor(t(m)); mean(cc[upper.tri(cc)])
}

## ---- load -------------------------------------------------------------------
aob <- readRDS(file.path("data/processed", paste0(stage, "_ATAC_OB_axes.rds")))
rob <- readRDS(file.path("data/processed", paste0(stage, "_RNA_OB.rds")))
DefaultAssay(rob) <- "RNA"
message(stage, ": ", ncol(aob), " ATAC OB nuclei, ", ncol(rob), " RNA OB cells")

imp <- GetAssayData(aob, assay = "ACTIVITY", layer = "data")
cnt <- GetAssayData(aob, assay = "ATAC", layer = "counts")
depth <- if ("nCount_ATAC" %in% colnames(aob@meta.data)) aob$nCount_ATAC else Matrix::colSums(cnt)

# per-cell annotation scores (imputed expression)
ann_genes <- intersect(unique(c(landmarks, unlist(annot_sets))), rownames(imp))
imp_ann <- as.matrix(imp[ann_genes, , drop = FALSE])
annot <- data.frame(row.names = colnames(aob))
for (g in intersect(landmarks, ann_genes)) annot[[g]] <- imp_ann[g, ]
for (nm in names(annot_sets)) annot[[nm]] <- set_score(imp_ann, annot_sets[[nm]])
annot$maturation        <- annot$osteo_late - annot$osteo_early
annot$hox2_vs_hox3      <- annot$hox2 - annot$hox3
annot$dorsal_vs_ventral <- annot$dorsal - annot$ventral
annot$log_depth         <- log10(depth)
aob$cc_score <- annot$cell_cycle

## ---- A. chromatin partition -------------------------------------------------
emb_all <- Embeddings(aob, "lsi_ob")
dcor <- as.vector(cor(emb_all, log10(depth)))
axes_idx <- which(abs(dcor) <= depth_cor_max)[seq_len(n_axes)]
axes_idx <- axes_idx[!is.na(axes_idx)]

DefaultAssay(aob) <- "ATAC"
aob <- FindNeighbors(aob, reduction = "lsi_ob", dims = axes_idx,
                     graph.name = c("chrom_nn", "chrom_snn"), verbose = FALSE)
for (r in res_grid) aob <- FindClusters(aob, graph.name = "chrom_snn", algorithm = 3,
                                        resolution = r, verbose = FALSE)

## ---- B. imputed-expression partition ----------------------------------------
DefaultAssay(aob) <- "ACTIVITY"
mu <- Matrix::rowMeans(imp); vv <- Matrix::rowMeans(imp^2) - mu^2
cand <- names(mu)[mu >= min_mean_imputed & vv > 0 & !grepl(techn_pattern, names(mu))]
hvg  <- head(cand[order(vv[cand], decreasing = TRUE)], n_hvg)
VariableFeatures(aob[["ACTIVITY"]]) <- hvg

aob <- ScaleData(aob, assay = "ACTIVITY", features = hvg,
                 vars.to.regress = if (regress_cell_cycle) "cc_score" else NULL, verbose = FALSE)
aob <- RunPCA(aob, assay = "ACTIVITY", features = hvg, npcs = 30,
              reduction.name = "pca_imp", reduction.key = "PCIMP_", verbose = FALSE)
aob <- FindNeighbors(aob, reduction = "pca_imp", dims = imp_pcs,
                     graph.name = c("imp_nn", "imp_snn"), verbose = FALSE)
for (r in res_grid) aob <- FindClusters(aob, graph.name = "imp_snn", resolution = r, verbose = FALSE)
aob <- RunUMAP(aob, reduction = "pca_imp", dims = imp_pcs,
               reduction.name = "umap_imp", reduction.key = "UMAPIMP_", verbose = FALSE)

## ---- resolution sweep + agreement -------------------------------------------
sweep <- bind_rows(lapply(res_grid, function(r) {
  a <- aob@meta.data[[paste0("chrom_snn_res.", r)]]
  b <- aob@meta.data[[paste0("imp_snn_res.", r)]]
  data.frame(resolution = r,
             n_chrom = length(unique(a)), min_size_chrom = min(table(a)),
             n_imp   = length(unique(b)), min_size_imp   = min(table(b)),
             ARI_chrom_vs_imp = ari(a, b))
}))
print(sweep)
write.csv(sweep, file.path(tab_dir, "resolution_sweep.csv"), row.names = FALSE)

aob$sub_chrom <- as.character(aob@meta.data[[paste0("chrom_snn_res.", chrom_res_use)]])
aob$sub_imp   <- as.character(aob@meta.data[[paste0("imp_snn_res.",   imp_res_use)]])

ct <- as.data.frame(table(chromatin = aob$sub_chrom, imputed = aob$sub_imp))
p_ct <- ggplot(ct, aes(imputed, chromatin, fill = Freq)) + geom_tile() +
  geom_text(aes(label = ifelse(Freq > 0, Freq, "")), size = 3) +
  scale_fill_gradient(low = "white", high = "firebrick") + theme_bw() +
  ggtitle(sprintf("Chromatin vs imputed subclusters (ARI = %.2f)", ari(aob$sub_chrom, aob$sub_imp)))

p_u <- (DimPlot(aob, reduction = "umap_ob",  group.by = "sub_chrom", label = TRUE) + NoLegend() +
          ggtitle("Chromatin UMAP, chromatin subclusters")) |
  (DimPlot(aob, reduction = "umap_ob",  group.by = "sub_imp",   label = TRUE) + NoLegend() +
     ggtitle("Chromatin UMAP, imputed subclusters")) |
  (DimPlot(aob, reduction = "umap_imp", group.by = "sub_imp",   label = TRUE) + NoLegend() +
     ggtitle("Imputed UMAP, imputed subclusters"))
ggsave(file.path(fig_dir, "01_subclusters_overview.pdf"), p_u / p_ct, width = 16, height = 11)

## ---- per-partition characterization -----------------------------------------
rdat <- as.matrix(GetAssayData(rob, assay = "RNA", layer = "data"))
rna_pct <- Matrix::rowMeans(GetAssayData(rob, assay = "RNA", layer = "counts") > 0)
detected <- names(rna_pct)[rna_pct >= rna_min_pct & !grepl(techn_pattern, names(rna_pct))]
rdat_d <- rdat[detected, , drop = FALSE]
bins <- cut(rank(Matrix::rowMeans(rdat_d)), breaks = 20, labels = FALSE); names(bins) <- detected
null_coh <- function(genes, n) {
  b <- table(bins[genes])
  replicate(n, {
    rg <- unlist(lapply(names(b), function(k) sample(names(bins)[bins == as.integer(k)], b[[k]])))
    coherence(rdat_d, rg)
  })
}

marker_genes_pool <- names(mu)[mu >= min_mean_imputed & vv > 0 & !grepl(techn_pattern, names(mu))]

characterize <- function(col, label) {
  Idents(aob) <- col
  grp <- aob@meta.data[[col]]
  
  # annotation means per subcluster (z-scored across subclusters)
  am <- annot |> mutate(subcluster = grp) |>
    group_by(subcluster) |> summarise(across(everything(), ~ mean(.x, na.rm = TRUE)),
                                      n_cells = n(), .groups = "drop")
  write.csv(am, file.path(tab_dir, paste0(label, "_subcluster_annotation.csv")), row.names = FALSE)
  am_long <- am |> select(-n_cells) |>
    pivot_longer(-subcluster, names_to = "annotation", values_to = "mean") |>
    group_by(annotation) |> mutate(z = (mean - mean(mean)) / sd(mean)) |> ungroup()
  p_am <- ggplot(am_long, aes(subcluster, factor(annotation, levels = rev(colnames(annot))), fill = z)) +
    geom_tile() + scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick") +
    theme_bw() + labs(y = NULL) + ggtitle(paste(label, "subcluster annotation (z across subclusters)"))
  
  # marker genes (imputed; p-values inflated by smoothing, use for ranking)
  mk <- FindAllMarkers(aob, assay = "ACTIVITY", features = marker_genes_pool, only.pos = TRUE,
                       min.pct = 0, logfc.threshold = 0.25, verbose = FALSE) |>
    mutate(rna_pct_detected = unname(rna_pct[gene]))
  write.csv(mk, file.path(tab_dir, paste0(label, "_markers_imputed.csv")), row.names = FALSE)
  
  top5 <- mk |> group_by(cluster) |> slice_max(avg_log2FC, n = 5) |> ungroup()
  hm <- GetAssayData(aob, assay = "ACTIVITY", layer = "data")[unique(top5$gene), , drop = FALSE]
  hm <- sapply(split(seq_len(ncol(aob)), grp), function(i) Matrix::rowMeans(hm[, i, drop = FALSE]))
  if (is.null(dim(hm))) hm <- matrix(hm, nrow = 1, dimnames = list(unique(top5$gene), names(hm)))
  hm_df <- as.data.frame(t(scale(t(hm)))) |> mutate(gene = rownames(hm)) |>
    pivot_longer(-gene, names_to = "subcluster", values_to = "z")
  p_mk <- ggplot(hm_df, aes(subcluster, factor(gene, levels = rev(unique(top5$gene))), fill = z)) +
    geom_tile() + scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick") +
    theme_bw() + labs(y = NULL) + ggtitle(paste(label, "top imputed markers"))
  
  # motifs
  mo <- FindAllMarkers(aob, assay = "chromvar", only.pos = TRUE, mean.fxn = rowMeans,
                       fc.name = "avg_diff", logfc.threshold = 0, min.pct = 0, verbose = FALSE) |>
    mutate(motif = motif_label(gene))
  write.csv(mo, file.path(tab_dir, paste0(label, "_motifs.csv")), row.names = FALSE)
  
  ggsave(file.path(fig_dir, paste0("02_", label, "_characterization.pdf")),
         p_am | p_mk, width = 15, height = max(7, 0.2 * nrow(top5)))
  
  # RNA check: do each subcluster's top markers co-vary in raw RNA counts?
  val <- bind_rows(lapply(sort(unique(grp)), function(k) {
    g <- mk |> filter(cluster == k, gene %in% detected) |>
      slice_max(avg_log2FC, n = module_n) |> pull(gene)
    coh <- coherence(rdat_d, g)
    nul <- if (length(g) >= 3) null_coh(g, n_perm) else NA
    data.frame(subcluster = k, n_cells = sum(grp == k), n_marker_genes = length(g),
               coherence = coh,
               perm_p = if (is.na(coh)) NA else (sum(nul >= coh, na.rm = TRUE) + 1) / (n_perm + 1),
               genes = paste(head(g, 12), collapse = ", "))
  }))
  write.csv(val, file.path(tab_dir, paste0(label, "_rna_validation.csv")), row.names = FALSE)
  message("\n", label, " subclusters: RNA check of marker modules")
  print(val |> select(-genes), digits = 3)
  invisible(val)
}

val_chrom <- characterize("sub_chrom", "chromatin")
val_imp   <- characterize("sub_imp",   "imputed")

## ---- save -------------------------------------------------------------------
saveRDS(aob, file.path("data/processed", paste0(stage, "_ATAC_OB_reclustered.rds")))
writeLines(capture.output(sessionInfo()), file.path(tab_dir, "sessionInfo.txt"))
message("Done. Tables in ", tab_dir, "; figures in ", fig_dir)