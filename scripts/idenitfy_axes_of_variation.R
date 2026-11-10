## 05_ob_axes_5dpf.R
## Step 2 (5 dpf): axes of variation within the osteoblast-lineage (OB) population.
##
##  1. Recompute TF-IDF/LSI on the OB nuclei only, so components capture variation
##     *within* the OB population; drop components that track sequencing depth.
##  2. Attach genes to each axis: correlation of imputed expression (ACTIVITY) with each
##     component. Structure comes from measured chromatin; genes only name it.
##  3. Attach TF motifs: correlation of chromVAR deviations with each component.
##  4. Annotate axes against landmarks (runx2a, sp7, nfatc1, abcc6a), maturation,
##     positional programs, cell cycle, stress, and depth.
##  5. Validate in raw RNA (OB cluster): do each axis's top genes co-vary in real counts
##     more than expression-matched random genes, and do opposite ends anticorrelate?
##
## Inputs (from 04_define_ob_clusters.R): data/processed/5dpf_ATAC_OB.rds, 5dpf_RNA_OB.rds
## Note: imputed expression is smoothed, so gene-axis correlations are inflated and the
## p-values are for ranking only. The RNA check (step 5) is the real test.

suppressPackageStartupMessages({
  library(Seurat); library(Signac); library(Matrix)
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(JASPAR2020); library(TFBSTools)
})

setwd("~/Desktop/Hayden/git_repos/sox10_sc_analysis")
stage   <- "5dpf"
fig_dir <- file.path("results/figures/ob_axes", stage)
tab_dir <- file.path("results/ob_axes", stage)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(1)

## ---- parameters -------------------------------------------------------------
atac_ob_file <- file.path("data/processed", paste0(stage, "_ATAC_OB.rds"))
rna_ob_file  <- file.path("data/processed", paste0(stage, "_RNA_OB.rds"))

min_peak_frac <- 0.02   # keep peaks accessible in >= 2% of OB nuclei
n_svd         <- 30
n_axes        <- 8      # number of non-depth components to interpret
depth_cor_max <- 0.5    # drop components with |r| > this vs log depth

min_mean_imputed <- 0.05   # genes considered for gene-axis correlation
top_n_genes      <- 30     # per axis end, written to tables
module_n         <- 25     # genes per axis end used for the RNA check
rna_min_pct      <- 0.05   # module genes must be detected in >= 5% of RNA OB cells
n_perm           <- 200    # random gene sets for the coherence null

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
cor_p <- function(r, n) { t <- r * sqrt((n - 2) / pmax(1 - r^2, 1e-12)); 2 * pt(-abs(t), n - 2) }

set_score <- function(mat, genes) {
  # mean of per-gene z-scores across cells; mat = genes x cells (dense)
  g <- intersect(genes, rownames(mat))
  if (length(g) == 0) return(rep(NA_real_, ncol(mat)))
  m <- mat[g, , drop = FALSE]
  m <- m[apply(m, 1, sd) > 0, , drop = FALSE]
  if (nrow(m) == 0) return(rep(NA_real_, ncol(mat)))
  colMeans(t(scale(t(m))))
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
  cc <- cor(t(m))
  mean(cc[upper.tri(cc)])
}

## ---- 1. LSI on OB nuclei only -----------------------------------------------
aob <- readRDS(atac_ob_file)
DefaultAssay(aob) <- "ATAC"
n_cells <- ncol(aob)
message(stage, " ATAC OB nuclei: ", n_cells)

cnt <- GetAssayData(aob, assay = "ATAC", layer = "counts")
pk  <- rownames(cnt)[Matrix::rowSums(cnt > 0) >= max(5, ceiling(min_peak_frac * n_cells))]
message("Peaks used for OB LSI: ", length(pk))

aob <- RunTFIDF(aob, assay = "ATAC", verbose = FALSE)
aob <- RunSVD(aob, assay = "ATAC", features = pk, n = n_svd,
              reduction.name = "lsi_ob", reduction.key = "LSIOB_", verbose = FALSE)

emb_all <- Embeddings(aob, "lsi_ob")
depth   <- if ("nCount_ATAC" %in% colnames(aob@meta.data)) aob$nCount_ATAC else Matrix::colSums(cnt)
dcor    <- as.vector(cor(emb_all, log10(depth)))
names(dcor) <- colnames(emb_all)
axes_idx <- which(abs(dcor) <= depth_cor_max)[seq_len(n_axes)]
axes_idx <- axes_idx[!is.na(axes_idx)]
axes     <- colnames(emb_all)[axes_idx]
emb      <- emb_all[, axes, drop = FALSE]
message("Depth-correlated components dropped: ",
        paste(names(dcor)[abs(dcor) > depth_cor_max], collapse = ", "))
message("Axes interpreted: ", paste(axes, collapse = ", "))

aob <- RunUMAP(aob, reduction = "lsi_ob", dims = axes_idx, reduction.name = "umap_ob",
               reduction.key = "UMAPOB_", verbose = FALSE)
aob <- FindNeighbors(aob, reduction = "lsi_ob", dims = axes_idx, verbose = FALSE)
aob <- FindClusters(aob, algorithm = 3, resolution = 0.5, verbose = FALSE)
aob$ob_subcluster <- as.character(Idents(aob))

p_depth <- data.frame(component = factor(names(dcor), levels = names(dcor)), r = dcor,
                      used = names(dcor) %in% axes) |>
  ggplot(aes(component, r, fill = used)) + geom_col() +
  geom_hline(yintercept = c(-1, 1) * depth_cor_max, linetype = 2) +
  scale_fill_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey70")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  ggtitle("OB-only LSI components: correlation with log depth")
p_umap <- DimPlot(aob, reduction = "umap_ob", group.by = "ob_subcluster", label = TRUE) +
  NoLegend() + ggtitle("OB nuclei, OB-only LSI")
p_axes <- FeaturePlot(aob, features = axes, reduction = "umap_ob", ncol = 4,
                      min.cutoff = "q2", max.cutoff = "q98")
ggsave(file.path(fig_dir, "01_lsi_ob_overview.pdf"),
       (p_depth | p_umap) / p_axes, width = 16, height = 14)

## ---- 2. genes per axis (imputed expression) ---------------------------------
imp <- GetAssayData(aob, assay = "ACTIVITY", layer = "data")
mu  <- Matrix::rowMeans(imp)
vv  <- Matrix::rowMeans(imp^2) - mu^2
genes_use <- rownames(imp)[mu >= min_mean_imputed & vv > 0]
imp_d <- as.matrix(imp[genes_use, ])
message("Genes tested against axes: ", length(genes_use))

gc <- cor(t(imp_d), emb, method = "spearman")       # genes x axes

rob <- readRDS(rna_ob_file)
DefaultAssay(rob) <- "RNA"
rcnt <- GetAssayData(rob, assay = "RNA", layer = "counts")
rna_pct <- Matrix::rowMeans(rcnt > 0)

gene_long <- as.data.frame(gc) |>
  mutate(gene = rownames(gc)) |>
  pivot_longer(-gene, names_to = "axis", values_to = "rho") |>
  mutate(p = cor_p(rho, n_cells)) |>
  group_by(axis) |> mutate(padj = p.adjust(p, "BH")) |> ungroup() |>
  mutate(rna_pct_detected = unname(rna_pct[gene]))
write.csv(gene_long, file.path(tab_dir, "gene_axis_correlations_all.csv"), row.names = FALSE)

top_genes <- gene_long |>
  group_by(axis) |>
  group_modify(~ bind_rows(
    slice_max(.x, rho, n = top_n_genes) |> mutate(end = "positive"),
    slice_min(.x, rho, n = top_n_genes) |> mutate(end = "negative"))) |>
  ungroup()
write.csv(top_genes, file.path(tab_dir, "gene_axis_top.csv"), row.names = FALSE)

p_genes <- top_genes |>
  group_by(axis, end) |> slice_max(abs(rho), n = 12) |> ungroup() |>
  mutate(axis = factor(axis, levels = axes),
         lab = paste(gene, axis, sep = "___"),
         lab = reorder(lab, rho)) |>
  ggplot(aes(rho, lab, fill = end)) + geom_col() +
  facet_wrap(~axis, scales = "free_y", ncol = 4) +
  scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
  scale_fill_manual(values = c(positive = "firebrick", negative = "steelblue")) +
  theme_bw() + labs(y = NULL, x = "Spearman rho (imputed expression vs axis)") +
  ggtitle(paste(stage, "top genes at each end of each axis"))
ggsave(file.path(fig_dir, "02_top_genes_per_axis.pdf"), p_genes, width = 16, height = 14)

## ---- 3. motifs per axis -----------------------------------------------------
cv <- as.matrix(GetAssayData(aob, assay = "chromvar", layer = "data"))
cv <- cv[apply(cv, 1, sd) > 0, , drop = FALSE]
mc <- cor(t(cv), emb, method = "spearman")
motif_long <- as.data.frame(mc) |>
  mutate(motif_id = rownames(mc), motif = motif_label(rownames(mc))) |>
  pivot_longer(-c(motif_id, motif), names_to = "axis", values_to = "rho") |>
  mutate(p = cor_p(rho, n_cells)) |>
  group_by(axis) |> mutate(padj = p.adjust(p, "BH")) |> ungroup()
write.csv(motif_long, file.path(tab_dir, "motif_axis_correlations_all.csv"), row.names = FALSE)

top_motifs <- motif_long |>
  group_by(axis) |>
  group_modify(~ bind_rows(slice_max(.x, rho, n = 5), slice_min(.x, rho, n = 5))) |>
  ungroup() |> pull(motif_id) |> unique()
p_motif <- motif_long |>
  filter(motif_id %in% top_motifs) |>
  mutate(axis = factor(axis, levels = axes)) |>
  ggplot(aes(axis, motif, fill = rho)) + geom_tile() +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle(paste(stage, "chromVAR motifs most correlated with each axis"))
ggsave(file.path(fig_dir, "03_motifs_per_axis.pdf"), p_motif,
       width = 9, height = max(6, 0.22 * length(top_motifs)))

## ---- 4. annotate axes ---------------------------------------------------------
imp_all <- as.matrix(imp[intersect(unique(c(landmarks, unlist(annot_sets))), rownames(imp)), , drop = FALSE])
annot <- data.frame(row.names = colnames(aob))
for (g in intersect(landmarks, rownames(imp_all))) annot[[g]] <- imp_all[g, ]
for (nm in names(annot_sets)) annot[[nm]] <- set_score(imp_all, annot_sets[[nm]])
annot$maturation   <- annot$osteo_late - annot$osteo_early
annot$hox2_vs_hox3 <- annot$hox2 - annot$hox3
annot$dorsal_vs_ventral <- annot$dorsal - annot$ventral
annot$log_depth    <- log10(depth)
annot <- annot[, colSums(!is.na(annot)) > 0, drop = FALSE]

ac <- cor(annot, emb, method = "spearman", use = "pairwise.complete.obs")
annot_long <- as.data.frame(ac) |> mutate(annotation = rownames(ac)) |>
  pivot_longer(-annotation, names_to = "axis", values_to = "rho")
write.csv(annot_long, file.path(tab_dir, "axis_annotation_atac.csv"), row.names = FALSE)

p_annot <- annot_long |>
  mutate(axis = factor(axis, levels = axes),
         annotation = factor(annotation, levels = rev(colnames(annot)))) |>
  ggplot(aes(axis, annotation, fill = rho)) + geom_tile() +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 2.8) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", limits = c(-1, 1)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle(paste(stage, "what each axis tracks (ATAC, imputed expression)"))
ggsave(file.path(fig_dir, "04_axis_annotation.pdf"), p_annot, width = 9, height = 7)

## ---- 5. RNA validation ------------------------------------------------------
rdat <- as.matrix(GetAssayData(rob, assay = "RNA", layer = "data"))
detected <- names(rna_pct)[rna_pct >= rna_min_pct]
rdat_d   <- rdat[detected, , drop = FALSE]
message("RNA OB cells: ", ncol(rob), "; genes detected in >= ",
        100 * rna_min_pct, "%: ", length(detected))

# expression bins for matched random gene sets
bins <- cut(rank(Matrix::rowMeans(rdat_d)), breaks = 20, labels = FALSE)
names(bins) <- detected
null_coh <- function(genes, n) {
  b <- table(bins[genes])
  replicate(n, {
    rg <- unlist(lapply(names(b), function(k) sample(names(bins)[bins == as.integer(k)], b[[k]])))
    coherence(rdat_d, rg)
  })
}

rna_land <- sapply(intersect(landmarks, rownames(rdat)), function(g) rdat[g, ])
rna_mat  <- set_score(rdat, annot_sets$osteo_late) - set_score(rdat, annot_sets$osteo_early)

val_rows <- list(); module_genes <- list()
for (ax in axes) {
  ranked <- gene_long |> filter(axis == ax, gene %in% detected)
  pos <- ranked |> filter(rho > 0) |> slice_max(rho, n = module_n) |> pull(gene)
  neg <- ranked |> filter(rho < 0) |> slice_min(rho, n = module_n) |> pull(gene)
  module_genes[[ax]] <- data.frame(axis = ax,
                                   end = c(rep("positive", length(pos)), rep("negative", length(neg))), gene = c(pos, neg))
  
  coh_pos <- coherence(rdat_d, pos); coh_neg <- coherence(rdat_d, neg)
  nul_pos <- if (length(pos) >= 3) null_coh(pos, n_perm) else NA
  nul_neg <- if (length(neg) >= 3) null_coh(neg, n_perm) else NA
  p_pos <- if (!is.na(coh_pos)) (sum(nul_pos >= coh_pos, na.rm = TRUE) + 1) / (n_perm + 1) else NA
  p_neg <- if (!is.na(coh_neg)) (sum(nul_neg >= coh_neg, na.rm = TRUE) + 1) / (n_perm + 1) else NA
  
  s_pos <- set_score(rdat_d, pos); s_neg <- set_score(rdat_d, neg)
  ax_rna <- s_pos - s_neg
  
  row <- data.frame(axis = ax, n_pos = length(pos), n_neg = length(neg),
                    coherence_pos = coh_pos, perm_p_pos = p_pos,
                    coherence_neg = coh_neg, perm_p_neg = p_neg,
                    cor_pos_vs_neg = suppressWarnings(cor(s_pos, s_neg, method = "spearman")),
                    rna_cor_maturation = suppressWarnings(cor(ax_rna, rna_mat, method = "spearman")))
  for (g in colnames(rna_land)) {
    row[[paste0("rna_cor_", g)]]  <- suppressWarnings(cor(ax_rna, rna_land[, g], method = "spearman"))
    row[[paste0("atac_cor_", g)]] <- if (g %in% rownames(ac)) ac[g, ax] else NA
  }
  val_rows[[ax]] <- row
}
validation <- bind_rows(val_rows)
write.csv(validation, file.path(tab_dir, "axis_rna_validation.csv"), row.names = FALSE)
write.csv(bind_rows(module_genes), file.path(tab_dir, "axis_module_genes.csv"), row.names = FALSE)
print(validation |> select(axis, coherence_pos, perm_p_pos, coherence_neg, perm_p_neg,
                           cor_pos_vs_neg, rna_cor_maturation), digits = 2)

# landmark agreement: ATAC axis-landmark correlation vs RNA axis-score-landmark correlation
lm_df <- bind_rows(lapply(colnames(rna_land), function(g) data.frame(
  axis = validation$axis, landmark = g,
  atac = validation[[paste0("atac_cor_", g)]], rna = validation[[paste0("rna_cor_", g)]])))
p_lm <- ggplot(lm_df, aes(atac, rna, colour = landmark, label = axis)) +
  geom_hline(yintercept = 0, colour = "grey60") + geom_vline(xintercept = 0, colour = "grey60") +
  geom_point(size = 2.5) + geom_text(size = 2.5, vjust = -0.8, show.legend = FALSE) +
  theme_bw() + labs(x = "ATAC: axis vs landmark (imputed)", y = "RNA: axis module score vs landmark") +
  ggtitle("Do axis-landmark relationships agree between modalities?")
ggsave(file.path(fig_dir, "05_rna_validation_landmarks.pdf"), p_lm, width = 8, height = 6)

## ---- save -------------------------------------------------------------------
saveRDS(aob, file.path("data/processed", paste0(stage, "_ATAC_OB_axes.rds")))
writeLines(capture.output(sessionInfo()), file.path(tab_dir, "sessionInfo.txt"))
message("Done. Tables in ", tab_dir, "; figures in ", fig_dir)
