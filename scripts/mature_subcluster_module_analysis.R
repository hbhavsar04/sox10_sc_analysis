## 08_mature_ob_srg_modules.R
## Candidate spatially restricted genes (SRGs) within MATURE osteoblasts (5 dpf).
##
##  1. Mature nuclei = chromatin subcluster 4 AND imputed subcluster 1 (script 06); "union" optional
##  2. Regress maturation, cell cycle and depth out of each gene's imputed expression
##  3. Most variable residual genes -> gene-gene correlation -> co-expression modules
##  4. Per module: how restricted (on in a subset of cells?), what it tracks (landmarks,
##     position, maturation, cell cycle, stress, depth), which motifs track it (chromVAR)
##  5. RNA check: module coherence in RNA mature-like cells and in all RNA OB cells,
##     against expression-matched random gene sets
##  6. Hub genes per module (correlation with module score), RNA detection, atlas SRG overlap
##
## Inputs: data/processed/5dpf_ATAC_OB_reclustered.rds (script 06), data/processed/5dpf_RNA_OB.rds
## Optional: data/reference/atlas_srgs.txt (one gene per line)

suppressPackageStartupMessages({
  library(Seurat); library(Signac); library(Matrix)
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(JASPAR2020); library(TFBSTools)
})

setwd("~/Desktop/Hayden/git_repos/sox10_sc_analysis")
stage   <- "5dpf"
fig_dir <- file.path("results/figures/mature_srg", stage)
tab_dir <- file.path("results/mature_srg", stage)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(1)

## ---- parameters -------------------------------------------------------------
mature_def  <- list(sub_chrom = "4", sub_imp = "1")
mature_mode <- "intersect"     # "intersect" (both partitions) or "union" (either)
min_cells   <- 25

min_mean_imputed  <- 0.05
n_var_genes       <- 400       # most variable residual genes used to build modules
module_cut_height <- 0.6       # cut dendrogram of (1 - Spearman r) here; lower = tighter modules
min_module_size   <- 8
regress_vars <- character(0)
exclude_file      <- ""   # suspected doublets, one per line

rna_mature_quantile <- 0.5     # RNA "mature-like" cells = top 50% of RNA OB cells by mature score
rna_min_pct         <- 0.05
n_perm              <- 500
techn_pattern       <- "^(rp[ls]|RP[LS]|mt-|MT-)"
atlas_srg_file      <- "data/reference/atlas_srgs.txt"

mature_genes <- c("col10a1a", "ifitm5", "panx3", "bgna", "col11a2", "slc8a4b")
landmarks    <- c("runx2a", "sp7", "nfatc1", "abcc6a")
annot_sets <- list(
  osteo_early = c("runx2a", "runx2b"),
  osteo_late  = c("sp7", "col10a1a", "ifitm5", "spp1", "bglap"),
  hox2        = c("hoxa2b", "hoxb2a"),
  hox3        = c("hoxa3a", "hoxb3a", "hoxd3a"),
  dorsal      = c("pou3f3b", "jag1b", "hey1"),
  ventral     = c("hand2", "dlx3b", "dlx5a"),
  cell_cycle  = c("mki67", "pcna", "top2a", "ccnb1", "cdk1", "mcm2"),
  stress      = c("fosab", "fosb", "junba", "jun", "egr1", "atf3", "ier2a",
                  "hsp70.3", "hsp70.2", "gadd45ba", "btg1"),
  translation = c("eef1g", "eef2b", "eef1a1l1", "faua", "tpt1", "uba52", "ppiaa")
)

## ---- helpers ----------------------------------------------------------------
zrows <- function(m) { m <- m[apply(m, 1, sd) > 0, , drop = FALSE]; t(scale(t(m))) }
set_score <- function(mat, genes) {
  g <- intersect(genes, rownames(mat))
  if (length(g) == 0) return(rep(NA_real_, ncol(mat)))
  z <- zrows(as.matrix(mat[g, , drop = FALSE]))
  if (nrow(z) == 0) return(rep(NA_real_, ncol(mat)))
  colMeans(z)
}
coherence <- function(mat, genes) {
  m <- as.matrix(mat[intersect(genes, rownames(mat)), , drop = FALSE])
  m <- m[apply(m, 1, sd) > 0, , drop = FALSE]
  if (nrow(m) < 3) return(NA_real_)
  cc <- cor(t(m)); mean(cc[upper.tri(cc)])
}
make_null <- function(mat) {
  bins <- cut(rank(Matrix::rowMeans(mat)), breaks = 20, labels = FALSE); names(bins) <- rownames(mat)
  function(genes, n) {
    genes <- intersect(genes, names(bins)); b <- table(bins[genes])
    replicate(n, {
      rg <- unlist(lapply(names(b), function(k) sample(names(bins)[bins == as.integer(k)], b[[k]])))
      coherence(mat, rg)
    })
  }
}
jaspar_names <- local({
  pfm <- getMatrixSet(JASPAR2020, opts = list(collection = "CORE", all_versions = TRUE))
  setNames(sapply(pfm, name), names(pfm))
})
motif_label <- function(ids) {
  ma <- sub("^.*?(MA[0-9]+[.-][0-9]+).*$", "\\1", ids, perl = TRUE)
  ma <- sub("-", ".", ma, fixed = TRUE)
  lab <- unname(jaspar_names[ma]); ifelse(is.na(lab), ids, lab)
}

## ---- 1. mature nuclei -------------------------------------------------------
aob <- readRDS(file.path("data/processed", paste0(stage, "_ATAC_OB_reclustered.rds")))
in_chrom <- aob$sub_chrom == mature_def$sub_chrom
in_imp   <- aob$sub_imp   == mature_def$sub_imp
message("Mature: chromatin ", sum(in_chrom), ", imputed ", sum(in_imp),
        ", both ", sum(in_chrom & in_imp), ", either ", sum(in_chrom | in_imp))
keep <- if (mature_mode == "intersect") in_chrom & in_imp else in_chrom | in_imp
if (file.exists(exclude_file)) {
  excl <- trimws(readLines(exclude_file)); excl <- excl[nzchar(excl)]
  message("Excluding ", sum(colnames(aob)[keep] %in% excl), " listed nuclei (", exclude_file, ")")
  keep <- keep & !(colnames(aob) %in% excl)
}
if (sum(keep) < min_cells) stop("Only ", sum(keep), " mature nuclei; try mature_mode <- 'union'")
mob <- subset(aob, cells = colnames(aob)[keep])
n <- ncol(mob)
message("Using ", n, " mature nuclei (", mature_mode, ")")

## ---- 2. covariates and residualization --------------------------------------
imp <- GetAssayData(mob, assay = "ACTIVITY", layer = "data")
mu <- Matrix::rowMeans(imp); vv <- Matrix::rowMeans(imp^2) - mu^2
genes_use <- names(mu)[mu >= min_mean_imputed & vv > 0 & !grepl(techn_pattern, names(mu))]
X <- as.matrix(imp[genes_use, ])

ann_mat <- as.matrix(imp[intersect(unique(c(landmarks, unlist(annot_sets))), rownames(imp)), , drop = FALSE])
annot <- data.frame(row.names = colnames(mob))
for (g in intersect(landmarks, rownames(ann_mat))) annot[[g]] <- ann_mat[g, ]
for (nm in names(annot_sets)) annot[[nm]] <- set_score(ann_mat, annot_sets[[nm]])
annot$maturation        <- annot$osteo_late - annot$osteo_early
annot$hox2_vs_hox3      <- annot$hox2 - annot$hox3
annot$dorsal_vs_ventral <- annot$dorsal - annot$ventral
annot$cc_score          <- annot$cell_cycle
depth <- if ("nCount_ATAC" %in% colnames(mob@meta.data)) mob$nCount_ATAC else
  Matrix::colSums(GetAssayData(mob, assay = "ATAC", layer = "counts"))
annot$log_depth <- log10(depth)

if (length(regress_vars) > 0) {
  cov <- as.matrix(annot[, regress_vars, drop = FALSE])
  cov[is.na(cov)] <- 0
  C <- cbind(1, scale(cov))
  C <- C[, colSums(is.na(C)) == 0, drop = FALSE]
  message("Regressing out: ", paste(regress_vars, collapse = ", "))
} else {
  C <- matrix(1, nrow = ncol(X), ncol = 1)               # intercept only = center each gene
  message("No covariates regressed (genes centered only)")
}
R <- X - (X %*% C) %*% solve(crossprod(C), t(C))      # residuals, genes x cells

## maturation-associated genes (center/edge candidates), reported separately
mat_cor <- cor(t(X), annot$maturation, method = "spearman")[, 1]
mat_tab <- data.frame(gene = names(mat_cor), rho_maturation = mat_cor, mean_imputed = mu[names(mat_cor)]) |>
  arrange(desc(abs(rho_maturation)))

## ---- 3. modules on residual variation ---------------------------------------
rvar <- apply(R, 1, var)
top  <- names(sort(rvar, decreasing = TRUE))[seq_len(min(n_var_genes, length(rvar)))]
cg   <- cor(t(R[top, ]), method = "spearman")
hc   <- hclust(as.dist(1 - cg), method = "average")
cut  <- cutree(hc, h = module_cut_height)
sizes <- sort(table(cut), decreasing = TRUE)
big   <- names(sizes)[sizes >= min_module_size]
if (length(big) == 0) stop("No modules >= ", min_module_size, " genes; raise module_cut_height")
mod_map <- setNames(paste0("M", seq_along(big)), big)
modules <- lapply(big, function(k) names(cut)[cut == k]); names(modules) <- unname(mod_map)
message("Modules: ", paste(sprintf("%s (%d genes)", names(modules), lengths(modules)), collapse = ", "))

# module scores (from un-residualized expression, so links to maturation etc. stay visible)
Xz <- zrows(X)
mscore <- sapply(modules, function(g) colMeans(Xz[intersect(g, rownames(Xz)), , drop = FALSE]))
mscore_res <- sapply(modules, function(g) colMeans(zrows(R[g, , drop = FALSE])))

restrict <- data.frame(
  module   = names(modules),
  n_genes  = lengths(modules),
  frac_cells_high = apply(scale(mscore_res), 2, function(s) mean(s > 1)),
  skewness = apply(mscore_res, 2, function(s) { z <- (s - mean(s)) / sd(s); mean(z^3) })
)

## hub genes (kME): correlation of each gene with its module's residual score
gene_tab <- bind_rows(lapply(names(modules), function(m) {
  g <- modules[[m]]
  data.frame(gene = g, module = m,
             kME = as.vector(cor(t(R[g, , drop = FALSE]), mscore_res[, m], method = "spearman")),
             residual_var = rvar[g], mean_imputed = mu[g],
             rho_maturation = mat_cor[g])
}))

## ---- 4. annotation + motifs -------------------------------------------------
ac <- cor(annot, mscore, method = "spearman", use = "pairwise.complete.obs")
annot_long <- as.data.frame(ac) |> mutate(annotation = rownames(ac)) |>
  pivot_longer(-annotation, names_to = "module", values_to = "rho")

cv <- as.matrix(GetAssayData(mob, assay = "chromvar", layer = "data"))
cv <- cv[apply(cv, 1, sd) > 0, , drop = FALSE]
mc <- cor(t(cv), mscore_res, method = "spearman")
motif_long <- as.data.frame(mc) |> mutate(motif_id = rownames(mc), motif = motif_label(rownames(mc))) |>
  pivot_longer(-c(motif_id, motif), names_to = "module", values_to = "rho")
top_motifs <- motif_long |> group_by(module) |> slice_max(rho, n = 5) |> ungroup()

## ---- 5. RNA check -----------------------------------------------------------
rob <- readRDS(file.path("data/processed", paste0(stage, "_RNA_OB.rds")))
DefaultAssay(rob) <- "RNA"
rdat <- GetAssayData(rob, assay = "RNA", layer = "data")
rna_pct <- Matrix::rowMeans(GetAssayData(rob, assay = "RNA", layer = "counts") > 0)
detected <- names(rna_pct)[rna_pct >= rna_min_pct & !grepl(techn_pattern, names(rna_pct))]

mat_rna <- set_score(rdat, mature_genes)
rna_mature_cells <- colnames(rob)[mat_rna >= quantile(mat_rna, rna_mature_quantile, na.rm = TRUE)]
r_all <- as.matrix(rdat[detected, ])
r_mat <- r_all[, rna_mature_cells, drop = FALSE]
r_mat <- r_mat[apply(r_mat, 1, sd) > 0, , drop = FALSE]
null_all <- make_null(r_all); null_mat <- make_null(r_mat)
message("RNA: ", ncol(r_all), " OB cells, ", ncol(r_mat), " mature-like cells")

val <- bind_rows(lapply(names(modules), function(m) {
  g_all <- intersect(modules[[m]], rownames(r_all))
  g_mat <- intersect(modules[[m]], rownames(r_mat))
  c_all <- coherence(r_all, g_all); c_mat <- coherence(r_mat, g_mat)
  p_all <- if (!is.na(c_all)) (sum(null_all(g_all, n_perm) >= c_all, na.rm = TRUE) + 1) / (n_perm + 1) else NA
  p_mat <- if (!is.na(c_mat)) (sum(null_mat(g_mat, n_perm) >= c_mat, na.rm = TRUE) + 1) / (n_perm + 1) else NA
  data.frame(module = m, n_genes = length(modules[[m]]), n_rna_detected = length(g_all),
             coherence_rna_mature = c_mat, perm_p_rna_mature = p_mat,
             coherence_rna_allOB = c_all, perm_p_rna_allOB = p_all)
}))

## ---- 6. gene table, atlas overlap, outputs -----------------------------------
atlas <- if (file.exists(atlas_srg_file)) { a <- trimws(readLines(atlas_srg_file)); a[nzchar(a)] } else NULL
gene_tab <- gene_tab |>
  mutate(rna_pct_detected = unname(rna_pct[gene]),
         in_atlas_srgs = if (is.null(atlas)) NA else gene %in% atlas) |>
  arrange(module, desc(kME))

summary_tab <- restrict |> left_join(val |> select(-n_genes), by = "module") |>
  left_join(annot_long |> filter(annotation %in% c("maturation", "cell_cycle", "stress", "log_depth",
                                                   "hox2_vs_hox3", "dorsal_vs_ventral", landmarks)) |>
              pivot_wider(names_from = annotation, values_from = rho, names_prefix = "rho_"),
            by = "module") |>
  left_join(gene_tab |> group_by(module) |> slice_max(kME, n = 8) |>
              summarise(hub_genes = paste(gene, collapse = ", ")), by = "module")

write.csv(summary_tab, file.path(tab_dir, "module_summary.csv"), row.names = FALSE)
write.csv(gene_tab,    file.path(tab_dir, "module_genes.csv"),   row.names = FALSE)
write.csv(mat_tab,     file.path(tab_dir, "maturation_associated_genes.csv"), row.names = FALSE)
write.csv(motif_long,  file.path(tab_dir, "module_motif_correlations.csv"), row.names = FALSE)
print(tibble::as_tibble(summary_tab) |>
        select(module, n_genes, frac_cells_high, skewness,
               perm_p_rna_mature, perm_p_rna_allOB, hub_genes), width = Inf)

## ---- plots ------------------------------------------------------------------
ord  <- unlist(modules)
cg_o <- cg[ord, ord]
cg_df <- as.data.frame(as.table(cg_o)) |> setNames(c("g1", "g2", "r"))
p_cor <- ggplot(cg_df, aes(g1, g2, fill = r)) + geom_raster() +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", limits = c(-1, 1)) +
  theme_void() + ggtitle("Residual gene-gene correlation, ordered by module")

p_ann <- annot_long |>
  mutate(annotation = factor(annotation, levels = rev(colnames(annot)))) |>
  ggplot(aes(module, annotation, fill = rho)) + geom_tile() +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 2.6) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", limits = c(-1, 1)) +
  theme_bw() + labs(y = NULL) + ggtitle("What each module tracks (mature nuclei)")

p_mot <- motif_long |> filter(motif_id %in% top_motifs$motif_id) |>
  ggplot(aes(module, motif, fill = rho)) + geom_tile() +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick") +
  theme_bw() + labs(y = NULL) + ggtitle("Motifs tracking each module")

ggsave(file.path(fig_dir, "01_modules.pdf"), (p_cor | p_ann) / p_mot, width = 15, height = 14)

for (m in names(modules)) mob@meta.data[[paste0("score_", m)]] <- mscore_res[, m]
p_umap <- FeaturePlot(mob, paste0("score_", names(modules)), reduction = "umap_ob",
                      ncol = 4, min.cutoff = "q2", max.cutoff = "q98")
ggsave(file.path(fig_dir, "02_module_scores_umap.pdf"), p_umap,
       width = 16, height = 4 * ceiling(length(modules) / 4))

saveRDS(mob, file.path("data/processed", paste0(stage, "_ATAC_mature_OB_modules.rds")))
writeLines(capture.output(sessionInfo()), file.path(tab_dir, "sessionInfo.txt"))
message("Done. Tables in ", tab_dir, "; figures in ", fig_dir)
