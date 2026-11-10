## 07_core_plus_partner_genes.R
## Find genes expressed in the core (mature) osteoblast subcluster AND one other specific
## subcluster, but low everywhere else. Such "core + partner" genes may mark a region
## shared by different cell types (positional genes), an interface, or a shared state.
##
## For each partition (chromatin / imputed subclusters from script 06):
##  1. mean imputed expression per subcluster
##  2. correlate each gene's subcluster profile with "core + partner k" templates; keep best k
##  3. require separation (core and partner both above every other subcluster) and balance
##     (core and partner expression comparable)
##  4. cell-level Wilcoxon tests: core vs rest and partner vs rest (rest excludes both)
##  5. RNA check: detection in RNA OB cells and correlation with the mature-OB program
##  6. genes found in both partitions = most robust
##
## Inputs: data/processed/5dpf_ATAC_OB_reclustered.rds (script 06), data/processed/5dpf_RNA_OB.rds
## Note: imputed values are smoothed, so Wilcoxon p-values are inflated; use them as filters,
## not as evidence. Imputation can also bleed expression between neighbouring subclusters.

suppressPackageStartupMessages({
  library(Seurat); library(Signac); library(Matrix)
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})

setwd("~/Desktop/Hayden/git_repos/sox10_sc_analysis")
stage   <- "5dpf"
fig_dir <- file.path("results/figures/core_partner", stage)
tab_dir <- file.path("results/core_partner", stage)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(1)

## ---- parameters -------------------------------------------------------------
# partition column -> core (mature OB) subcluster; check these against your heatmaps
partitions <- list(sub_chrom = "4", sub_imp = "1")

min_cluster_mean <- 0.1    # gene must reach this mean imputed value in at least one subcluster
min_balance      <- 0.4    # weaker of core/partner must be >= 40% of the stronger
top_n_test       <- 300    # best template matches carried into cell-level tests
padj_max         <- 0.05
n_heatmap        <- 40     # genes shown in the heatmap
techn_pattern    <- "^(rp[ls]|RP[LS]|mt-|MT-)"
mature_genes     <- c("col10a1a", "ifitm5", "panx3", "bgna", "col11a2", "slc8a4b")

## ---- load -------------------------------------------------------------------
aob <- readRDS(file.path("data/processed", paste0(stage, "_ATAC_OB_reclustered.rds")))
rob <- readRDS(file.path("data/processed", paste0(stage, "_RNA_OB.rds")))
DefaultAssay(rob) <- "RNA"

imp     <- GetAssayData(aob, assay = "ACTIVITY", layer = "data")
rdat    <- GetAssayData(rob, assay = "RNA", layer = "data")
rna_pct <- Matrix::rowMeans(GetAssayData(rob, assay = "RNA", layer = "counts") > 0)

set_score <- function(mat, genes) {
  g <- intersect(genes, rownames(mat))
  m <- as.matrix(mat[g, , drop = FALSE])
  m <- m[apply(m, 1, sd) > 0, , drop = FALSE]
  colMeans(t(scale(t(m))))
}
mature_rna <- set_score(rdat, mature_genes)

## ---- core function ----------------------------------------------------------
find_core_partner <- function(col, core) {
  grp <- as.character(aob@meta.data[[col]])
  if (!core %in% grp) stop("Core subcluster ", core, " not found in ", col)
  cl <- as.character(sort(as.numeric(unique(grp))))
  
  M <- sapply(cl, function(k) Matrix::rowMeans(imp[, grp == k, drop = FALSE]))   # genes x subclusters
  M <- M[apply(M, 1, max) >= min_cluster_mean & !grepl(techn_pattern, rownames(M)), , drop = FALSE]
  Z <- t(scale(t(M)))
  ok <- rowSums(is.na(Z)) == 0
  M <- M[ok, , drop = FALSE]; Z <- Z[ok, , drop = FALSE]
  
  partners <- setdiff(cl, core)
  Tm <- sapply(partners, function(k) as.numeric(cl %in% c(core, k)))
  rownames(Tm) <- cl
  R <- cor(t(Z), Tm)                                              # genes x partners
  
  best_i  <- max.col(R, ties.method = "first")
  partner <- partners[best_i]
  i_core  <- match(core, cl)
  i_part  <- match(partner, cl)
  mean_core    <- M[, i_core]
  mean_partner <- M[cbind(seq_len(nrow(M)), i_part)]
  max_other    <- vapply(seq_len(nrow(M)), function(i) max(M[i, -c(i_core, i_part[i])]), numeric(1))
  
  res <- data.frame(gene = rownames(M), partner = partner,
                    r_template = R[cbind(seq_len(nrow(R)), best_i)],
                    mean_core = mean_core, mean_partner = mean_partner, max_other = max_other) |>
    mutate(separation    = pmin(mean_core, mean_partner) - max_other,
           balance       = pmin(mean_core, mean_partner) / pmax(mean_core, mean_partner),
           fold_vs_other = pmin(mean_core, mean_partner) / (max_other + 1e-3))
  
  cand <- res |>
    filter(separation > 0, balance >= min_balance) |>
    arrange(desc(r_template)) |>
    head(top_n_test)
  
  w_greater <- function(x, a, b) suppressWarnings(wilcox.test(x[a], x[b], alternative = "greater")$p.value)
  cand$p_core <- mapply(function(g, k) {
    x <- imp[g, ]; w_greater(x, grp == core, !(grp %in% c(core, k)))
  }, cand$gene, cand$partner)
  cand$p_partner <- mapply(function(g, k) {
    x <- imp[g, ]; w_greater(x, grp == k, !(grp %in% c(core, k)))
  }, cand$gene, cand$partner)
  
  cand <- cand |>
    mutate(padj_core    = p.adjust(p_core, "BH"),
           padj_partner = p.adjust(p_partner, "BH"),
           passes       = padj_core < padj_max & padj_partner < padj_max,
           rna_pct_detected = unname(rna_pct[gene]),
           rna_cor_mature   = vapply(gene, function(g) {
             if (!g %in% rownames(rdat)) return(NA_real_)
             x <- as.numeric(rdat[g, ]); if (sd(x) == 0) return(NA_real_)
             suppressWarnings(cor(x, mature_rna, method = "spearman"))
           }, numeric(1)),
           partition = col) |>
    arrange(partner, desc(r_template))
  
  list(all = res, cand = cand, M = M, cl = cl, grp = grp, core = core, col = col)
}

## ---- plots ------------------------------------------------------------------
plot_results <- function(out) {
  hits <- out$cand |> filter(passes)
  label <- out$col
  
  p_count <- hits |> count(partner) |>
    ggplot(aes(factor(partner, levels = out$cl), n)) + geom_col(fill = "firebrick") +
    theme_bw() + labs(x = "Partner subcluster", y = "Genes shared with core",
                      title = paste0(label, ": core = ", out$core))
  
  top <- hits |> group_by(partner) |>
    slice_max(r_template, n = ceiling(n_heatmap / max(1, n_distinct(hits$partner)))) |>
    ungroup() |> arrange(partner, desc(r_template))
  if (nrow(top) == 0) { message(label, ": no genes passed"); return(invisible(NULL)) }
  
  Mz <- t(scale(t(out$M[top$gene, , drop = FALSE])))
  hm <- as.data.frame(Mz) |> mutate(gene = rownames(Mz)) |>
    pivot_longer(-gene, names_to = "subcluster", values_to = "z") |>
    left_join(top |> select(gene, partner), by = "gene") |>
    mutate(row = factor(paste0(gene, "  [", partner, "]"),
                        levels = rev(paste0(top$gene, "  [", top$partner, "]"))),
           subcluster = factor(subcluster, levels = out$cl))
  p_hm <- ggplot(hm, aes(subcluster, row, fill = z)) + geom_tile(colour = "grey90") +
    geom_vline(xintercept = match(out$core, out$cl) + c(-0.5, 0.5), linewidth = 0.8) +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick") +
    theme_bw() + labs(y = "gene  [partner]", x = "Subcluster (core outlined)",
                      title = "Core + partner genes (mean imputed, z across subclusters)")
  
  ggsave(file.path(fig_dir, paste0(label, "_core_partner.pdf")),
         p_count / p_hm + plot_layout(heights = c(1, 4)),
         width = 9, height = max(8, 0.22 * nrow(top) + 4))
  
  # top gene per partner on the chromatin UMAP
  feat <- hits |> group_by(partner) |> slice_max(r_template, n = 1) |> pull(gene)
  if (length(feat)) {
    DefaultAssay(aob) <- "ACTIVITY"
    p_f <- FeaturePlot(aob, feat, reduction = "umap_ob", order = TRUE, ncol = 3) |
      (DimPlot(aob, reduction = "umap_ob", group.by = label, label = TRUE) + NoLegend())
    ggsave(file.path(fig_dir, paste0(label, "_top_gene_per_partner_umap.pdf")), p_f,
           width = 16, height = 4 * ceiling(length(feat) / 3) + 1)
  }
}

## ---- run --------------------------------------------------------------------
results <- lapply(names(partitions), function(col) find_core_partner(col, partitions[[col]]))
names(results) <- names(partitions)

for (out in results) {
  write.csv(out$cand, file.path(tab_dir, paste0(out$col, "_core_partner_candidates.csv")), row.names = FALSE)
  plot_results(out)
  message("\n", out$col, " (core = ", out$core, "): ", sum(out$cand$passes), " genes pass; by partner:")
  print(table(partner = out$cand$partner[out$cand$passes]))
}

## genes passing in both partitions (most robust)
both <- inner_join(
  results$sub_chrom$cand |> filter(passes) |> select(gene, partner_chrom = partner, r_chrom = r_template,
                                                     rna_pct_detected, rna_cor_mature),
  results$sub_imp$cand   |> filter(passes) |> select(gene, partner_imp = partner, r_imp = r_template),
  by = "gene") |>
  arrange(desc(r_chrom + r_imp))
write.csv(both, file.path(tab_dir, "core_partner_both_partitions.csv"), row.names = FALSE)
message("\nGenes passing in both partitions: ", nrow(both))
print(head(both, 30))

# which chromatin partner corresponds to which imputed partner (by shared cells)
message("\nCell overlap between chromatin and imputed subclusters:")
print(table(chromatin = aob$sub_chrom, imputed = aob$sub_imp))

writeLines(capture.output(sessionInfo()), file.path(tab_dir, "sessionInfo.txt"))