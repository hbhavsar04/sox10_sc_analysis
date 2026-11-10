## 02_opercle_osteo_subset.R
## Identify putative opercle osteoblast-lineage cells (osteoblasts + precursors in the
## dorsal hyoid / arch 2 condensation) in the Fabian et al. 2022 scRNA-seq and
## pseudo-multiome (snATAC + imputed expression) objects, count them per stage and
## replicate, and score maturation.
##
## Strategy: module scores for (1) osteoblast lineage, (2) arch 2 identity
## (hox PG2 high, PG3 low), (3) dorsal identity (dorsal minus ventral program).
## A cell is a candidate if it passes all three gates. Gates are crude per-cell
## filters, so the script also reports cluster-level mean scores: choosing whole
## clusters (especially in ATAC space, where position separates best) is usually
## more robust than per-cell gating for sparse markers.
##
## Inspect the histograms in results/figures/opercle/ and adjust `thresholds`
## and the gene panels before trusting the subsets.

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
})

setwd("~/Desktop/Hayden/git_repos/sox10_sc_analysis")
fig_dir <- "results/figures/opercle"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
set.seed(1)

## ---- parameters -------------------------------------------------------------
stages <- c("3dpf", "5dpf")   # stages without files in data/raw are skipped

# expected files (gunzip first): data/raw/<stage>_scRNA.rds,
#                                data/raw/<stage>_snATAC_SnapATAC.rds
rna_file  <- function(s) file.path("data/raw", paste0(s, "_scRNA.rds"))
atac_file <- function(s) file.path("data/raw", paste0(s, "_snATAC_SnapATAC.rds"))

gene_sets <- list(
  osteo_early = c("runx2a", "runx2b"),
  osteo_late  = c("sp7", "col10a1a", "ifitm5", "spp1", "bglap"),
  hox2        = c("hoxa2b", "hoxb2a"),             # arch 2 and posterior
  hox3        = c("hoxa3a", "hoxb3a", "hoxd3a"),   # arch 3 and posterior
  dorsal      = c("pou3f3b", "jag1b", "hey1"),
  ventral     = c("hand2", "dlx3b", "dlx5a")
)
gene_sets$osteo <- unique(c(gene_sets$osteo_early, gene_sets$osteo_late))

# gates on module scores (scores are relative to random control genes, so
# > 0 means above background). Tune after looking at the histograms.
thresholds <- list(
  osteo  = 0,   # osteo score
  arch2  = 0,   # hox2 - hox3
  dorsal = 0    # dorsal - ventral
)

## ---- helpers ----------------------------------------------------------------
score_sets <- function(obj, assay, sets) {
  for (nm in names(sets)) {
    g <- intersect(sets[[nm]], rownames(obj[[assay]]))
    miss <- setdiff(sets[[nm]], g)
    if (length(miss)) message("  [", nm, "] not found in ", assay, ": ", paste(miss, collapse = ", "))
    if (length(g) == 0) { obj@meta.data[[nm]] <- NA_real_; next }
    obj <- AddModuleScore(obj, features = list(g), name = paste0(nm, "_tmp"),
                          assay = assay, seed = 1)
    obj@meta.data[[nm]] <- obj@meta.data[[paste0(nm, "_tmp1")]]
    obj@meta.data[[paste0(nm, "_tmp1")]] <- NULL
  }
  obj$maturation   <- obj$osteo_late - obj$osteo_early
  obj$arch2_score  <- obj$hox2 - obj$hox3
  obj$dorsal_score <- obj$dorsal - obj$ventral
  obj$is_osteo  <- obj$osteo        > thresholds$osteo
  obj$is_arch2  <- obj$arch2_score  > thresholds$arch2
  obj$is_dorsal <- obj$dorsal_score > thresholds$dorsal
  obj$opercle_candidate <- obj$is_osteo & obj$is_arch2 & obj$is_dorsal
  obj
}

score_cols <- c("osteo", "maturation", "arch2_score", "dorsal_score")

plot_scores <- function(obj, label, cluster_col) {
  # score feature plots + candidate highlight
  p_feat <- FeaturePlot(obj, score_cols, reduction = "umap", order = TRUE, ncol = 2)
  cand <- colnames(obj)[obj$opercle_candidate %in% TRUE]
  p_cand <- DimPlot(obj, reduction = "umap", cells.highlight = cand,
                    cols.highlight = "firebrick", sizes.highlight = 0.6) +
    ggtitle(paste0(label, ": candidates (n = ", length(cand), ")")) + NoLegend()
  p_clu <- DimPlot(obj, reduction = "umap", group.by = cluster_col, label = TRUE) + NoLegend()
  ggsave(file.path(fig_dir, paste0(label, "_scores_umap.pdf")),
         p_feat / (p_clu | p_cand), width = 11, height = 15)

  # histograms with gate lines
  md <- obj@meta.data[, c("osteo", "arch2_score", "dorsal_score")] |>
    pivot_longer(everything(), names_to = "score", values_to = "value")
  gates <- data.frame(score = c("osteo", "arch2_score", "dorsal_score"),
                      value = unlist(thresholds[c("osteo", "arch2", "dorsal")]))
  p_hist <- ggplot(md, aes(value)) + geom_histogram(bins = 80) +
    geom_vline(data = gates, aes(xintercept = value), colour = "firebrick") +
    facet_wrap(~score, scales = "free") + theme_bw() + ggtitle(paste0(label, ": gates"))
  ggsave(file.path(fig_dir, paste0(label, "_score_histograms.pdf")), p_hist, width = 12, height = 4)
}

cluster_summary <- function(obj, cluster_col) {
  obj@meta.data |>
    group_by(cluster = .data[[cluster_col]]) |>
    summarise(n = n(),
              across(all_of(c("osteo", "maturation", "arch2_score", "dorsal_score")),
                     ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}"),
              frac_candidate = mean(opercle_candidate, na.rm = TRUE),
              .groups = "drop") |>
    arrange(desc(frac_candidate))
}

## ---- run per stage ----------------------------------------------------------
counts <- list()

for (s in stages) {
  if (!file.exists(rna_file(s)) || !file.exists(atac_file(s))) {
    message("Skipping ", s, ": files not found in data/raw (gunzip them first?)")
    next
  }
  message("== ", s, " ==")

  rna  <- UpdateSeuratObject(readRDS(rna_file(s)))
  atac <- UpdateSeuratObject(readRDS(atac_file(s)))
  DefaultAssay(rna)  <- "RNA"
  DefaultAssay(atac) <- "ACTIVITY"   # authors' imputed expression

  rna_clu  <- grep("^RNA_snn_res",  colnames(rna@meta.data),  value = TRUE)[1]
  atac_clu <- grep("^ATAC_snn_res", colnames(atac@meta.data), value = TRUE)[1]

  rna  <- score_sets(rna,  "RNA",      gene_sets)
  atac <- score_sets(atac, "ACTIVITY", gene_sets)

  for (x in list(list(obj = rna,  mod = "RNA",  clu = rna_clu),
                 list(obj = atac, mod = "ATAC", clu = atac_clu))) {
    label <- paste0(s, "_", x$mod)
    plot_scores(x$obj, label, x$clu)
    write.csv(cluster_summary(x$obj, x$clu),
              file.path("results", paste0(label, "_opercle_cluster_scores.csv")), row.names = FALSE)

    md <- x$obj@meta.data
    rep_col <- if ("bio.rep" %in% colnames(md)) md$bio.rep else "all"
    counts[[label]] <- data.frame(
      stage = s, modality = x$mod, bio.rep = rep_col,
      osteo = md$is_osteo, arch2 = md$is_arch2, dorsal = md$is_dorsal,
      candidate = md$opercle_candidate
    ) |>
      group_by(stage, modality, bio.rep) |>
      summarise(total = n(), osteo = sum(osteo, na.rm = TRUE),
                osteo_arch2 = sum(osteo & arch2, na.rm = TRUE),
                candidates = sum(candidate, na.rm = TRUE), .groups = "drop")

    cand <- colnames(x$obj)[x$obj$opercle_candidate %in% TRUE]
    if (length(cand) > 0) {
      saveRDS(subset(x$obj, cells = cand),
              file.path("data/processed", paste0(label, "_opercle_osteo_candidates.rds")))
    }
  }

  # within-candidate view: is there structure beyond maturation?
  cand <- colnames(atac)[atac$opercle_candidate %in% TRUE]
  if (length(cand) >= 10) {
    p_mat <- ggplot(atac@meta.data[cand, ], aes(maturation, dorsal_score, colour = arch2_score)) +
      geom_point() + scale_colour_viridis_c() + theme_bw() +
      ggtitle(paste0(s, " ATAC candidates: maturation vs dorsal identity"))
    ggsave(file.path(fig_dir, paste0(s, "_ATAC_candidates_maturation.pdf")), p_mat, width = 6, height = 5)
  }
}

counts_df <- bind_rows(counts)
print(counts_df)
write.csv(counts_df, "results/opercle_candidate_counts.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "results/sessionInfo_02_opercle_osteo_subset.txt")
