## 02_opercle_osteo_subset.R
## Identify putative opercle osteoblast-lineage cells (osteoblasts + precursors in the
## dorsal hyoid / arch 2 condensation) in the Fabian et al. 2022 scRNA-seq and
## pseudo-multiome (snATAC + imputed expression) objects, count them per stage and
## replicate, and score maturation.
##
## Gates (all must pass):
##   osteo  : osteoblast-lineage module score
##   arch 2 : hox PG2 score expressed AND hox PG2 > hox PG3
##   dorsal : dorsal minus ventral program
##
## RNA gates use fixed thresholds on module scores. Imputed (ATAC) expression is
## smoothed, which inflates module scores, so ATAC gates are calibrated: each ATAC
## threshold is set so the same fraction of cells passes as in the RNA data for
## that stage.
##
## Per-cell gates are crude. The cluster-level score tables are usually the
## better way to choose the population, especially in ATAC space.

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

# RNA thresholds on module scores (> 0 = above random control genes).
# Tune after looking at the histograms; set a gate at the dip if a score is bimodal.
rna_thresholds <- list(
  osteo    = 0,     # osteo score
  hox2     = 0.1,   # hox PG2 must be actually expressed, not just background
  dorsal   = 0      # dorsal - ventral
)
arch2_margin <- 0   # hox2 - hox3 must exceed this (both modalities)

gate_vars <- c(osteo = "osteo", hox2 = "hox2", dorsal = "dorsal_score")

## ---- helpers ----------------------------------------------------------------
add_scores <- function(obj, assay, sets) {
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
  obj
}

apply_gates <- function(obj, thr) {
  obj$is_osteo  <- obj$osteo > thr$osteo
  obj$is_arch2  <- obj$hox2 > thr$hox2 & obj$arch2_score > arch2_margin
  obj$is_dorsal <- obj$dorsal_score > thr$dorsal
  obj$opercle_candidate <- obj$is_osteo & obj$is_arch2 & obj$is_dorsal
  obj@misc$gates <- thr
  obj
}

# fraction of cells passing each underlying gate variable in the RNA data
pass_fractions <- function(obj, thr) {
  sapply(names(gate_vars), function(k) mean(obj@meta.data[[gate_vars[k]]] > thr[[k]], na.rm = TRUE))
}

# thresholds that let the same fractions pass in another object
calibrated_thresholds <- function(obj, fracs) {
  lapply(setNames(names(gate_vars), names(gate_vars)), function(k)
    unname(quantile(obj@meta.data[[gate_vars[k]]], 1 - fracs[[k]], na.rm = TRUE)))
}

score_cols <- c("osteo", "maturation", "hox2", "arch2_score", "dorsal_score")

plot_scores <- function(obj, label, cluster_col) {
  p_feat <- FeaturePlot(obj, score_cols, reduction = "umap", order = TRUE, ncol = 2)
  cand <- colnames(obj)[obj$opercle_candidate %in% TRUE]
  p_cand <- DimPlot(obj, reduction = "umap", cells.highlight = cand,
                    cols.highlight = "firebrick", sizes.highlight = 0.6) +
    ggtitle(paste0(label, ": candidates (n = ", length(cand), ")")) + NoLegend()
  p_clu <- DimPlot(obj, reduction = "umap", group.by = cluster_col, label = TRUE) + NoLegend()
  ggsave(file.path(fig_dir, paste0(label, "_scores_umap.pdf")),
         p_feat / (p_clu | p_cand), width = 11, height = 18)

  thr <- obj@misc$gates
  md <- obj@meta.data[, unname(gate_vars)] |>
    pivot_longer(everything(), names_to = "score", values_to = "value")
  gates <- data.frame(score = unname(gate_vars), value = unlist(thr[names(gate_vars)]))
  p_hist <- ggplot(md, aes(value)) + geom_histogram(bins = 80) +
    geom_vline(data = gates, aes(xintercept = value), colour = "firebrick") +
    facet_wrap(~score, scales = "free") + theme_bw() + ggtitle(paste0(label, ": gates"))
  ggsave(file.path(fig_dir, paste0(label, "_score_histograms.pdf")), p_hist, width = 12, height = 4)
}

cluster_summary <- function(obj, cluster_col) {
  obj@meta.data |>
    group_by(cluster = .data[[cluster_col]]) |>
    summarise(n_cells = n(),
              across(all_of(score_cols), ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}"),
              frac_candidate = mean(opercle_candidate, na.rm = TRUE),
              .groups = "drop") |>
    arrange(desc(frac_candidate))
}

count_table <- function(obj, s, mod) {
  md <- obj@meta.data
  data.frame(
    stage = s, modality = mod,
    bio.rep = if ("bio.rep" %in% colnames(md)) as.character(md$bio.rep) else "all",
    is_osteo = md$is_osteo, is_arch2 = md$is_arch2,
    is_dorsal = md$is_dorsal, is_cand = md$opercle_candidate
  ) |>
    group_by(stage, modality, bio.rep) |>
    summarise(n_total       = n(),
              n_osteo       = sum(is_osteo, na.rm = TRUE),
              n_arch2       = sum(is_arch2, na.rm = TRUE),
              n_osteo_arch2 = sum(is_osteo & is_arch2, na.rm = TRUE),
              n_candidates  = sum(is_cand, na.rm = TRUE),
              .groups = "drop") |>
    mutate(pct_candidates = round(100 * n_candidates / n_total, 2))
}

## ---- run per stage ----------------------------------------------------------
counts <- list()
gates_used <- list()

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

  rna  <- add_scores(rna,  "RNA",      gene_sets)
  atac <- add_scores(atac, "ACTIVITY", gene_sets)

  rna   <- apply_gates(rna, rna_thresholds)
  fracs <- pass_fractions(rna, rna_thresholds)
  atac  <- apply_gates(atac, calibrated_thresholds(atac, fracs))

  gates_used[[s]] <- bind_rows(
    data.frame(stage = s, modality = "RNA",  gate = names(gate_vars),
               threshold = unlist(rna@misc$gates[names(gate_vars)]),  frac_pass = fracs),
    data.frame(stage = s, modality = "ATAC", gate = names(gate_vars),
               threshold = unlist(atac@misc$gates[names(gate_vars)]), frac_pass = fracs)
  )

  for (x in list(list(obj = rna,  mod = "RNA",  clu = rna_clu),
                 list(obj = atac, mod = "ATAC", clu = atac_clu))) {
    label <- paste0(s, "_", x$mod)
    plot_scores(x$obj, label, x$clu)
    write.csv(cluster_summary(x$obj, x$clu),
              file.path("results", paste0(label, "_opercle_cluster_scores.csv")), row.names = FALSE)
    counts[[label]] <- count_table(x$obj, s, x$mod)

    cand <- colnames(x$obj)[x$obj$opercle_candidate %in% TRUE]
    if (length(cand) > 0) {
      saveRDS(subset(x$obj, cells = cand),
              file.path("data/processed", paste0(label, "_opercle_osteo_candidates.rds")))
    }
  }

  cand <- colnames(atac)[atac$opercle_candidate %in% TRUE]
  if (length(cand) >= 10) {
    p_mat <- ggplot(atac@meta.data[cand, ], aes(maturation, dorsal_score, colour = arch2_score)) +
      geom_point() + scale_colour_viridis_c() + theme_bw() +
      ggtitle(paste0(s, " ATAC candidates: maturation vs dorsal identity"))
    ggsave(file.path(fig_dir, paste0(s, "_ATAC_candidates_maturation.pdf")), p_mat, width = 6, height = 5)
  }
}

counts_df <- bind_rows(counts)
gates_df  <- bind_rows(gates_used)
print(counts_df)
print(gates_df)
write.csv(counts_df, "results/opercle_candidate_counts.csv", row.names = FALSE)
write.csv(gates_df,  "results/opercle_gates_used.csv",       row.names = FALSE)
writeLines(capture.output(sessionInfo()), "results/sessionInfo_02_opercle_osteo_subset.txt")
