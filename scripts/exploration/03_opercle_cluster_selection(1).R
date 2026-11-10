## 03_opercle_cluster_selection.R
## Choose opercle osteoblast-lineage cells at the cluster level instead of per-cell gates.
##
## RNA : subset osteo-lineage clusters -> recluster -> % of cells detecting positional
##       genes per subcluster (robust to dropout at the group level).
## ATAC: subset osteo-rich clusters -> recluster on LSI -> imputed positional expression
##       and chromVAR motif deviations (HOX/PBX/MEIS, POU3F3, HAND, DLX...) per subcluster.
##       Chromatin can retain positional marks after the TF stops being transcribed.
## Atlas: optional module score from opercle markers in your voxel atlas.
##
## Target profile: osteo-high, hox PG2 detected, hox PG3 ~0, pou3f3b / dorsal motifs up,
## HAND/DLX (ventral) low.

suppressPackageStartupMessages({
  library(Seurat); library(Signac); library(Matrix)
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(JASPAR2020); library(TFBSTools)
})

setwd("~/Desktop/Hayden/git_repos/sox10_sc_analysis")
fig_dir <- "results/figures/opercle_clusters"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(1)

## ---- parameters -------------------------------------------------------------
stages <- c("3dpf", "5dpf")
rna_file  <- function(s) file.path("data/raw", paste0(s, "_scRNA.rds"))
atac_file <- function(s) file.path("data/raw", paste0(s, "_snATAC_SnapATAC.rds"))

osteo_genes       <- c("runx2a", "runx2b", "sp7")
osteo_late_genes  <- c("sp7", "col10a1a", "ifitm5", "spp1")
osteo_select_mode <- "cells"   # "cells": cells detecting any osteo_genes (use when osteoblasts
                               #          don't form their own cluster, e.g. 3 dpf)
                               # "clusters": whole clusters with >= osteo_cluster_min detecting
osteo_cluster_min <- 0.30

# set these after looking at the printed tables, e.g. list("3dpf" = c("4", "11"))
manual_rna_osteo_clusters  <- list()
manual_atac_osteo_clusters <- list()
atac_n_osteo_clusters      <- 3   # otherwise: top N ATAC clusters by imputed osteo score

subcluster_res <- 1.0
rna_pcs        <- 1:15
atac_lsi       <- 2:30

positional_genes <- c(
  "runx2a", "runx2b", "sp7", "col10a1a", "ifitm5",          # osteo lineage
  "hoxa2b", "hoxb2a",                                       # arch 2+
  "hoxa3a", "hoxb3a", "hoxd3a",                             # arch 3+
  "pou3f3b", "jag1b", "hey1",                               # dorsal
  "dlx5a", "dlx3b", "hand2",                                # intermediate / ventral
  "scxa", "tnmd", "nfatc1",                                 # attachment / enthesis
  "col2a1a", "sox9a"                                        # cartilage (contamination check)
)
motif_names <- c("HOXA2", "HOXB2", "PBX1", "PBX3", "MEIS1", "MEIS2", "POU3F3",
                 "HAND1", "HAND2", "DLX5", "DLX6", "PITX1", "NR5A2", "RUNX2", "SP7")

# optional: one gene per line, e.g. opercle-domain markers or top SRGs from your voxel atlas
atlas_marker_file <- "data/reference/opercle_atlas_markers.txt"

## ---- helpers ----------------------------------------------------------------
pct_detected <- function(obj, genes, group, assay = "RNA") {
  genes <- intersect(genes, rownames(obj[[assay]]))
  det <- GetAssayData(obj, assay = assay, layer = "counts")[genes, , drop = FALSE] > 0
  grp <- as.character(obj@meta.data[[group]])
  out <- sapply(split(seq_len(ncol(obj)), grp),
                function(i) Matrix::rowMeans(det[, i, drop = FALSE]))
  out <- as.data.frame(t(matrix(out, nrow = length(genes),
                                dimnames = list(genes, names(split(seq_len(ncol(obj)), grp))))))
  out$n_cells <- as.integer(table(grp)[rownames(out)])
  out$cluster <- rownames(out)
  out
}

mean_by_group <- function(mat, grp) {
  grp <- as.character(grp)
  sapply(split(seq_len(ncol(mat)), grp), function(i) Matrix::rowMeans(mat[, i, drop = FALSE]))
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

read_atlas_markers <- function() {
  if (!file.exists(atlas_marker_file)) return(NULL)
  g <- trimws(readLines(atlas_marker_file))
  g[nzchar(g)]
}

## ---- run per stage ----------------------------------------------------------
for (s in stages) {
  if (!file.exists(rna_file(s)) || !file.exists(atac_file(s))) {
    message("Skipping ", s, ": files not found"); next
  }
  message("\n================ ", s, " ================")

  ## ----- RNA -----
  rna <- UpdateSeuratObject(readRDS(rna_file(s)))
  DefaultAssay(rna) <- "RNA"
  clu <- grep("^RNA_snn_res", colnames(rna@meta.data), value = TRUE)[1]

  og  <- intersect(osteo_genes, rownames(rna))
  any_osteo <- Matrix::colSums(GetAssayData(rna, layer = "counts")[og, , drop = FALSE] > 0) > 0
  osteo_tab <- data.frame(cluster = rna@meta.data[[clu]], osteo = any_osteo) |>
    group_by(cluster) |>
    summarise(n_cells = n(), pct_osteo = mean(osteo), .groups = "drop") |>
    arrange(desc(pct_osteo))
  message("RNA clusters by % cells detecting ", paste(og, collapse = "/"), ":")
  print(osteo_tab, n = Inf)

  if (osteo_select_mode == "clusters") {
    sel <- manual_rna_osteo_clusters[[s]]
    if (is.null(sel)) sel <- as.character(osteo_tab$cluster[osteo_tab$pct_osteo >= osteo_cluster_min])
    message("RNA osteo clusters used: ", paste(sel, collapse = ", "))
    cells_use <- colnames(rna)[as.character(rna@meta.data[[clu]]) %in% sel]
  } else {
    cells_use <- colnames(rna)[any_osteo]
  }
  frac_osteo <- length(cells_use) / ncol(rna)
  message("RNA osteo-lineage cells selected: ", length(cells_use),
          " (", round(100 * frac_osteo, 2), "%)")
  if (length(cells_use) < 30) {
    message("Too few osteo-lineage cells at ", s, " to subcluster; skipping."); next
  }

  sub <- subset(rna, cells = cells_use)
  sub <- FindVariableFeatures(sub, nfeatures = 2000) |> ScaleData(verbose = FALSE) |>
    RunPCA(npcs = 30, verbose = FALSE)
  sub <- FindNeighbors(sub, dims = rna_pcs, verbose = FALSE) |>
    FindClusters(resolution = subcluster_res, verbose = FALSE) |>
    RunUMAP(dims = rna_pcs, verbose = FALSE)
  sub$subcluster <- as.character(Idents(sub))

  det_tab <- pct_detected(sub, positional_genes, "subcluster")
  max_of <- function(df, cols) {
    cols <- intersect(cols, colnames(df))
    if (length(cols) == 0) return(NA_real_)
    do.call(pmax, df[cols])
  }
  det_tab$hox2_pct <- max_of(det_tab, c("hoxa2b", "hoxb2a"))
  det_tab$hox3_pct <- max_of(det_tab, c("hoxa3a", "hoxb3a", "hoxd3a"))
  det_tab <- relocate(det_tab, cluster, n_cells, hox2_pct, hox3_pct) |> arrange(desc(hox2_pct))
  write.csv(det_tab, file.path("results", paste0(s, "_RNA_osteo_subcluster_detection.csv")), row.names = FALSE)

  rep_tab <- if ("bio.rep" %in% colnames(sub@meta.data)) table(sub$subcluster, sub$bio.rep) else NULL
  if (!is.null(rep_tab)) { message("Subcluster x replicate:"); print(rep_tab) }

  p1 <- DimPlot(sub, group.by = "subcluster", label = TRUE) + NoLegend() + ggtitle(paste(s, "RNA osteo subclusters"))
  p2 <- if ("bio.rep" %in% colnames(sub@meta.data)) DimPlot(sub, group.by = "bio.rep") else plot_spacer()
  p3 <- DotPlot(sub, features = intersect(positional_genes, rownames(sub)), group.by = "subcluster") +
    RotatedAxis() + ggtitle("% detecting (dot size) and mean expression")
  ggsave(file.path(fig_dir, paste0(s, "_RNA_osteo_subclusters.pdf")),
         (p1 | p2) / p3, width = 14, height = 11)

  atlas <- read_atlas_markers()
  if (!is.null(atlas)) {
    ag <- intersect(atlas, rownames(sub))
    message("Atlas markers found in RNA: ", length(ag), "/", length(atlas))
    if (length(ag) >= 3) {
      sub <- AddModuleScore(sub, features = list(ag), name = "atlas_opercle", seed = 1)
      p_at <- VlnPlot(sub, "atlas_opercle1", group.by = "subcluster", pt.size = 0) + NoLegend() +
        ggtitle("Opercle atlas marker score")
      ggsave(file.path(fig_dir, paste0(s, "_RNA_atlas_score.pdf")), p_at, width = 8, height = 4)
    }
  }
  saveRDS(sub, file.path("data/processed", paste0(s, "_RNA_osteo_subclustered.rds")))

  ## ----- ATAC -----
  atac <- UpdateSeuratObject(readRDS(atac_file(s)))
  aclu <- grep("^ATAC_snn_res", colnames(atac@meta.data), value = TRUE)[1]
  aog  <- intersect(osteo_genes, rownames(atac[["ACTIVITY"]]))   # core genes only; col10a1a/spp1 also mark cartilage
  atac <- AddModuleScore(atac, features = list(aog), name = "osteo_imp", assay = "ACTIVITY", seed = 1)

  atab <- atac@meta.data |>
    group_by(cluster = .data[[aclu]]) |>
    summarise(n_cells = n(), mean_osteo_imputed = mean(osteo_imp1), .groups = "drop") |>
    arrange(desc(mean_osteo_imputed))
  message("ATAC clusters by mean imputed osteo score:")
  print(atab, n = Inf)

  if (osteo_select_mode == "clusters") {
    asel <- manual_atac_osteo_clusters[[s]]
    if (is.null(asel)) asel <- as.character(head(atab$cluster, atac_n_osteo_clusters))
    message("ATAC osteo clusters used: ", paste(asel, collapse = ", "))
    acells <- colnames(atac)[as.character(atac@meta.data[[aclu]]) %in% asel]
  } else {
    # imputed expression is smoothed, so take the same fraction of cells as in RNA
    cut <- quantile(atac$osteo_imp1, 1 - frac_osteo)
    acells <- colnames(atac)[atac$osteo_imp1 >= cut]
  }
  message("ATAC osteo-lineage cells selected: ", length(acells))

  asub <- subset(atac, cells = acells)
  DefaultAssay(asub) <- "ATAC"
  asub <- FindNeighbors(asub, reduction = "lsi", dims = atac_lsi, verbose = FALSE) |>
    FindClusters(algorithm = 3, resolution = subcluster_res, verbose = FALSE) |>
    RunUMAP(reduction = "lsi", dims = atac_lsi, verbose = FALSE)
  asub$subcluster <- as.character(Idents(asub))

  # imputed positional expression per subcluster
  pg  <- intersect(positional_genes, rownames(asub[["ACTIVITY"]]))
  imp <- mean_by_group(GetAssayData(asub, assay = "ACTIVITY", layer = "data")[pg, , drop = FALSE],
                       asub$subcluster)
  write.csv(t(imp), file.path("results", paste0(s, "_ATAC_osteo_subcluster_imputed_expr.csv")))

  # chromVAR motif deviations per subcluster
  cv  <- GetAssayData(asub, assay = "chromvar", layer = "data")
  lab <- motif_label(rownames(cv))
  keep <- which(toupper(lab) %in% toupper(motif_names))
  if (length(keep) > 0) {
    mot <- mean_by_group(cv[keep, , drop = FALSE], asub$subcluster)
    rownames(mot) <- make.unique(lab[keep])
    write.csv(t(mot), file.path("results", paste0(s, "_ATAC_osteo_subcluster_motifs.csv")))
    mot_df <- as.data.frame(mot) |> mutate(motif = rownames(mot)) |>
      pivot_longer(-motif, names_to = "subcluster", values_to = "deviation")
    p_mot <- ggplot(mot_df, aes(subcluster, motif, fill = deviation)) + geom_tile() +
      scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick") +
      theme_bw() + ggtitle(paste(s, "chromVAR deviations, ATAC osteo subclusters"))
  } else {
    message("No chromVAR rows matched motif_names; check rownames(atac[['chromvar']])")
    p_mot <- plot_spacer()
  }

  imp_df <- as.data.frame(imp) |> mutate(gene = rownames(imp)) |>
    pivot_longer(-gene, names_to = "subcluster", values_to = "mean_imputed")
  p_imp <- ggplot(imp_df, aes(subcluster, factor(gene, levels = rev(pg)), fill = mean_imputed)) +
    geom_tile() + scale_fill_viridis_c() + theme_bw() + labs(y = NULL) +
    ggtitle("Mean imputed expression")
  pa <- DimPlot(asub, group.by = "subcluster", label = TRUE) + NoLegend() +
    ggtitle(paste(s, "ATAC osteo subclusters"))
  ggsave(file.path(fig_dir, paste0(s, "_ATAC_osteo_subclusters.pdf")),
         pa / (p_imp | p_mot), width = 14, height = 12)

  saveRDS(asub, file.path("data/processed", paste0(s, "_ATAC_osteo_subclustered.rds")))
}

writeLines(capture.output(sessionInfo()), "results/sessionInfo_03_opercle_cluster_selection.txt")
