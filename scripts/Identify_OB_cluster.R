library(Seurat)
library(Signac)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Matrix)
library(tidyr)

setwd("~/Desktop/Hayden/git_repos/sox10_sc_analysis")
fig_dir <- "results/figures/OB_cluster"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
set.seed(1)

stages <-c("3dpf", "5dpf")

# read in data
rna_file  <- function(s) file.path("data/raw", paste0(s, "_scRNA.rds"))
atac_file <- function(s) file.path("data/raw", paste0(s, "_snATAC_SnapATAC.rds"))

ob_genes   <- c("runx2a", "runx2b", "sp7", "col10a1a", "ifitm5", "abcc6a")
cart_genes <- c("col2a1a", "sox9a", "matn1", "col9a1a")
dot_genes  <- c(ob_genes, "nfatc1", cart_genes,
                "hoxa2b", "hoxb2a", "hoxb3a", "pou3f3b", "hand2")

ob_z_min <- 1.5   # suggest clusters whose OB index is >= 1.5 SD above the cluster mean

# override after inspecting the tables/plots, e.g. list("3dpf" = c("14"), "5dpf" = c("3", "9"))
manual_rna_ob  <- list()
manual_atac_ob <- list()

## ---- helpers ----------------------------------------------------------------
zs <- function(x) { s <- sd(x, na.rm = TRUE); if (is.na(s) || s == 0) x * 0 else (x - mean(x, na.rm = TRUE)) / s }

group_means <- function(mat, grp) {
  m <- sapply(split(seq_len(ncol(mat)), grp), function(i) Matrix::rowMeans(mat[, i, drop = FALSE]))
  if (is.null(dim(m))) m <- matrix(m, nrow = 1, dimnames = list(rownames(mat), names(m)))
  m
}

rank_rna <- function(rna, clu) {
  cnt <- GetAssayData(rna, assay = "RNA", layer = "counts")
  grp <- as.character(rna@meta.data[[clu]])
  idx <- function(genes) {
    g <- intersect(genes, rownames(cnt))
    colMeans(group_means((cnt[g, , drop = FALSE] > 0) * 1, grp))
  }
  ob <- idx(ob_genes); ca <- idx(cart_genes)
  data.frame(cluster = names(ob), n_cells = as.integer(table(grp)[names(ob)]),
             ob_index = ob, cart_index = ca[names(ob)]) |>
    mutate(ob_z = zs(ob_index), cart_z = zs(cart_index)) |>
    arrange(desc(ob_index))
}

rank_atac <- function(atac, clu) {
  dat <- GetAssayData(atac, assay = "ACTIVITY", layer = "data")
  grp <- as.character(atac@meta.data[[clu]])
  idx <- function(genes) {
    g <- intersect(genes, rownames(dat))
    m <- group_means(dat[g, , drop = FALSE], grp)          # genes x clusters
    colMeans(t(apply(m, 1, zs)), na.rm = TRUE)             # mean z across genes
  }
  ob <- idx(ob_genes); ca <- idx(cart_genes)
  data.frame(cluster = names(ob), n_cells = as.integer(table(grp)[names(ob)]),
             ob_index = ob, cart_index = ca[names(ob)]) |>
    mutate(ob_z = zs(ob_index), cart_z = zs(cart_index)) |>
    arrange(desc(ob_index))
}

pick_clusters <- function(tab, manual) {
  sel <- if (!is.null(manual)) as.character(manual) else tab$cluster[tab$ob_z >= ob_z_min]
  cart_flag <- tab$cluster[tab$cluster %in% sel & tab$cart_z >= 1]
  if (length(cart_flag)) message("  ! selected clusters also high for cartilage: ",
                                 paste(cart_flag, collapse = ", "), " (check dot plot)")
  sel
}

rep_counts <- function(obj) {
  if ("bio.rep" %in% colnames(obj@meta.data)) table(bio.rep = obj$bio.rep) else table(all = rep("all", ncol(obj)))
}

## ---- run per stage ----------------------------------------------------------
summary_rows <- list()

for (s in stages) {
  if (!file.exists(rna_file(s)) || !file.exists(atac_file(s))) {
    message("Skipping ", s, ": files not found"); next
  }
  message("\n================ ", s, " ================")
  
  rna  <- UpdateSeuratObject(readRDS(rna_file(s)))
  atac <- UpdateSeuratObject(readRDS(atac_file(s)))
  DefaultAssay(rna)  <- "RNA"
  DefaultAssay(atac) <- "ACTIVITY"
  rclu <- grep("^RNA_snn_res",  colnames(rna@meta.data),  value = TRUE)[1]
  aclu <- grep("^ATAC_snn_res", colnames(atac@meta.data), value = TRUE)[1]
  
  ## rank and select
  rtab <- rank_rna(rna, rclu);   message("RNA clusters:");  print(rtab, row.names = FALSE)
  atab <- rank_atac(atac, aclu); message("ATAC clusters:"); print(atab, row.names = FALSE)
  write.csv(rtab, file.path("results", paste0(s, "_RNA_ob_cluster_ranking.csv")),  row.names = FALSE)
  write.csv(atab, file.path("results", paste0(s, "_ATAC_ob_cluster_ranking.csv")), row.names = FALSE)
  
  rsel <- pick_clusters(rtab, manual_rna_ob[[s]])
  asel <- pick_clusters(atab, manual_atac_ob[[s]])
  message("RNA OB clusters:  ", paste(rsel, collapse = ", "))
  message("ATAC OB clusters: ", paste(asel, collapse = ", "))
  
  rna$ob_cluster  <- as.character(rna@meta.data[[rclu]])  %in% rsel
  atac$ob_cluster <- as.character(atac@meta.data[[aclu]]) %in% asel
  
  ## dot plots + UMAP highlight
  dg_r <- intersect(dot_genes, rownames(rna))
  dg_a <- intersect(dot_genes, rownames(atac[["ACTIVITY"]]))
  p_r <- DotPlot(rna, features = dg_r, group.by = rclu) + RotatedAxis() +
    ggtitle(paste(s, "RNA: % detecting (size), mean expression (colour)"))
  p_a <- DotPlot(atac, features = dg_a, group.by = aclu) + RotatedAxis() +
    ggtitle(paste(s, "ATAC imputed: read colour only (sizes inflated by smoothing)"))
  u_r <- DimPlot(rna,  cells.highlight = colnames(rna)[rna$ob_cluster],
                 cols.highlight = "firebrick", sizes.highlight = 0.4) + NoLegend() + ggtitle("RNA OB clusters")
  u_a <- DimPlot(atac, cells.highlight = colnames(atac)[atac$ob_cluster],
                 cols.highlight = "firebrick", sizes.highlight = 0.4) + NoLegend() + ggtitle("ATAC OB clusters")
  ggsave(file.path(fig_dir, paste0(s, "_ob_clusters.pdf")),
         (u_r | u_a) / p_r / p_a, width = 14, height = 18)
  
  ## RNA-vs-ATAC cluster correspondence on RNA variable genes
  hvg <- intersect(VariableFeatures(rna), rownames(atac[["ACTIVITY"]]))
  if (length(hvg) == 0) hvg <- intersect(VariableFeatures(FindVariableFeatures(rna, nfeatures = 2000)),
                                         rownames(atac[["ACTIVITY"]]))
  rm_ <- group_means(GetAssayData(rna,  assay = "RNA",      layer = "data")[hvg, ], as.character(rna@meta.data[[rclu]]))
  am_ <- group_means(GetAssayData(atac, assay = "ACTIVITY", layer = "data")[hvg, ], as.character(atac@meta.data[[aclu]]))
  rz <- t(apply(rm_, 1, zs)); az <- t(apply(am_, 1, zs))   # z-score genes across clusters
  keep <- rowSums(is.na(rz)) == 0 & rowSums(is.na(az)) == 0
  cc <- cor(rz[keep, ], az[keep, ])
  cc_df <- as.data.frame(as.table(cc)) |> setNames(c("rna_cluster", "atac_cluster", "r")) |>
    mutate(rna_ob = rna_cluster %in% rsel, atac_ob = atac_cluster %in% asel)
  p_cc <- ggplot(cc_df, aes(atac_cluster, rna_cluster, fill = r)) + geom_tile() +
    geom_tile(data = subset(cc_df, rna_ob & atac_ob), colour = "black", linewidth = 0.8, fill = NA) +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick") +
    theme_bw() + ggtitle(paste(s, "RNA vs ATAC cluster correspondence (OB pairs outlined)"))
  ggsave(file.path(fig_dir, paste0(s, "_rna_atac_correspondence.pdf")), p_cc, width = 9, height = 8)
  for (r in rsel) {
    best <- names(sort(cc[r, ], decreasing = TRUE))[1:min(3, ncol(cc))]
    message("  RNA OB cluster ", r, " best ATAC matches: ", paste(best, collapse = ", "),
            if (!any(best %in% asel)) "  <- none is an ATAC OB cluster" else "")
  }
  
  ## save subsets + counts
  rob <- subset(rna,  cells = colnames(rna)[rna$ob_cluster])
  aob <- subset(atac, cells = colnames(atac)[atac$ob_cluster])
  saveRDS(rob, file.path("data/processed", paste0(s, "_RNA_OB.rds")))
  saveRDS(aob, file.path("data/processed", paste0(s, "_ATAC_OB.rds")))
  
  message("RNA OB cells per replicate:");  print(rep_counts(rob))
  message("ATAC OB cells per replicate:"); print(rep_counts(aob))
  summary_rows[[s]] <- data.frame(stage = s,
                                  rna_clusters = paste(rsel, collapse = ";"),  rna_cells = ncol(rob),
                                  atac_clusters = paste(asel, collapse = ";"), atac_cells = ncol(aob))
}

ob_summary <- bind_rows(summary_rows)
print(ob_summary)
write.csv(ob_summary, "results/ob_cluster_summary.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "results/sessionInfo_04_define_ob_clusters.txt")













