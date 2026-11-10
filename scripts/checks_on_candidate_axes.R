library(dplyr); library(readr)
d <- "results/ob_axes/5dpf"

# 1. Technical/state contamination for candidate axes
read_csv(file.path(d, "axis_annotation_atac.csv")) |>
  filter(axis %in% paste0("LSIOB_", 2:6),
         annotation %in% c("cell_cycle", "stress", "log_depth", "maturation",
                           "hox2_vs_hox3", "dorsal_vs_ventral")) |>
  tidyr::pivot_wider(names_from = annotation, values_from = rho) |> print()

# 2. Genes in each coherent module
coherent <- data.frame(axis = paste0("LSIOB_", c(2, 3, 4, 5, 6)),
                       end  = c("negative", "negative", "positive", "positive", "negative"))
read_csv(file.path(d, "axis_module_genes.csv")) |>
  semi_join(coherent, by = c("axis", "end")) |>
  group_by(axis, end) |> summarise(genes = paste(gene, collapse = ", "), .groups = "drop") |>
  print(width = Inf)


library(Seurat); library(dplyr); library(readr)
rob <- readRDS("data/processed/5dpf_RNA_OB.rds"); DefaultAssay(rob) <- "RNA"
mods <- read_csv("results/ob_axes/5dpf/axis_module_genes.csv", show_col_types = FALSE)

for (ax in c("LSIOB_3", "LSIOB_5")) {
  end <- if (ax == "LSIOB_3") "negative" else "positive"
  g <- mods |> filter(axis == ax, end == !!end) |> pull(gene)
  rob <- AddModuleScore(rob, features = list(intersect(g, rownames(rob))), name = paste0(ax, "_"), seed = 1)
  sc <- rob@meta.data[[paste0(ax, "_1")]]
  rob@meta.data[[paste0(ax, "_group")]] <- ifelse(sc >= quantile(sc, 0.75), "module-high", "rest")
}

genes <- c("runx2a", "runx2b", "sp7", "col10a1a", "ifitm5", "abcc6a", "nfatc1",   # OB / landmarks
           "dcn", "slit3", "efemp2a", "fstb", "cebpd", "nfia", "mfap2",           # module genes
           "cxcl12a", "osr1", "pah", "hpdb", "krt4", "col2a1a")                     # stroma / dermal fib / epithelium / cartilage
genes <- intersect(genes, rownames(rob))
DotPlot(rob, features = genes, group.by = "LSIOB_3_group") + RotatedAxis() + ggtitle("LSIOB_3 module")
DotPlot(rob, features = genes, group.by = "LSIOB_5_group") + RotatedAxis() + ggtitle("LSIOB_5 module")
