library(SnapATAC)
library(Seurat)
library(Signac)
sc_data_3dpf <- readRDS("/Users/msbl/Desktop/Hayden/git_repos/sox10_sc_analysis/data/raw/3dpf_scRNA.rds")
atac_data_3dpf <- readRDS("/Users/msbl/Desktop/Hayden/git_repos/sox10_sc_analysis/data/raw/3dpf_snATAC_SnapATAC.rds")
inspect <- function(obj) {
  cat("class:", class(obj), "\n")
  if (inherits(obj, "Seurat")) {
    obj <- UpdateSeuratObject(obj)   # old Seurat v3 objects need this
    print(obj)
    cat("assays:", Assays(obj), "\n")
    cat("reductions:", Reductions(obj), "\n")
    print(head(obj@meta.data, 3))
  } else if (inherits(obj, "snap")) {
    print(obj)                       # shows bmat/pmat/gmat dimensions
    cat("embedding dims:", dim(obj@smat@dmat), "\n")
  }
  invisible(obj)
}
rna  <- inspect(sc_data_3dpf)
atac <- inspect(atac_data_3dpf)
