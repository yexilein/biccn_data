
library(rhdf5)
library(Matrix)

ZENG_DIR <- "~/data/biccn_180521/Zeng"

zeng_10x_cells <- function() {
  return(zeng_10x("10X_cells_MOp"))
}

zeng_10x_nuclei <- function() {
  return(zeng_10x("10x_nuclei_MOp"))
}

zeng_10x <- function(ten_x_dir) {
  h5_data <- h5read(file.path(ZENG_DIR, ten_x_dir, "umi_counts.h5"),
  	     	    "/mm10_mrna")
  result <- sparseMatrix(i = h5_data$indices + 1,
                         j = column_indices(h5_data$indptr),
			 x = h5_data$data)
  rownames(result) <- h5_data$genes
  colnames(result) <- h5_data$barcodes
  return(result)
}

column_indices <- function(col_limits) {
  elements_per_col <- col_limits[-1] - col_limits[-length(col_limits)]
  return(rep(seq_along(elements_per_col), elements_per_col))
}
