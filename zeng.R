
library(rhdf5)
library(Matrix)

ZENG_DIR <- "~/data/biccn_180521/Zeng"
ZENG_META <- file.path(ZENG_DIR, "BICCN_analysis_AIBS")

zeng_10x_cells_counts <- function() {
  return(zeng_10x_counts("10X_cells_MOp"))
}

zeng_10x_nuclei_counts <- function() {
  return(zeng_10x_counts("10x_nuclei_MOp"))
}

zeng_10x_counts <- function(ten_x_dir) {
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

zeng_smart_nuclei_counts <- function() {
  result <- read.table(file.path(ZENG_DIR, "SMARTer_nuclei_MOp", "exon_counts.csv.gz"), header = TRUE, sep = ",")
  rownames(result) <- result$sample_id
  result <- subset(result, select = -sample_id)
  return(t(Matrix(as.matrix(result), sparse = TRUE)))
}

zeng_metadata <- function(subdir) {
  clusters <- read.table(
    file.path(ZENG_META, subdir, "cluster.membership.csv"),
    sep = ",", header = TRUE
  )
  colnames(clusters) <- c("cell_id", "cluster_id")
  clusters$cell_id <- simplify_barcode_suffix(clusters$cell_id)
  cluster_info <- read.table(
    file.path(ZENG_META, subdir, "cluster.annotation.csv"),
    sep = ",", header = TRUE
  )
  cluster_info <- cluster_info[cluster_info$category != "Noise", ]
  cluster_id_match <- match(clusters$cluster_id, cluster_info$cluster_id)
  clusters$cluster_label <- cluster_info$cluster_label[cluster_id_match]
  clusters <- clusters[!is.na(clusters$cluster_label), ]
  return(clusters)
}

simplify_barcode_suffix <- function(barcodes) {
  suffix <- substring(barcodes, 18)
  suffix_values <- unique(suffix)
  for (i in seq_along(suffix_values)) {
    suffix[suffix == suffix_values[i]] <- i
  }
  return(paste0(substring(barcodes, 1, 17), suffix))
}
