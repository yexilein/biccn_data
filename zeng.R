
library(rhdf5)
library(Matrix)

ZENG_DIR <- "/home/fischer/data/biccn/raw_data/Zeng"

zeng_10x_counts <- function(ten_x_dir) {
  h5_data <- h5read(file.path(ZENG_DIR, ten_x_dir, "umi_counts.h5"), "/")[[1]]
  result <- sparseMatrix(i = h5_data$indices + 1,
                         j = column_indices(h5_data$indptr),
                         x = h5_data$data)
  rownames(result) <- h5_data$gene
  colnames(result) <- h5_data$barcodes
  return(result)
}

column_indices <- function(col_limits) {
  elements_per_col <- col_limits[-1] - col_limits[-length(col_limits)]
  return(rep(seq_along(elements_per_col), elements_per_col))
}

zeng_smart_counts <- function(subdir) {
  result <- read.table(file.path(ZENG_DIR, subdir, "exon.counts.csv.gz"), header = TRUE, sep = ",", row.names = 1)
  return(Matrix(as.matrix(result), sparse = TRUE))
}

zeng_clusters <- function(subdir) {
  clusters <- read.table(
    file.path(ZENG_DIR, subdir, "cluster.membership.csv"),
    sep = ",", header = TRUE
  )
  colnames(clusters) <- c("cell_id", "cluster_id")
  cluster_info <- read.table(
    file.path(ZENG_DIR, subdir, "cluster.annotation.csv"),
    sep = ",", header = TRUE
  )
  cluster_id_match <- match(clusters$cluster_id, cluster_info$cluster_id)
  clusters <- cbind(
    clusters,
    cluster_info[cluster_id_match, c('cluster_label', 'subclass_label', 'class_label')]
  )
  clusters <- clusters[!is.na(clusters$cluster_label), ]
  return(clusters)
}
