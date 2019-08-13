
library(rhdf5)
library(Matrix)

ZENG_DIR = "/home/fischer/data/biccn/raw_data/Zeng"

zeng_10x_counts = function(ten_x_dir) {
  h5_data = h5read(file.path(ZENG_DIR, ten_x_dir, "umi_counts.h5"), "/")[[1]]
  result = Matrix::sparseMatrix(i = h5_data$indices+1,
                                p = h5_data$indptr,
                                x = h5_data$data,
                                dim = c(length(h5_data$gene_names), length(h5_data$barcodes)),
                                dimnames = list(h5_data$gene_names, h5_data$barcodes))
  return(result)
}

column_indices = function(col_limits) {
  elements_per_col = col_limits[-1] - col_limits[-length(col_limits)]
  return(rep(seq_along(elements_per_col), elements_per_col))
}

zeng_smart_counts = function(subdir) {
  result = read.table(file.path(ZENG_DIR, subdir, "exon.counts.csv.gz"), header = TRUE, sep = ",", row.names = 1)
  return(Matrix(as.matrix(result), sparse = TRUE))
}

zeng_split_seq = function() {
  h5_data = h5read(file.path(ZENG_DIR, "Splitseq_MOp", "raw_counts.h5"), "/")
  barcodes = sapply()
}

zeng_clusters = function(subdir, labels = c('cluster_label', 'subclass_label', 'class_label')) {
  clusters = read.table(
    file.path(ZENG_DIR, subdir, "cluster.membership.csv"),
    sep = ",", header = TRUE
  )
  colnames(clusters) = c("cell_id", "cluster_id")
  cluster_info = read.table(
    file.path(ZENG_DIR, subdir, "cluster.annotation.csv"),
    sep = ",", header = TRUE
  )
  cluster_id_match = match(clusters$cluster_id, cluster_info$cluster_id)
  clusters = cbind(clusters, cluster_info[cluster_id_match, labels])
  clusters = clusters[!is.na(clusters$cluster_label), ]
  return(clusters)
}

convert_barcode_to_nucleotides = function(integer_barcode) {
  result = rep(0, 16)
  for (i in 1:16) {
    result[16-i] = integer_barcode %% 4
    integer_barcode = as.integer(integer_barcode %/% 4)
  }
  result[result == 0] = 'A'
  result[result == 1] = 'C'
  result[result == 2] = 'G'
  result[result == 3] = 'T'
  return(paste(result, collapse=""))
}
