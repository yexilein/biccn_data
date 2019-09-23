
library(rhdf5)
library(Matrix)

zeng_dir = function() {
    "/home/fischer/data/biccn/dropbox/BICCN MiniBrain JointDataAnalysis/Zeng"
}

joint_dir = function() {
    "/home/fischer/data/biccn/dropbox/BICCN MiniBrain JointDataAnalysis/Transcriptome_jointAnalysis"
}

zeng_10x_counts = function(ten_x_dir) {
  h5_data = h5read(file.path(ZENG_DIR, ten_x_dir, "umi_counts.h5"), "/")[[1]]
  gene_names = h5_data$features$name
  barcodes = zeng_short_to_long_barcode(ten_x_dir)[h5_data$barcodes]
  result = Matrix::sparseMatrix(i = h5_data$indices+1,
                                p = h5_data$indptr,
                                x = h5_data$data,
                                dim = c(length(gene_names), length(barcodes)),
                                dimnames = list(gene_names, barcodes))
  result = result[, !is.na(barcodes)]
  return(result)
}

zeng_short_to_long_barcode = function(subdir) {
  result = read.table(
    file.path(ZENG_DIR, subdir, "cluster.membership.csv"),
    sep = ",", header = TRUE, stringsAsFactors = FALSE
  )$X
  names(result) = sapply(strsplit(result, split = "L"), "[", 1)
  return(result)
}

zeng_smart_counts = function(subdir) {
  result = read.table(file.path(zeng_dir(), subdir, "exon.counts.csv.gz"),
                      header = TRUE, sep = ",", row.names = 1, check.names = FALSE)
  return(Matrix(as.matrix(result), sparse = TRUE))
}

macosko_v3_counts = function() {
    subdir = file.path(joint_dir(), "10X_nuclei_v3_Broad")
    result = Matrix::readMM(file.path(subdir, "matrix.mtx.gz"))
    rownames(result) = read.csv(file.path(subdir, "features.csv"), stringsAsFactors=FALSE)$gene
    colnames(result) = read.csv(file.path(subdir, "barcodes.csv"), stringsAsFactors=FALSE)$barcode
    return(result)
}

zeng_split_seq = function() {
  h5_data = h5read(file.path(zeng_dir(), "Splitseq_MOp", "raw_counts.h5"), "/")
  barcodes = sapply()
}

zeng_clusters = function(subdir, labels = c('cluster_label', 'subclass_label', 'class_label'), main_dir = zeng_dir()) {
    data_dir = file.path(main_dir, subdir)
    clusters = read.csv(file.path(data_dir, "cluster.membership.csv"), stringsAsFactors = FALSE)
    colnames(clusters) = c("cell_id", "cluster_id")
    cluster_info = read.csv(file.path(data_dir, "cluster.annotation.csv"), stringsAsFactors = FALSE)
    cluster_id_match = match(clusters$cluster_id, cluster_info$cluster_id)
    clusters = cbind(clusters, cluster_info[cluster_id_match, labels])
    clusters = clusters[!is.na(clusters$cluster_label), ]

    clusters$is_qc = clusters$cell_id %in% qc_cells(subdir)
    return(clusters)
}

qc_cells = function(subdir, main_dir = zeng_dir()) {
    qc_cells = read.csv(file.path(main_dir, subdir, "QC.csv"), stringsAsFactors = FALSE)$x
    return(qc_cells)
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
