library(rhdf5)
library(Matrix)

BICCN_DIR <- "/home/fischer/data/biccn/raw_data_180521"
ATAC_DIR <- file.path(BICCN_DIR, "EckerCallaway_snATAC-Seq/v2_20180525")
SNMC_DIR <- file.path(BICCN_DIR, "EckerCallaway_snmC-Seq/v2_20180521")

ecker_atac_gene_counts <- function() {
  result <- read.table(file.path(ATAC_DIR, "gene_count_table", "MOp_raw_gene_count_table.txt.gz"), header=TRUE)
  rownames(result) <- result$barcode
  result <- subset(result, select = -barcode)
  return(t(Matrix(as.matrix(result), sparse = TRUE)))
}

ecker_atac_metadata <- function(level=4) {
  metadata <- read.table(file.path(ATAC_DIR, "MOp_multi_level_cluster.csv"), sep = ",", header = TRUE)
  return(data.frame(cell_id = metadata$barcode, cluster_label = paste0("zzz", metadata[, level+1])))
}

ecker_snmc_gene_counts <- function() {
  result <- cbind(
    parse_ecker_snmc_counts(
      file.path(SNMC_DIR, "5. count_table", "3C_gene.mCH.h5"), "/gene_2473_53379_mc"
    ),
    parse_ecker_snmc_counts(
      file.path(SNMC_DIR, "5. count_table", "4B_gene.mCH.h5"), "/gene_2881_53379_mc"
    )
  )
  return(result)
}

parse_ecker_snmc_counts <- function(filename, h5_path) {
  h5_data <- h5read(filename, h5_path)
  result <- Matrix(h5_data$block0_values, sparse=TRUE)
  gene_names <- h5_data$axis0_level0[h5_data$axis0_label0+1]
  rownames(result) <- sapply(strsplit(gene_names, "\\."), head, 1)
  colnames(result) <- h5_data$axis1
  return(result)
}

ecker_snmc_metadata <- function(level=7) {
  metadata <- read.table(file.path(SNMC_DIR, "4. cluster_martix.csv"), sep = ",", header = TRUE)
  return(data.frame(cell_id = metadata$cell_id, cluster_label = metadata[, level+1]))
}
