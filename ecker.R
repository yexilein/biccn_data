library(Matrix)

BICCN_DIR <- "/home/fischer/data/biccn/raw_data_180521"
ATAC_DIR <- file.path(BICCN_DIR, "EckerCallaway_snATAC-Seq/v2_20180525")

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
