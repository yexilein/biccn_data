
ZENG_DIR <- "~/data/biccn_180521/Zeng/BICCN_analysis_AIBS/"

read_zeng_annotation <- function(subdir) {
  clusters <- read.table(
    file.path(ZENG_DIR, subdir, "cluster.membership.csv"),
    sep = ",", header = TRUE
  )
  colnames(clusters) <- c("cell_id", "cluster_id")
  clusters$cell_id <- simplify_barcode_suffix(clusters$cell_id)
  cluster_info <- read.table(
    file.path(ZENG_DIR, subdir, "cluster.annotation.csv"),
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