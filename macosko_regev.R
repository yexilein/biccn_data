
library(cellrangerRkit)

MACOSKO_DIR <- "/home/fischer/data/biccn/raw_data_180521/Macosko_Regev"

macosko_10x <- function(data_dir = MACOSKO_DIR) {
  return(fuse_10x_dataset(file.path(data_dir, "macosko_motor_cortex_1")))
}

regev_10x <- function(data_dir = MACOSKO_DIR) {
  return(fuse_10x_dataset(file.path(data_dir, "regev_motor_cortex_2")))
}

fuse_10x_dataset <- function(data_dir) {
  subdirectories <- list.dirs(data_dir, full.names = TRUE, recursive = FALSE)
  gbm <- lapply(subdirectories,
                function(dir) cellrangerRkit::load_cellranger_matrix(dir))
  for (i in seq_along(gbm)) {
    colnames(gbm[[i]]) <- change_barcode_suffix(gbm[[i]], i)
  }
  check_rowname_consistency(gbm)
  count_matrix <- do.call(cbind, lapply(gbm, Biobase::exprs))
  return(count_matrix)
}

change_barcode_suffix <- function(dataset, suffix) {
  return(paste0(
      substring(colnames(dataset), 1, 17), suffix
  ))
}

check_rowname_consistency <- function(dataset_list) {
  ref <- rownames(dataset_list[[1]])
  is_consistent<- all(sapply(dataset_list,
                             function(d) all(rownames(d) == ref)))
  if (!is_consistent) {
    stop("Row names are not consistent across datasets.")
  }
}

if (!interactive()) {
  args = commandArgs(trailingOnly = TRUE)
  saveRDS(fuse_10x_dataset(args[1]), args[2])
}
