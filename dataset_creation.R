library(MetaNeighbor)
library(SingleCellExperiment)

source("ecker.R")
source("zeng.R")

create_zeng_10x_cells <- function() {
  data_ <- match_count_and_cluster_info(zeng_10x_counts("10X_cells_MOp"),
                                        zeng_clusters("10X_cells_MOp"))
  return(create_sce(data_$counts, data_$clusters, "zeng_10x_cells"))
}

match_count_and_cluster_info <- function(counts, clusters, simplify_barcode = FALSE) {
  if (simplify_barcode) {
    colnames(counts) <- simplify_barcode(colnames(counts))
    clusters$cell_id <- simplify_barcode(clusters$cell_id)
  }
  cell_id_match <- match(colnames(counts), clusters$cell_id)
  counts <- counts[, !is.na(cell_id_match)]
  clusters <- clusters[cell_id_match[!is.na(cell_id_match)], colnames(clusters) != 'cell_id']
  return(list(counts = counts, clusters = clusters))
}

simplify_barcode_suffix <- function(barcodes) {
  suffix <- substring(barcodes, 18)
  suffix_values <- unique(suffix)
  for (i in seq_along(suffix_values)) {
    suffix[suffix == suffix_values[i]] <- i
  }
  return(paste0(substring(barcodes, 1, 17), suffix))
}

create_sce <- function(counts, clusters, study_id) {
  rownames(counts) <- as.character(rownames(counts))
  colnames(counts) <- as.character(colnames(counts))
  result <- SingleCellExperiment(counts)
  result@colData <- cbind(result@colData, clusters)
  result$study_id <- rep(study_id, ncol(counts))
  return(result)
}

create_zeng_10x_nuclei <- function() {
  data_ <- match_count_and_cluster_info(zeng_10x_counts("10x_nuclei_MOp"),
                                        zeng_clusters("10x_nuclei_MOp"))
  return(create_sce(data_$counts, data_$clusters, "zeng_10x_nuclei"))
}

create_zeng_smart_nuclei <- function() {
  counts <- zeng_smart_counts("SMARTer_nuclei_MOp")
  colnames(counts) <- gsub("\\.", "-", colnames(counts))
  clusters <- zeng_clusters("SMARTer_nuclei_MOp")
  data_ <- match_count_and_cluster_info(counts, clusters)
  return(create_sce(data_$counts, data_$clusters, "zeng_smart_nuclei"))
}

create_zeng_smart_cells <- function() {
  counts <- zeng_smart_counts("SMARTer_cells_MOp")
  colnames(counts) <- gsub("\\.", "-", colnames(counts))
  clusters <- zeng_clusters("SMARTer_cells_MOp")
  data_ <- match_count_and_cluster_info(counts, clusters)
  return(create_sce(data_$counts, data_$clusters, "zeng_smart_cells"))
}

filter_10x_genes <- function(dataset, ensembl=TRUE) {
  gene_mapping <- readRDS("sara_mapping.rds")
  gene_ids <- if (ensembl) gene_mapping$ensembl else gene_mapping$name
  found_genes <- gene_ids %in% rownames(dataset)
  result <- dataset[as.character(gene_ids)[found_genes],]
  rownames(result) <- as.character(gene_mapping$ensembl)[found_genes]
  return(result)
}

macosko_10x <- function() {
  return(counts_with_metadata(readRDS("macosko.rds"),
                              zeng_metadata("10X_nuclei_Macosko"),
                              "macosko_10x"))
}

regev_10x <- function() {
  return(counts_with_metadata(readRDS("regev_10x.rds"),
                              zeng_metadata("10X_nuclei_Regev"),
                              "regev_10x"))
}

ecker_atac <- function() {
  return(counts_with_metadata(
    filter_10x_genes(readRDS("ecker_atac_gene_counts.rds")), ensembl = FALSE,
    ecker_atac_metadata(level = 2), "ecker_atac", simplify_barcode = FALSE
  ))
}

ecker_snmc <- function() {
  return(counts_with_metadata(
    filter_10x_genes(-ecker_snmc_gene_counts()),
    ecker_snmc_metadata(level = 1), "ecker_snmc", simplify_barcode = FALSE
  ))
}

fuse_datasets <- function(datasets) {
  counts <- do.call(cbind, lapply(datasets, function(df) assay(df)))
  col_data <- do.call(rbind, lapply(datasets, function(df) colData(df)))
  cell_names <- paste(col_data$study_id, rownames(col_data), sep = "_")
  colnames(counts) <- cell_names
  rownames(col_data) <- cell_names
  result <- SingleCellExperiment(counts, colData = col_data)
  return(result)
}

sample_based_variable_genes <- function(dataset) {
  dat <- dataset[, sample.int(ncol(dataset), 10000)]
  return(variableGenes(dat, exp_labels = dat$study_id))
}
