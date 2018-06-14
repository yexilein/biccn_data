library(SummarizedExperiment)
source("zeng.R")
source("/home/fischer/code/metaneighbor_v2/variableGene.R")

create_dataset <- function() {
  return(fuse_datasets(list(
    zeng_10x_nuclei(), zeng_smart_nuclei(), macosko_10x(), regev_10x()
  )))
}

zeng_10x_cells <- function() {
  return(counts_with_metadata(readRDS("zeng_10x_cells.rds"),
                              zeng_metadata("10X_cells_AIBS"),
			      "zeng_10x_cells"))
}

zeng_10x_nuclei <- function() {
  return(counts_with_metadata(zeng_10x_nuclei_counts(),
                              zeng_metadata("10X_nuclei_AIBS"),
			      "zeng_10x_nuclei"))
}

zeng_smart_nuclei <- function() {
  return(counts_with_metadata(filter_10x_genes(readRDS("zeng_smart_nuclei.rds")),
                              zeng_metadata("SmartSeq_nuclei_AIBS"),
			      "zeng_smart_nuclei"))
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

counts_with_metadata <- function(counts, notes, study_id) {
  colnames(counts) <- simplify_barcode_suffix(colnames(counts))
  cell_id_match <- match(colnames(counts), notes$cell_id)
  counts <- counts[, !is.na(cell_id_match)]
  labels <- notes$cluster_label[cell_id_match[!is.na(cell_id_match)]]
  return(list(data = counts, cell_type = labels,
              study_id = rep(study_id, length(labels))))
}

filter_10x_genes <- function(matrix_) {
  gene_mapping <- readRDS("sara_mapping.rds")
  match_10x_genes <- match(gene_mapping$name, rownames(matrix_))
  unfound_genes <- is.na(match_10x_genes)
  match_10x_genes[unfound_genes] <- 1
  result <- matrix_[match_10x_genes,]
  result[unfound_genes,] <- 0
  rownames(result) <- gene_mapping$ensembl
  return(result)
}

fuse_datasets <- function(datasets) {
  fused_data <- do.call(cbind, lapply(datasets, function(df) df$data))
  study_id <- do.call(c, lapply(datasets, function(df) df$study_id))
  cell_type <- do.call(c, lapply(datasets,
                                 function(df) as.character(df$cell_type)))
  colnames(fused_data) <- paste(study_id, cell_type, sep = "|")
  return(list(data = fused_data, study_id = study_id, cell_type = cell_type))
}

filter_variable_genes <- function(dataset) {
  centroids <- compute_centroids(dataset$data)
  rownames(centroids) <- as.character(rownames(centroids))
  variable_genes <- variableGenes(SummarizedExperiment(centroids), exp_labels=get_study_id(colnames(centroids)))
  dataset$data <- dataset$data[variable_genes, ]
  return(dataset)
}

