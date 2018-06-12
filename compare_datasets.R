library(SummarizedExperiment)
source("annotations.R")
source("my_metaneighbor.R")

test_metaneighbor <- function(dataset = NULL) {
  if (is.null(dataset)) { dataset <- filter_genes(fused_data()) }
  runtime <- system.time(
    all_cells <- MetaNeighborUS(dataset$data,
				dataset$study_id,
				as.character(dataset$cell_type),
				n_centroids = 1)
  )
  print(runtime)
  return(all_cells)
}

test_plots <- function(results) {
  pdf('results/macosko_zeng_whole.pdf')
  plot_NV_heatmap((results + t(results))/2, as_dist = TRUE)
  dev.off()
  macosko_ind <- get_study_id(rownames(results)) == 1
  zeng_ind <- get_study_id(rownames(results)) == 2
  macosko_votes <- results[zeng_ind, macosko_ind]
  zeng_votes <- results[macosko_ind, zeng_ind]
  pdf('results/macosko_vs_zeng.pdf')
  plot_NV_heatmap(macosko_votes)
  dev.off()
  pdf('results/zeng_vs_macosko.pdf')
  plot_NV_heatmap(zeng_votes)
  dev.off()
  pdf('results/macosko_vs_zeng_summary.pdf')
  plot_NV_heatmap((macosko_votes + t(zeng_votes))/2)
  dev.off()
}

fused_data <- function() {
  zeng <- zeng_10x_cells()
  macosko <- macosko_10x()
  fused_data <- cbind(macosko$data, zeng$data)
  study_id <- rep(1:2, sapply(c(macosko$data, zeng$data), ncol))
  cell_type = c(as.character(macosko$cell_type),
                as.character(zeng$cell_type))
  colnames(fused_data) <- paste(study_id, cell_type, sep = "|")
  return(list(data = fused_data, study_id = study_id, cell_type = cell_type))
}

zeng_10x_cells <- function() {
  zeng <- readRDS("zeng.rds")
  colnames(zeng) <- simplify_barcode_suffix(colnames(zeng))
  notes <- read_zeng_annotation("10X_cells_AIBS")
  cell_id_match <- match(colnames(zeng), notes$cell_id)
  zeng <- zeng[, !is.na(cell_id_match)]
  return(list(data = zeng, cell_type = notes$cluster_label))
}

macosko_10x <- function() {
  macosko <- readRDS("macosko.rds")
  notes <- read_zeng_annotation("10X_nuclei_Macosko")
  cell_id_match <- match(colnames(macosko), notes$cell_id)
  macosko <- macosko[, !is.na(cell_id_match)]
  labels <- notes$cluster_label[cell_id_match[!is.na(cell_id_match)]]
  return(list(data = macosko, cell_type = labels))
}

filter_genes <- function(dataset) {
  centroids <- compute_centroids(dataset$data)
  variable_genes <- variableGenes(SummarizedExperiment(centroids), exp_labels=get_study_id(colnames(centroids)))
  dataset$data <- dataset$data[variable_genes, ]
  return(dataset)
}
