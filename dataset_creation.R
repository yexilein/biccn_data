
library(SingleCellExperiment)

source("ecker.R", local=TRUE)
source("zeng.R", local=TRUE)
#source("macosko_regev.R", local=TRUE) ## TODO update this file

create_zeng_datasets = function() {
    saveRDS(create_zeng_10x("10X_cells_v2_MOp", "zeng_10x_cells_v2"), "zeng_10x_cells_v2.rds")
    saveRDS(create_zeng_10x("10x_cells_v3_MOp", "zeng_10x_cells_v3"), "zeng_10x_cells_v3.rds")
    saveRDS(create_zeng_10x("10x_nuclei_v2_MOp", "zeng_10x_nuclei_v2"), "zeng_10x_nuclei_v2.rds")
    saveRDS(create_zeng_10x("10x_nuclei_v3_MOp", "zeng_10x_nuclei_v3"), "zeng_10x_nuclei_v3.rds")
    saveRDS(create_zeng_smart("SMARTer_cells_MOp", "zeng_smart_cells"), "zeng_smart_cells.rds")
    saveRDS(create_zeng_smart("SMARTer_nuclei_MOp", "zeng_smart_nuclei"), "zeng_smart_nuclei.rds")
}

create_zeng_10x = function(input_dir, study_id) {
  data_ = match_count_and_cluster_info(zeng_10x_counts(input_dir), zeng_clusters(input_dir))
  return(create_sce(data_$counts, data_$clusters, study_id))
}

match_count_and_cluster_info = function(counts, clusters) {
  cell_id_match = match(colnames(counts), clusters$cell_id)
  clusters = clusters[cell_id_match, colnames(clusters) != 'cell_id', drop=FALSE]
  return(list(counts = counts, clusters = clusters))
}

create_sce = function(counts, clusters, study_id) {
  rownames(counts) = as.character(rownames(counts))
  colnames(counts) = as.character(colnames(counts))
  result = SingleCellExperiment(list(counts = counts))
  result@colData = cbind(result@colData, clusters)
  result$study_id = rep(study_id, ncol(counts))
  return(result)
}

create_zeng_smart = function(input_dir, study_id) {
  data_ = match_count_and_cluster_info(zeng_smart_counts(input_dir), zeng_clusters(input_dir))
  return(create_sce(data_$counts, data_$clusters, study_id))
}

create_macosko_10x = function() {
  counts = macosko_10x_counts()
  colnames(counts) = paste0(colnames(counts), substring(colnames(counts), 18))
  clusters = zeng_clusters("BICCN_analysis_AIBS/10X_nuclei_Macosko",
                            c("cluster_label", "category_label"))
  clusters = dplyr::rename(clusters, class_label = category_label)
  clusters$subclass_label = clusters$cluster_label
  data_ = match_count_and_cluster_info(counts, clusters)
  return(create_sce(data_$counts, data_$clusters, "macosko_10x"))
}

create_macosko_10x_scportal = function() {
  counts = read.table("~/data/biccn/scportal/Macosko_snRNA_MC_ExpressionFile.txt.gz", header=TRUE, row.names=1)
  counts = Matrix(as.matrix(counts, sparse=TRUE))
  clusters = read.table("~/data/biccn/scportal/Macosko_snRNA_MC_metadat.txt", header=TRUE, skip=1)
  clusters = dplyr::rename(clusters, cluster_label = group, cell_id = TYPE)
  data_ = match_count_and_cluster_info(counts, clusters)
  return(create_sce(data_$counts, data_$clusters, "macosko_10x_scportal"))
}

create_regev_10x = function() {
  counts = regev_10x_counts()
  clusters = zeng_clusters("BICCN_analysis_AIBS/10X_nuclei_Regev",
                            c("cluster_label", "category_label"))
  clusters = dplyr::rename(clusters, class_label = category_label)
  clusters$subclass_label = clusters$cluster_label
  data_ = match_count_and_cluster_info(counts, clusters)
  return(create_sce(data_$counts, data_$clusters, "regev_10x"))
}

create_regev_10x_scportal = function() {
  counts = read.table("~/data/biccn/scportal/Regev_snRNA_MC_ExpressionFile.txt.gz", header=TRUE, row.names=1)
  counts = Matrix(as.matrix(counts, sparse=TRUE))
  clusters = read.table("~/data/biccn/scportal/Regev_snRNA_MC_metadat.txt", header=TRUE, skip=1)
  clusters = dplyr::rename(clusters, cluster_label = group, cell_id = TYPE)
  data_ = match_count_and_cluster_info(counts, clusters)
  return(create_sce(data_$counts, data_$clusters, "regev_10x_scportal"))
}

create_ecker_atac = function() {
    data_ = match_count_and_cluster_info(ecker_atac_gene_counts(),
                                         ecker_atac_metadata())
    result = create_sce(data_$counts, data_$clusters, "ecker_atac")
    assay(result, "fpkm") = ecker_atac_gene_fpkm()
    return(result)
}

create_ecker_snmc = function() {
    data_ = match_count_and_cluster_info(ecker_snmc_gene_counts(),
                                         ecker_snmc_metadata())
    return(create_sce(data_$counts, data_$clusters, "ecker_snmc"))
}
