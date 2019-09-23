
library(SingleCellExperiment)

current_dir_ = getwd()
common_dir_ = "~/projects/common"

setwd(common_dir_)
source("brain_datasets.R")
source("variable_genes.R")
setwd(current_dir_)
source("ecker.R", local=TRUE)
source("zeng.R", local=TRUE)
#source("macosko_regev.R", local=TRUE) ## TODO update this file


main = function() {
    create_combined_dataset()
}

create_combined_dataset = function() {
    dataset = lapply(set_names(biccn_datasets()), load_biccn_dataset)
    gc()
    dataset = restrict_to_common_genes(dataset)
    gc()
    dataset = fuse_datasets(dataset)
    saveRDS(dataset, "full_data.rds")
}

create_hvg_list = function() {
    dataset = readRDS("full_data.rds")
    hvg = variable_genes_by_dataset(dataset)
    recurrence = table(unlist(hvg))
    write(names(recurrence)[recurrence >= 6], "hvg_6.txt")
    write(names(recurrence)[recurrence >= 7], "hvg_7.txt")
    write(names(recurrence)[recurrence >= 8], "hvg_8.txt")
}

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

update_broad_dataset = function(input_file = "old/190813/macosko_10x.rds") {
    dataset = readRDS(input_file)
    rownames(dataset) = convert_to_mgi_symbols_from_10x(rownames(dataset))
    dataset$subclass_label = infer_subclass_labels(dataset)
    names(assays(dataset))[1] = "counts"
    result = SingleCellExperiment(assays = assays(dataset),
                                  colData = colData(dataset))
    return(result)
}

infer_subclass_labels = function(dataset) {
    # Special rules
    # L5 NP -> L5/6 NP
    # Lamp5 Sncg -> Sncg
    exceptions = c("L5 NP", "L5 PT")
    labels = as.character(dataset$cluster_label)
    labels[labels == "Lamp5 Sncg"] = "Sncg Lamp5"
    labels = reduce_to_prefix(labels, c(valid_subclasses(), exceptions))
    labels[labels == "L5 NP"] = "L5/6 NP"
    labels[labels == "L5 PT"] = "L5 ET"
    return(labels)
}

create_macosko_v3 = function() {
    input_dir = "10X_nuclei_v3_Broad"
    data_ = match_count_and_cluster_info(macosko_v3_counts(), zeng_clusters(input_dir, main_dir = joint_dir()))
    result = create_sce(data_$counts, data_$clusters, "macosko_10x_nuclei_v3")
    passed_initial_qc = !is.na(sce$is_qc)
    return(result[, passed_initial_qc])
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

if (!interactive()) {
    main()
}