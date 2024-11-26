
args<-commandArgs(TRUE)
library(Seurat)
library(SCAPE)
library(magrittr)
library(rhdf5)
library(patchwork)
library(ggplot2)
library(dplyr)
library(reshape2)
library(tibble)
library(sceasy)
library(reticulate)
library(Matrix)
#library(jsonlite)
library(SingleCellExperiment)
library(ComplexHeatmap)
#library(Scillus)
library(gridExtra)
library(clusterProfiler)
library(org.Hs.eg.db)
library(RColorBrewer)

## usage: Rscript scape_prepecess.r GSE197023 GSM5907094

####  指定样本id and IN_PATH
res_path = args[1]
GSE = args[2]
GSM = args[3]

###### 数据集和样本号
dataset = GSE
samp = GSM
###### 数据集和样本号

spatial_dir <- file.path(res_path, GSE, GSM)
gene_obj <- Load10X_Spatial(spatial_dir)

counts <- GetAssayData(gene_obj, slot = "counts")
gene_totals <- rowSums(counts)
non_zero_genes <- names(gene_totals[gene_totals > 0])
gene_obj <- subset(gene_obj, features = non_zero_genes)

gene_obj <- NormalizeData(gene_obj)
gene_obj <- FindVariableFeatures(gene_obj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(gene_obj)
gene_obj <- ScaleData(gene_obj, features = all.genes)
gene_obj <- RunPCA(gene_obj, features = VariableFeatures(object = gene_obj))
gene_obj <- RunUMAP(gene_obj, dims = 1:30, reduction.name = "umap_Spatial", reduction.key = "SpatialUMAP_") ### umap
gene_obj <- FindNeighbors(gene_obj, dims = 1:30)
gene_obj <- FindClusters(gene_obj, resolution = 0.5, graph.name = 'Spatial_snn', cluster.name="Clusters_Spatial")  ### cluster 信息
all_markers_rna <- FindAllMarkers(gene_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_features_rna = all_markers_rna %>% group_by(cluster) %>% filter(p_val < 0.05)
dir.create(file.path(res_path,"tool","diff_information",dataset,samp), recursive = TRUE, showWarnings = FALSE)
write.csv(top_features_rna,file.path(res_path,"tool","diff_information",dataset,samp,"markers_top_Spatial.csv"),row.names=FALSE)   ###################################  write


########  cluster for apa
annot_info = read.csv(file.path(res_path,"spatial_h5ad/apa/rds",dataset,samp,"APA_annotations.csv"),row.names=1)
pa_mtx_fi = readRDS(file.path(res_path,GSE,GSM,"pa_mtx_fi.rds"))
gene_obj[['apa']] <- CreateAssayObject(pa_mtx_fi[, colnames(gene_obj)])
DefaultAssay(gene_obj) <- "apa"

gene_obj <- FindVariableFeatures(gene_obj, selection.method = "vst", nfeatures = 2000)
all.apas <- rownames(gene_obj)
gene_obj <- ScaleData(gene_obj, features = all.apas)
gene_obj <- RunPCA(gene_obj,verbose = FALSE, assay = "apa")

pca_embeddings <- Embeddings(gene_obj[["pca"]])
if(dim(pca_embeddings)[2]<30){
    gene_obj <- RunUMAP(gene_obj, dims = 1:dim(pca_embeddings)[2], reduction.name = "umap_APA", reduction.key = "APAUMAP_")
    gene_obj <- FindNeighbors(gene_obj, dims = 1:dim(pca_embeddings)[2])
}else {
    gene_obj <- RunUMAP(gene_obj, dims = 1:30, reduction.name = "umap_APA", reduction.key = "APAUMAP_")
    gene_obj <- FindNeighbors(gene_obj, dims = 1:30)
}

gene_obj <- FindClusters(gene_obj, resolution = 0.5, graph.name = 'apa_snn', cluster.name="Clusters_apa")
all_markers_apa <- FindAllMarkers(gene_obj, only.pos = TRUE)
top_features_apa = all_markers_apa %>% group_by(cluster) %>% filter(p_val < 0.05)
write.csv(top_features_apa,file.path(res_path,"tool","diff_information",dataset,samp,"markers_top_apa.csv"),row.names=FALSE)   ###################################  write


##### 转json对象 top APA events  #####
pa_mtx_summary = pa_mtx_fi[, colnames(gene_obj)]
non_zero_counts <- rowSums(pa_mtx_summary != 0)
non_zero_counts_sorted <- sort(non_zero_counts, decreasing = TRUE)
non_zero_counts_list <- as.list(non_zero_counts_sorted)
# 将每个值从数组转换为单一值
non_zero_counts_list <- lapply(non_zero_counts_list, function(x) x[1])
# # 转换为JSON对象
# #json_output <- toJSON(non_zero_counts_list, pretty = TRUE)
# json_output <- toJSON(non_zero_counts_list, pretty = TRUE, auto_unbox = TRUE)
# # 写入JSON文件
# dir.create(file.path(res_path,"spatial_h5ad/apa/overview",dataset,samp), recursive = TRUE, showWarnings = FALSE)
# write(json_output, file = file.path(res_path,"spatial_h5ad/apa/overview", dataset,samp,"APA_counts.json"))     ###################################  write
# ##### 转json对象 top APA events  #####

## calculate and classify psi
gene_obj <-
  psi(gene_obj,
      annot = annot_info,
      chunk = 4000,
      cores = 4)

###### figure ==  overview - APA_Category  ######
pa_cate <- SCAPE::psiCate(gene_obj, annot_info)

cate_colors = ggsci::pal_nejm()(7)
names(cate_colors) = c(
  "Multimodal",
  "NonExpr",
  "NonAPA",
  "L_shape",
  "J_shape",
  "OverDispersed",
  "UnderDispersed"
)

p1 <- ggplot(pa_cate, aes(mean_psi, SD_psi, color = PA_category)) +
        geom_point(size = 3, alpha = 0.6) +
        scale_color_manual(values = cate_colors) +
        labs(y = 'SD of PAratio', x = 'Mean of PAratio') +
        theme(
            plot.title = element_text(hjust = 0.5,
                                    face = 'bold',
                                    size = 16),
            legend.position = 'top',
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            axis.text.x = element_text(size = 15),
            axis.title.x = element_text(size = 15),
            axis.title.y = element_text(size = 15),
            axis.text.y  = element_text(size = 15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black")
        )

ggsave(file.path(res_path, "spatial_h5ad/apa/overview", dataset,samp, "overview_7_apa_catego.pdf"), plot = p1,width=10,height=6)     ###################################  write
###### figure ==  overview - APA_Category  ######


#######  cluster for apapsi
DefaultAssay(gene_obj) <- "apapsi"
gene_obj <- FindVariableFeatures(gene_obj, selection.method = "vst", nfeatures = 2000)
all.apapsi <- rownames(gene_obj)
gene_obj <- ScaleData(gene_obj, features = all.apapsi)
gene_obj <- RunPCA(gene_obj,verbose = FALSE, assay = "apapsi")

pca_embeddings <- Embeddings(gene_obj[["pca"]])
if(dim(pca_embeddings)[2]<30){
  gene_obj <- RunUMAP(gene_obj, dims = 1:dim(pca_embeddings)[2], reduction.name = "umap_PSI", reduction.key = "PSIUMAP_")
  gene_obj <- FindNeighbors(gene_obj, dims = 1:dim(pca_embeddings)[2])
}else {
  gene_obj <- RunUMAP(gene_obj, dims = 1:30, reduction.name = "umap_PSI", reduction.key = "PSIUMAP_")
  gene_obj <- FindNeighbors(gene_obj, dims = 1:30)
}

gene_obj <- FindClusters(gene_obj, resolution = 0.5, graph.name = 'apapsi_snn', cluster.name="Clusters_apapsi")
all_markers_apapsi <- FindAllMarkers(gene_obj, only.pos = TRUE)
top_features_apapsi = all_markers_apapsi %>% group_by(cluster) %>% filter(p_val < 0.05)
write.csv(top_features_apapsi,file.path(res_path,"tool","diff_information",dataset,samp,"markers_top_apapsi.csv"),row.names=FALSE)   ###################################  write


### 添加graph信息
gene_obj[["apa_snn_res.0.5"]] <- gene_obj$Clusters_apa
gene_obj[["apapsi_snn_res.0.5"]] <- gene_obj$Clusters_apapsi

##### top features
Idents(gene_obj) <- "Clusters_Spatial"

#######  cluster for apalen
pa_len_mtx <- read.csv(file.path(res_path, GSE, GSM, "gene_apa_length_mtx.csv"), check.names=FALSE)

# 检查第一列是否有缺失值，并去除包含缺失值的行
pa_len_mtx <- pa_len_mtx[!is.na(pa_len_mtx[, 1]), ]

if (any(duplicated(pa_len_mtx[, 1]))) {
  warning("First column contains duplicated values. Generating unique row names.")
  pa_len_mtx[, 1] <- make.unique(as.character(pa_len_mtx[, 1]))
}

# 将第一列设置为行名，并删除该列
rownames(pa_len_mtx) <- pa_len_mtx[, 1]
pa_len_mtx <- pa_len_mtx[, -1]

# 检查 NA 值并处理
# 检查并处理 Inf 和 NA 值
for (col in colnames(pa_len_mtx)) {
  pa_len_mtx[is.infinite(pa_len_mtx[, col]), col] <- NA
  pa_len_mtx[is.na(pa_len_mtx[, col]), col] <- 0
}

# common_spot = intersect(colnames(gene_obj),colnames(pa_len_mtx))
# 重新排列 pa_len_mtx 的列，并添加缺失的列
missing_spots <- setdiff(colnames(gene_obj), colnames(pa_len_mtx))
if (length(missing_spots) > 0) {
  # 创建一个全零的矩阵用于缺失的列
  zero_matrix <- matrix(0, nrow = nrow(pa_len_mtx), ncol = length(missing_spots))
  colnames(zero_matrix) <- missing_spots
  rownames(zero_matrix) <- rownames(pa_len_mtx)
  
  # 将零矩阵合并到 pa_len_mtx 中
  pa_len_mtx <- cbind(pa_len_mtx, zero_matrix)
}

non_zero_ratio <- rowMeans(pa_len_mtx > 0)
# 设置阈值，比如保留至少在90%的样本中有表达的基因
threshold <- 0.1  # 至少10%的样本中有表达的基因

# 过滤掉非零值比例低于阈值的基因
filtered_pa_len_mtx <- pa_len_mtx[non_zero_ratio >= threshold, ]
if(dim(filtered_pa_len_mtx)[1]!=0){
  gene_obj[["percent.mt"]] <- PercentageFeatureSet(gene_obj, pattern = "^MT-")
  gene_obj[['apalen']] <- CreateAssayObject(filtered_pa_len_mtx[, colnames(gene_obj)])
  DefaultAssay(gene_obj) <- "apalen"

  gene_obj <- FindVariableFeatures(gene_obj, selection.method = "vst", nfeatures = 2000)
  all.apalen <- rownames(gene_obj)
  gene_obj <- ScaleData(gene_obj, features = all.apalen)
  gene_obj <- RunPCA(gene_obj,verbose = FALSE, assay = "apalen")

  pca_embeddings <- Embeddings(gene_obj[["pca"]])
  if(dim(pca_embeddings)[2]<30){
    gene_obj <- RunUMAP(gene_obj, dims = 1:dim(pca_embeddings)[2], reduction.name = "umap_WUL", reduction.key = "WULUMAP_")
    gene_obj <- FindNeighbors(gene_obj, dims = 1:dim(pca_embeddings)[2])
  }else {
    gene_obj <- RunUMAP(gene_obj, dims = 1:30, reduction.name = "umap_WUL", reduction.key = "WULUMAP_")
    gene_obj <- FindNeighbors(gene_obj, dims = 1:30)
  }

  gene_obj <- FindClusters(gene_obj, resolution = 0.5, graph.name = 'apalen_snn', cluster.name="Clusters_apalen")
  all_markers_apalen <- FindAllMarkers(gene_obj, only.pos = TRUE)
  top_features_apalen = all_markers_apalen %>% group_by(cluster) %>% filter(p_val < 0.05)
  write.csv(top_features_apalen,file.path(res_path,"tool","diff_information",dataset,samp,"markers_top_apalen.csv"),row.names=FALSE)   ###################################  write

  gene_obj[["apalen_snn_res.0.5"]] <- gene_obj$Clusters_apalen
  select_assay = "apalen"
  gene_obj[[select_assay]]$scale.data <- NULL
  gene_obj[[select_assay]]$data <- NULL
  gene_obj[[select_assay]] <- as(gene_obj[[select_assay]], "Assay")
  sceasy::convertFormat(gene_obj, from="seurat", to="anndata", assay = select_assay,
                      outFile=file.path(res_path, "spatial_h5ad/apa/rds",dataset,samp, paste0("anndata_",select_assay,".h5ad"))) ###################################  write (checked)
}


gene_obj[["percent.mt"]] <- PercentageFeatureSet(gene_obj, pattern = "^MT-")
assays = c("Spatial","apa","apapsi")
for (select_assay in assays){
    DefaultAssay(gene_obj) <- select_assay

    gene_obj[[select_assay]]$scale.data <- NULL
    gene_obj[[select_assay]]$data <- NULL
    gene_obj[[select_assay]] <- as(gene_obj[[select_assay]], "Assay")
    sceasy::convertFormat(gene_obj, from="seurat", to="anndata", assay = select_assay,
                       outFile=file.path(res_path, "spatial_h5ad/apa/rds",dataset,samp, paste0("anndata_",select_assay,".h5ad"))) ###################################  write (checked)
}


DefaultAssay(gene_obj) <- "Spatial"
saveRDS(gene_obj,file.path(res_path,"spatial_h5ad/apa/rds",dataset,samp,"spatial_data.rds"))   ###################################  write (checked)

RDS_file = file.path(res_path,"spatial_h5ad/apa/rds",dataset,samp, "spatial_data.rds")  ################################################
gene_obj <- readRDS(RDS_file)
if ("apalen" %in% names(gene_obj@assays)) {
  filtered_pa_len_mtx <- matrix(0, nrow = 2, ncol = 2)
} else {
  filtered_pa_len_mtx <- matrix(0, nrow = 0, ncol = 1)
}

############################################################ overview_1 QC ############################################################
DefaultAssay(gene_obj) <- "Spatial"
plot_a <- SpatialFeaturePlot(gene_obj, features = "nCount_Spatial") + 
  labs(title="Gene") + 
  theme(legend.position = "right",plot.title = element_text(hjust = 0.5))

DefaultAssay(gene_obj) <- "apa"
plot_b <- SpatialFeaturePlot(gene_obj, features = "nCount_apa") + 
  labs(title="APA") + 
  theme(legend.position = "right",plot.title = element_text(hjust = 0.5))

DefaultAssay(gene_obj) <- "apapsi"
plot_c <- SpatialFeaturePlot(gene_obj, features = "nCount_apapsi") + 
  labs(title="PSI") + 
  theme(legend.position = "right",plot.title = element_text(hjust = 0.5))

DefaultAssay(gene_obj) <- "Spatial"
plot_e <- SpatialFeaturePlot(gene_obj, features = "nFeature_Spatial") + 
  # labs(title="Gene") + 
  theme(legend.position = "right",plot.title = element_text(hjust = 0.5))

DefaultAssay(gene_obj) <- "apa"
plot_f <- SpatialFeaturePlot(gene_obj, features = "nFeature_apa") + 
  # labs(title="APA") + 
  theme(legend.position = "right",plot.title = element_text(hjust = 0.5))

DefaultAssay(gene_obj) <- "apapsi"
plot_g <- SpatialFeaturePlot(gene_obj, features = "nFeature_apapsi") + 
  # labs(title="PSI") + 
  theme(legend.position = "right",plot.title = element_text(hjust = 0.5))

if(dim(filtered_pa_len_mtx)[1]!=0){
  DefaultAssay(gene_obj) <- "apalen"
  plot_d <- SpatialFeaturePlot(gene_obj, features = "nCount_apalen") + 
    labs(title="WUL") + 
    theme(legend.position = "right",plot.title = element_text(hjust = 0.5))
  plot_h <- SpatialFeaturePlot(gene_obj, features = "nFeature_apalen") + 
    # labs(title="WUL") + 
    theme(legend.position = "right",plot.title = element_text(hjust = 0.5))
  plot_all <- plot_a + plot_b + plot_c + plot_d + plot_e + plot_f + plot_g + plot_h + plot_layout(nrow=2, ncol = 4) + 
    theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm"))
  ggsave(file.path(res_path, "spatial_h5ad/apa/overview", dataset, samp, "overview_1_QC.pdf"), plot = plot_all, width = 12, height = 6) ###################################  write
}else{
  plot_all <- plot_a + plot_b + plot_c + plot_e + plot_f + plot_g + plot_layout(nrow=2, ncol = 3) + 
    theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm"))
  ggsave(file.path(res_path, "spatial_h5ad/apa/overview", dataset, samp, "overview_1_QC.pdf"), plot = plot_all, width = 10, height = 6)
}

############################################################ overview_1 QC ############################################################


############################################################ overview_2 cluster ############################################################
DefaultAssay(gene_obj) <- "Spatial"
Idents(gene_obj) = "Clusters_Spatial"
plot_1 <- SpatialDimPlot(gene_obj, label = TRUE, label.size = 3) +
  theme(plot.margin = unit(c(0.5, 0.05, 0.5, 0.05), "cm")) +
  labs(title="Gene") +
  theme(plot.title = element_text(hjust = 0.5)) ### 标题居中

DefaultAssay(gene_obj) <- "apa"
Idents(gene_obj) = "Clusters_apa"
plot_2 <- SpatialDimPlot(gene_obj, label = TRUE, label.size = 3) +
  theme(plot.margin = unit(c(0.5, 0.05, 0.5, 0.05), "cm")) +
  labs(title="APA") +
  theme(plot.title = element_text(hjust = 0.5)) ### 标题居中

DefaultAssay(gene_obj) <- "apapsi"
Idents(gene_obj) = "Clusters_apapsi"
plot_3 <- SpatialDimPlot(gene_obj, label = TRUE, label.size = 3) +
  theme(plot.margin = unit(c(0.5, 0.05, 0.5, 0.05), "cm")) +
  labs(title="PSI") +
  theme(plot.title = element_text(hjust = 0.5)) ### 标题居中

# ## "umap_Spatial", "umap_APA", "umap_PSI"

DefaultAssay(gene_obj) <- "Spatial"
Idents(gene_obj) = "Clusters_Spatial"
plot_5 <- DimPlot(gene_obj, reduction = "umap_Spatial") + 
  theme(legend.position = "right") +
  theme(plot.margin = unit(c(1, 0.05, 1, 0.05), "cm")) +
  labs(x = "Gene_UMAP_1", y = "Gene_UMAP_2")  ##修改坐标轴名字

DefaultAssay(gene_obj) <- "apa"
Idents(gene_obj) = "Clusters_apa"
plot_6 <- DimPlot(gene_obj, reduction = "umap_APA") + 
  theme(legend.position = "right") +
  theme(plot.margin = unit(c(1, 0.05, 1, 0.05), "cm")) +
  labs(x = "APA_UMAP_1", y = "APA_UMAP_2")  ##修改坐标轴名字

DefaultAssay(gene_obj) <- "apapsi"
Idents(gene_obj) = "Clusters_apapsi"
plot_7 <- DimPlot(gene_obj, reduction = "umap_PSI") + 
  theme(legend.position = "right") +
  theme(plot.margin = unit(c(1, 0.05, 1, 0.05), "cm")) +
  labs(x = "PSI_UMAP_1", y = "PSI_UMAP_2")  ##修改坐标轴名字

if(dim(filtered_pa_len_mtx)[1]!=0){
  DefaultAssay(gene_obj) <- "apalen"
  Idents(gene_obj) = "Clusters_apalen"
  plot_4 <- SpatialDimPlot(gene_obj, label = TRUE, label.size = 3) +
    theme(plot.margin = unit(c(0.5, 0.05, 0.5, 0.05), "cm")) +
    labs(title="WUL") +
    theme(plot.title = element_text(hjust = 0.5)) ### 标题居中

  plot_8 <- DimPlot(gene_obj, reduction = "umap_WUL") + 
    theme(legend.position = "right") +
    theme(plot.margin = unit(c(1, 0.05, 1, 0.05), "cm")) +
    labs(x = "WUL_UMAP_1", y = "WUL_UMAP_2")  ##修改坐标轴名字
  
  plot_all <- plot_1 + plot_2 + plot_3 + plot_4 + plot_5 + plot_6 + plot_7 + plot_8 + plot_layout(nrow=2, ncol = 4)
  #theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm"))

  # ggsave(file.path(outdir, "overview_2_cluster.pdf"), plot = plot_all, width=10, height=6)
  ggsave(file.path(res_path, "spatial_h5ad/apa/overview", dataset, samp, "overview_2_cluster.pdf"), plot = plot_all, width = 12, height = 6) ###################################  write
}else{
  plot_all <- plot_1 + plot_2 + plot_3 + plot_5 + plot_6 + plot_7 + plot_layout(nrow=2, ncol = 3)
  ggsave(file.path(res_path, "spatial_h5ad/apa/overview", dataset, samp, "overview_2_cluster.pdf"), plot = plot_all, width = 10, height = 6) ###################################  write
}
############################################################ overview_2 cluster ############################################################


############################################################ overview_3 heatmap ############################################################
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                        '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                        '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                        '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                        '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                        '#968175'
)

#### 读取top100差异基因列表
top_features_path = file.path(res_path,"tool","diff_information",dataset,samp)

heatmap_all <- function(gene_obj, top_features_path, select_assay_fi){
  DefaultAssay(gene_obj) <- select_assay_fi
  select_cluster_matrix = "Clusters_Spatial"
  Idents(gene_obj) <- select_cluster_matrix #选择不同的cluster
  gene_obj <- NormalizeData(gene_obj)
  gene_obj <- FindVariableFeatures(gene_obj, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(gene_obj)
  gene_obj <- ScaleData(gene_obj, features = all.genes)
  top_features <- read.csv(file.path(top_features_path, paste0("markers_top_",select_assay_fi,".csv")))
  
  #限定cluster展示数目
  cluster_ids <- unique(Idents(gene_obj))
  # if(length(cluster_ids)<=5){slice_num=10} else{slice_num=8}
  top_genes_per_cluster <- top_features %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 10) %>%
    #filter(!is.infinite(avg_log2FC)) %>%  # 排除 avg_log2FC 为 Inf 的行
    ungroup() %>%
    slice_head(n = 50)

  # Extract expression data
  expr_data <- GetAssayData(gene_obj, slot = "scale.data")[top_genes_per_cluster$gene, ]
  
  # Order the columns based on cluster information
  barcode_order <- order(Idents(gene_obj))
  expr_data <- expr_data[, barcode_order]
  
  # 使用my36colors的前10个颜色
  cluster_colors <- setNames(my36colors[1:length(cluster_ids)], cluster_ids)
  
  # Create annotation for clusters
  cluster_annotation <- HeatmapAnnotation(
    Clusters = as.character(Idents(gene_obj)[barcode_order]),
    col = list(Clusters = cluster_colors)
  )
  
  # 指定标题
  if(select_assay_fi=="Spatial"){
    column_title = "Heatmap of gene expression"
  }else if(select_assay_fi=="apa"){
    column_title = "Heatmap of APA expression"
  }else if(select_assay_fi=="apapsi"){
    column_title = "Heatmap of PSI value"
  }else if(select_assay_fi=="apalen"){
    column_title = "Heatmap of Weighted 3' UTR length"
  }
  
  # Create heatmap
  p_tmp <- Heatmap(
    expr_data, 
    name = "Expression",
    top_annotation = cluster_annotation,
    column_title = column_title,
    column_title_gp = gpar(fontsize = 16),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = FALSE
  )
  
  return(p_tmp)
}

p_gene <- heatmap_all(gene_obj,top_features_path,"Spatial")
p_apa <- heatmap_all(gene_obj,top_features_path,"apa")
p_psi <- heatmap_all(gene_obj,top_features_path,"apapsi")
if(dim(filtered_pa_len_mtx)[1]!=0){
  p_wul <- heatmap_all(gene_obj,top_features_path,"apalen")
  pdf(file = file.path(res_path, "spatial_h5ad/apa/overview", dataset, samp, "overview_3_heatmap.pdf"), width = 18, height = 13) ###################################  write
  draw(p_gene, padding = unit(c(2, 2, 2, 2), "cm"), newpage = TRUE) # Adjust padding as needed
  draw(p_apa, padding = unit(c(2, 2, 2, 2), "cm"), newpage = TRUE) # Adjust padding as needed
  draw(p_psi, padding = unit(c(2, 2, 2, 2), "cm"), newpage = TRUE) # Adjust padding as needed
  draw(p_wul, padding = unit(c(2, 2, 2, 2), "cm"), newpage = TRUE) # Adjust padding as needed
  dev.off()
}else{
  pdf(file = file.path(res_path, "spatial_h5ad/apa/overview", dataset, samp, "overview_3_heatmap.pdf"), width = 18, height = 13) ###################################  write
  draw(p_gene, padding = unit(c(2, 2, 2, 2), "cm"), newpage = TRUE) # Adjust padding as needed
  draw(p_apa, padding = unit(c(2, 2, 2, 2), "cm"), newpage = TRUE) # Adjust padding as needed
  draw(p_psi, padding = unit(c(2, 2, 2, 2), "cm"), newpage = TRUE) # Adjust padding as needed
  dev.off()
}

############################################################ overview_3 heatmap ############################################################


############################################################ overview_4 dotplot ############################################################
dotplot_all <- function(gene_obj, top_features_path, select_assay_fi){
  #### 读取top100差异基因列表
  # top_features <- read.csv(file.path(top_features_path, paste0("markers_top100_",select_assay_fi,".csv")))
  top_features <- read.csv(file.path(top_features_path, paste0("markers_top_",select_assay_fi,".csv")))
  # 使用dplyr进行分组并提取每个cluster的前5个gene
  top_genes_per_cluster <- top_features %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 10) %>%
    ungroup()
  Idents(gene_obj) <- "Clusters_Spatial" 
  DefaultAssay(gene_obj) <- select_assay_fi
  if(select_assay_fi=="Spatial"){title = "Gene"}
  else if(select_assay_fi=="apa"){title = "APA"}
  else if(select_assay_fi=="apapsi"){title = "PSI"}
  else if(select_assay_fi=="apalen"){title = "WUL"}
  else {title = select_assay_fi}
  # 归一化数据
  gene_obj <- NormalizeData(gene_obj)
  # 对数转换
  gene_obj <- ScaleData(gene_obj)
  dot_plot_tmp <- DotPlot(gene_obj, features = unique(top_genes_per_cluster$gene), cols = "RdYlBu") + RotatedAxis() +
    labs(title=title,x="") + 
    theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.5), "cm")) +
    theme(axis.text.x = element_text(angle = 90)) +
    theme(plot.title = element_text(hjust = 0.5)) ### 标题居中
  return(dot_plot_tmp)
}

p_gene <- dotplot_all(gene_obj,top_features_path,"Spatial")
p_apa <- dotplot_all(gene_obj,top_features_path,"apa")
p_psi <- dotplot_all(gene_obj,top_features_path,"apapsi")
if(dim(filtered_pa_len_mtx)[1]!=0){
  p_wul <- dotplot_all(gene_obj,top_features_path,"apalen")
  plot_all <- p_gene + p_apa + p_psi  + p_wul + plot_layout(nrow=4, ncol = 1)
  ggsave(file.path(res_path, "spatial_h5ad/apa/overview", dataset, samp, "overview_4_dotplot.pdf"), plot = plot_all, width = 18, height = 20) ###################################  write
}else{
  plot_all <- p_gene + p_apa + p_psi + plot_layout(nrow=3, ncol = 1)
  ggsave(file.path(res_path, "spatial_h5ad/apa/overview", dataset, samp, "overview_4_dotplot.pdf"), plot = plot_all, width = 18, height = 15) ###################################  write
}

############################################################ overview_4 dotplot ############################################################


############################################################ overview_5 WUL_boxplot ############################################################
if(dim(filtered_pa_len_mtx)[1]!=0){
  RDS_file = file.path(res_path,"spatial_h5ad/apa/rds",dataset,samp, "spatial_data.rds")
  gene_obj <- readRDS(RDS_file)
  select_cluster_matrix = "Clusters_Spatial" 
  Idents(gene_obj) <- select_cluster_matrix #选择不同的cluster

  # 确保 'apalen' assay 被正确添加到 gene_obj 中
  DefaultAssay(gene_obj) <- "apalen"
  # 提取 apalen 数据
  apalen_data <- GetAssayData(gene_obj, slot = "data")
  # 转置数据矩阵，使行是 spot，列是 gene
  apalen_data_transposed <- t(apalen_data)
  # 将数据转换为长格式
  apalen_data_long <- as.data.frame(as.matrix(apalen_data_transposed))
  apalen_data_long$spot <- rownames(apalen_data_long)
  # 添加 cluster 信息
  apalen_data_long$cluster <- Idents(gene_obj)[rownames(apalen_data_transposed)]
  # 将数据转换为长格式以适合 ggplot2 绘图
  apalen_data_long <- melt(apalen_data_long, id.vars = c("spot", "cluster"))
  # 计算99百分位数
  threshold <- quantile(apalen_data_long$value, 0.99)
  # 过滤掉高于阈值的离群值
  apalen_data_filtered <- apalen_data_long[apalen_data_long$value <= threshold, ]
  # 将表达值进行 log2 处理
  apalen_data_filtered$log2_value <- log2(apalen_data_filtered$value + 1)
  # 绘制小提琴图
  p <- ggplot(apalen_data_filtered, aes(x = as.factor(cluster), y = log2_value, fill = as.factor(cluster))) +
    geom_violin(trim = FALSE, scale = "width") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    scale_fill_brewer(palette = "Set3") +
    labs(title = "Violin Plot of WUL by cluster",
        x = "Cluster",
        y = "Log2 WUL") +
    theme_minimal(base_size = 15) +
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_blank(),
          panel.grid.major = element_line(size = 0.5, linetype = 'dashed', colour = "gray"),
          panel.grid.minor = element_blank())
  ggsave(file.path(res_path, "spatial_h5ad/apa/overview", dataset, samp, "overview_5_APAlength.pdf"), plot = p, width = 10, height = 8)
}
############################################################ overview_5 WUL_boxplot ############################################################



############################################################  preprocess_1_cluster_genes ############################################################
RDS_file = file.path(res_path,"spatial_h5ad/apa/rds",dataset,samp, "spatial_data.rds")
gene_obj <- readRDS(RDS_file)
select_cluster_matrix = "Clusters_Spatial" 
Idents(gene_obj) <- select_cluster_matrix #选择不同的cluster

if(dim(filtered_pa_len_mtx)[1]!=0){
  assays = c("Spatial","apa","apalen")
}else{
  assays = c("Spatial","apa")
}

for (select_assay in assays) {
    DefaultAssay(gene_obj) <- select_assay
    # 提取表达矩阵
    expression_matrix <- GetAssayData(gene_obj, slot = "data")
    # 找到表达量不为0的基因
    non_zero_genes <- sort(rownames(expression_matrix)[rowSums(expression_matrix != 0) > 0])
    json_data <- toJSON(non_zero_genes, pretty = TRUE)
    if (select_assay=="Spatial"){
        output_dir <- file.path(res_path,"tool","gene_information",dataset)
        # output_dir <- file.path("/root/projects/spatialapadb/tool/gene_information", dataset)
    }else if(select_assay=="apa"){
        output_dir <- file.path(res_path,"tool","apa_information",dataset)
        # output_dir <- file.path("/root/projects/spatialapadb/tool/apa_information", dataset)
    }else if(select_assay=="apalen"){
        output_dir <- file.path(res_path,"tool","apalen_information",dataset)
    }
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    # 保存JSON文件
    write(json_data, file = file.path(output_dir, paste0(samp, ".json")))
}

### 提取cluster信息
clusters = sort(unique(gene_obj[[select_cluster_matrix]])[,1])
json_data_cluster <- toJSON(clusters, pretty = TRUE)
output_dir_cluster <- file.path(res_path,"tool","cluster_information", dataset)
dir.create(output_dir_cluster, recursive = TRUE, showWarnings = FALSE)
# 保存JSON文件
write(json_data_cluster, file = file.path(output_dir_cluster, paste0(samp, ".json")))
######


##### 转换 marker文件
# library(dplyr)
if(dim(filtered_pa_len_mtx)[1]!=0){
  assays = c("Spatial","apa","apapsi","apalen")
}else{
  assays = c("Spatial","apa","apapsi")
}

for (select_mtx in assays) {
    df <- read.csv(file.path(res_path,"tool","diff_information", dataset, samp, paste0("markers_top_", select_mtx, ".csv")))
    # Get unique clusters
    clusters <- unique(df$cluster)
    # Initialize an empty list to store JSON objects
    json_list <- list()
    for (each_cluster in clusters) {
        # Subset data for each cluster
        tmp_mtx <- subset(df, cluster == each_cluster)
        # Create a list for the current cluster
        cluster_list <- lapply(1:nrow(tmp_mtx), function(i) {
            list(
                Cluster = as.numeric(tmp_mtx$cluster[i]), ##
                Gene = tmp_mtx$gene[i],
                p_val = as.numeric(tmp_mtx$p_val[i]),
                p_val_adj = as.numeric(tmp_mtx$p_val_adj[i]),
                logfoldchange = as.numeric(tmp_mtx$avg_log2FC[i])
            )
        }) 
        json_list <- append(json_list, cluster_list)
    }
    json_data <- toJSON(json_list, pretty = TRUE, auto_unbox = TRUE)
    # Define output directory
    output_dir <- file.path(res_path,"tool","diff_information", dataset, samp)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    # Save JSON file
    if (select_mtx=="Spatial"){
        select_mtx_fi = "Gene"
    }else if(select_mtx=="apa"){
        select_mtx_fi = "APA"
    }else if(select_mtx=="apapsi"){
        select_mtx_fi = "PSI"
    }else if(select_mtx=="apalen"){
        select_mtx_fi = "WUL"
    }
    write(json_data, file = file.path(output_dir, paste0("markers_top_", select_mtx_fi, ".json")))
}

for(select_assay in assays){
  list_dat <- read.csv(file.path(res_path,"tool","diff_information", dataset, samp,paste0("markers_top_",select_assay,".csv")))
  # outpath = "/root/projects/spatialapadb/tool/enrich_information"
  outpath = file.path(res_path,"tool","enrich_information")
  clusters = unique(list_dat$cluster)
  if(select_assay=="Spatial" || select_assay=="apalen"){
    gene.df <- bitr(list_dat$gene, fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)
    Gene = list_dat$gene
  }else{
    genes <- sapply(strsplit(list_dat$gene, split = "-chr"), function(x) x[1])
    gene.df <- bitr(genes, fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)
    Gene = genes
  }
  newdat = data.frame(Cluster=list_dat$cluster, Gene = Gene, 
                      p_val=list_dat$p_val, p_val_adj=list_dat$p_val_adj, 
                      logfoldchange=list_dat$avg_log2FC)
  # 将 list_dat 和 gene.df 按照 SYMBOL 列进行合并，并将 ENTREZID 添加到 list_dat 中
  merged_dat <- merge(newdat, gene.df, by.x = "Gene", by.y = "SYMBOL", all.x = TRUE)
  merged_dat = na.omit(merged_dat)
  for (cluster in clusters){
    gene_list = subset(merged_dat,Cluster==cluster)$ENTREZID
    # 对基因列表进行GO富集分析
    go_enrich <- enrichGO(gene = gene_list, 
                          OrgDb = org.Hs.eg.db,  
                          keyType = "ENTREZID",  # 指定输入基因ID的类型
                          ont = "ALL",           # 选择GO分析的方向，这里选择"ALL"表示所有,可以是BP/MF/CC
                          pvalueCutoff = 0.05,   # 设定P值的阈值
                          pAdjustMethod = "BH",  # 多重假设检验校正方法，这里使用Benjamini & Hochberg的方法
                          qvalueCutoff = 0.1,    # 设定FDR的阈值
                          readable=TRUE)         # 将结果的entrezID转为symbol
    if (is.null(go_enrich) || nrow(go_enrich@result) == 0) {
      empty_result <- data.frame(ONTOLOGY = NA, ID = NA, Description = NA, GeneRatio = NA, BgRatio = NA, 
                                 pvalue = NA, p.adjust = NA, qvalue = NA, geneID = NA, Count = NA)
      go_enrich <- new("enrichResult", result = empty_result, pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                       organism = "UNKNOWN", ontology = "ALL", gene = gene_list, keytype = "ENTREZID")
    }
    # 对基因列表进行KEGG富集分析
    kegg_enrich <- enrichKEGG(gene = gene_list,
                              organism = 'hsa',     
                              pvalueCutoff = 0.2,   # 设定P值的阈值
                              pAdjustMethod = "BH",  # 多重假设检验校正方法，这里使用Benjamini & Hochberg的方法
                              qvalueCutoff = 0.5)    # 设定FDR的阈值
    if (is.null(kegg_enrich) || nrow(kegg_enrich@result) == 0) {
      empty_result <- data.frame(category = NA, subcategory = NA, ID = NA, Description = NA, 
                                 GeneRatio = NA, BgRatio = NA, pvalue = NA, p.adjust = NA, 
                                 qvalue = NA, geneID = NA, Count = NA)
      kegg_enrich <- new("enrichResult", result = empty_result, pvalueCutoff = 0.2, pAdjustMethod = "BH", 
                         organism = "hsa", ontology = "KEGG", gene = gene_list, keytype = "ENTREZID")
    }
    kegg_enrich <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    # 将enrichResult对象转换为数据框
    go_enrich_df <- as.data.frame(go_enrich)
    kegg_enrich_df <- as.data.frame(kegg_enrich)
    rownames(go_enrich_df) <- NULL
    rownames(kegg_enrich_df) <- NULL
    # # Remove "_row" column
    # go_enrich_df <- go_enrich_df[ , !names(go_enrich_df) %in% "_row"]
    # kegg_enrich_df <- kegg_enrich_df[ , !names(kegg_enrich_df) %in% "_row"]
    # 将数据框转换为JSON格式
    json_data_go <- toJSON(go_enrich_df, pretty = TRUE)
    json_data_kegg <- toJSON(kegg_enrich_df, pretty = TRUE)
    # 创建保存路径
    output_dir <- file.path(outpath, dataset, samp)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    # 保存JSON文件
    write(json_data_go, file = file.path(output_dir, paste0(cluster, "_enrich_go_", select_assay, ".json")))
    write(json_data_kegg, file = file.path(output_dir, paste0(cluster, "_enrich_kegg_", select_assay, ".json")))
    # 保存CSV文件
    write.csv(go_enrich_df, file.path(output_dir, paste0(cluster, "_enrich_go_", select_assay, ".csv")), row.names = FALSE)
    write.csv(kegg_enrich_df, file.path(output_dir, paste0(cluster, "_enrich_kegg_", select_assay, ".csv")), row.names = FALSE)
  }
}


for(select_assay in assays){
    list_dat <- read.csv(file.path(res_path,"tool","diff_information", dataset, samp,paste0("markers_top_",select_assay,".csv")))
    # outpath = "/root/projects/spatialapadb/tool/enrich_information"
    outpath = file.path(res_path,"tool","enrich_information")
    clusters = unique(list_dat$cluster)
    if(select_assay=="Spatial"){
        gene.df <- bitr(list_dat$gene, fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)
        Gene = list_dat$gene
    }else{
        genes <- sapply(strsplit(list_dat$gene, split = "-chr"), function(x) x[1])
        gene.df <- bitr(genes, fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)
        Gene = genes
    }
    newdat = data.frame(Cluster=list_dat$cluster, Gene = Gene, 
                        p_val=list_dat$p_val, p_val_adj=list_dat$p_val_adj, 
                        logfoldchange=list_dat$avg_log2FC)
    # 将 list_dat 和 gene.df 按照 SYMBOL 列进行合并，并将 ENTREZID 添加到 list_dat 中
    merged_dat <- merge(newdat, gene.df, by.x = "Gene", by.y = "SYMBOL", all.x = TRUE)
    merged_dat = na.omit(merged_dat)
    for (cluster in clusters){
        gene_list = subset(merged_dat,Cluster==cluster)$ENTREZID
        # 对基因列表进行GO富集分析
        go_enrich <- enrichGO(gene = gene_list, 
                            OrgDb = org.Hs.eg.db,  #
                            keyType = "ENTREZID",  # 指定输入基因ID的类型
                            ont = "ALL",           # 选择GO分析的方向，这里选择"ALL"表示所有,可以是BP/MF/CC
                            pvalueCutoff = 0.05,   # 设定P值的阈值
                            pAdjustMethod = "BH",  # 多重假设检验校正方法，这里使用Benjamini & Hochberg的方法
                            qvalueCutoff = 0.1,    # 设定FDR的阈值
                            readable=TRUE)         # 将结果的entrezID转为symbol
        # 对基因列表进行KEGG富集分析
        kegg_enrich <- enrichKEGG(gene = gene_list,
                            organism = 'hsa',      #
                            pvalueCutoff = 0.2,   # 设定P值的阈值
                            pAdjustMethod = "BH",  # 多重假设检验校正方法，这里使用Benjamini & Hochberg的方法
                            qvalueCutoff = 0.5)    # 设定FDR的阈值
        kegg_enrich <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
        # 将enrichResult对象转换为数据框
        go_enrich_df <- as.data.frame(go_enrich)
        kegg_enrich_df <- as.data.frame(kegg_enrich)
        # 将数据框转换为JSON格式
        json_data_go <- toJSON(go_enrich_df, pretty = TRUE)
        json_data_kegg <- toJSON(kegg_enrich_df, pretty = TRUE)
        # 创建保存路径
        output_dir <- file.path(outpath, dataset, samp)
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        # 保存JSON文件
        write(json_data_go, file = file.path(output_dir, paste0(cluster, "_enrich_go_", select_assay, ".json")))
        write(json_data_kegg, file = file.path(output_dir, paste0(cluster, "_enrich_kegg_", select_assay, ".json")))
        # 保存CSV文件
        write.csv(go_enrich_df, file.path(output_dir, paste0(cluster, "_enrich_go_", select_assay, ".csv")), row.names = FALSE)
        write.csv(kegg_enrich_df, file.path(output_dir, paste0(cluster, "_enrich_kegg_", select_assay, ".csv")), row.names = FALSE)
    }
}

############################################################  preprocess_1_cluster_genes ############################################################

cat("finish preprocess: ", GSM,"\n")

