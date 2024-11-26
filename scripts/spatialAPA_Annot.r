args<-commandArgs(TRUE)
library(Seurat)
library(SCAPE)
library(reticulate)
library(Matrix)

####  指定样本id
res_path = args[1]
GSE = args[2]
GSM = args[3]
gtf_file = args[4]
# /public/home/jiangzh/projects/sc_apa/GRCh38.v44.gtf.gz
ensg_symbol = args[5]
# 

###### 数据集和样本号
dataset = GSE
samp = GSM

exp_file = file.path(res_path, GSE, GSM, "pasite.csv.gz")
names(exp_file) = GSM
collapse_pa = file.path(res_path, GSE, GSM, "collapse_pa.tsv.gz")

pa_mtx <- loadData(
  fileList = exp_file,
  collapsePa = collapse_pa,
  matrix = TRUE,
  cores = 8
)
colnames(pa_mtx) <- sapply(strsplit(colnames(pa_mtx), split = '[.]'), function(x) x[2])

# Load pa matrix into Seurat object
# Only these pA sites whcih expressed in more than 50 cell were kept.
binary_filter <- Matrix::rowSums(+pa_mtx)
pa_mtx <- pa_mtx[binary_filter > 50, ]
### annotation of pA
# gtf_file <- "/public/home/jiangzh/projects/sc_apa/GRCh38.v44.gtf.gz"
# It will consume a lot of time if it is the first time to annotate.
annot_info <- AnnotationSite(rownames(pa_mtx), gtf_file, 'Hg38',cores = 10)
gtf_v44_file = read.csv(ensg_symbol,sep="\t",row.names=1)
# gtf_v44_file = read.csv("/public/home/jiangzh/resource/ensg_symbol_type.txt",sep="\t",row.names=1)
annot_info$annot.symbol = gtf_v44_file[annot_info$annot.gene_id,"Symbol"]
#rownames(annot_info) = paste0(annot_info$annot.symbol,"|",rownames(annot_info))
annot_info$apa_id = paste0(annot_info$annot.symbol,"-",rownames(annot_info))

pa_mtx_fi = pa_mtx[rownames(annot_info),]
rownames(pa_mtx_fi) = annot_info$apa_id
rownames(annot_info) = annot_info$apa_id
annot_info$pa_site = annot_info$apa_id
# write.csv(annot_info,file.path(res_path, GSE, GSM, "APA_annotations.csv"))   ###################################  write
dir.create(file.path(res_path,"spatial_h5ad/apa/rds",dataset,samp), recursive = TRUE, showWarnings = FALSE)
write.csv(annot_info,file.path(res_path,"spatial_h5ad/apa/rds",dataset,samp,"APA_annotations.csv"))   ###################################  write

### 保存稀疏矩阵
saveRDS(pa_mtx_fi, file = file.path(res_path,GSE,GSM,"pa_mtx_fi.rds"))

