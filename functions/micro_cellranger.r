library(dplyr)
library(tidyverse)
library(optparse)
library(readr)
library(readxl)
library(scales)
library(Seurat)
library(reshape2)
library(parallel)
argv=commandArgs(T)
decontam=read.csv(argv[1],header = T,sep = '\t') #result/decontamn_sample_v5/decontamn.txt


#path <- paste0('result/step6_sample_v5/', i, '.cell_line_quantile_hits.tsv')
cell <- read.csv(argv[2], header = TRUE, sep = '\t')
cell$decontamn <- ifelse(cell$taxid %in% decontam$taxid, 
                         decontam[match(cell$taxid, decontam$taxid),]$p, NA)

#path <- paste('result/cellranger/', i, "/outs/filtered_feature_bc_matrix", sep = '')
data <- Read10X(data.dir = argv[3])
seurat_object <- CreateSeuratObject(counts = data, project = 'prj', min.cells = 0, min.features = 0)
#seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
express <- seurat_object@assays$RNA@counts


#path <- paste('result/step7_sample_v5/', i, '.counts.txt', sep = '')
step7 <- read.delim(argv[4], header = TRUE)

test <- function(tax, cou) {
  co_210 <- cou %>% filter(taxid == tax)
  if (nrow(co_210) == 0) return(c(NA, NA, NA, NA))
  
  
  Count_mtx1 <- acast(co_210, barcode ~ taxid, value.var = "counts", fill = 0)
  barcode_210 <- rownames(Count_mtx1)
  t1 <- rowSums(Count_mtx1, na.rm = TRUE)
  result_210 <- data.frame(barcode_210, umi_210 = t1)
  

  c1 <- cou %>% filter(rank == 'k') %>% select(barcode, counts, taxid)
  Count_mtx1 <- acast(c1, barcode ~ taxid, value.var = "counts", fill = 0)
  barcode <- rownames(Count_mtx1)
  t1 <- rowSums(Count_mtx1, na.rm = TRUE)
  result <- data.frame(barcode, TotalUMIs = t1)
  

  result$umis <- ifelse(result$barcode %in% result_210$barcode_210, 
                        result_210[match(result$barcode, result_210$barcode_210), ]$umi_210, NA)
  result$include_tax <- ifelse(result$barcode %in% co_210$barcode, 'Y', 'N')
  result$include_cellranger <- ifelse(paste0(result$barcode, '-1') %in% colnames(express), 'Y', 'N')
  
 
  s1=result %>% subset( include_tax == 'Y' )
  total_tax_umi=sum(s1$umis)
  
  inc=s1 %>% subset(include_cellranger == 'Y')
  outc=s1 %>% subset(include_cellranger == 'N')
  inc_umi=sum(inc$umis)
  outc_umi=sum(outc$umis)
  
  inc=result %>% subset(include_cellranger == 'Y' & include_tax == 'Y')
  outc=result %>% subset(include_cellranger == 'N' & include_tax == 'Y')
  inc_cell=dim(inc)[1]
  outc_cell=dim(outc)[1]
  
  s1$include_cellranger=factor(s1$include_cellranger,level=c('Y','N'))
  if (sum(s1$include_cellranger == 'Y') < 2 || sum(s1$include_cellranger == 'N') < 2 ||
      var(s1$umis[s1$include_cellranger == 'Y']) == 0 ||
      var(s1$umis[s1$include_cellranger == 'N']) == 0) {
    ttest_t=NA
    ttest_p=NA
  } else {
    ttest<-with(s1, t.test( umis ~ include_cellranger ) ) ###cellranger barcode.include tax vs not include tax
    ttest_t=ttest$statistic[[1]]
    ttest_p=ttest$p.value[[1]]}
  
  return(c(inc_umi/inc_cell,outc_umi/outc_cell,ttest_t,ttest_p)) 
}


n_cores <- detectCores() - 1 
now()
res <- mclapply(cell$taxid, function(x) { test(x, step7) }, mc.cores = 40)
now()
result_df <- do.call(rbind, res)
col_cell=colnames(cell)
cell=cbind(cell,result_df)
colnames(cell)=c(col_cell,'tax_umi_include_cellranger','tax_umi_exclude_cellranger','tax_umi_cellranger_t','tax_umi_cellranger_p')


mg=read.csv('mgnify_modified.txt',header = T,sep='\t')
cell$mgnify <- ifelse(cell$name %in% mg$taxonomy, 
                      mg[match(cell$name, mg$taxonomy),]$tax_study, '')
write.table(cell,argv[5],row.names=F,col.names=T,sep='\t',quote=F)
