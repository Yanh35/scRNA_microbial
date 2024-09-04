library(dplyr)
library(tidyverse)
library(optparse)
library(readr)
library(readxl)
library(scales)
library(repr)
library(ggplot2)
library("decontam")
library("reshape2")
library("stats")
argv=commandArgs(T)
sam=unlist(strsplit(argv[1], ","))


args <- list(
  resout = argv[2],  # 输出结果的文件名
  threshold = as.numeric(argv[3]) ,  # 阈值
  verbose = TRUE# 是否打印额外输出
)



PRJNA=data.frame()
for (i in sam){
  path=paste('result/step6_sample_v5/',i,'.cell_line_quantile_hits.tsv',sep='')
  rpmm=read.csv(path,header = T,sep='\t')
  rpmm$run=rep(i,dim(rpmm)[1])
  rpmm$study=rep('project',dim(rpmm)[1])
  PRJNA=rbind(PRJNA,rpmm)}


d1=PRJNA[,c('rpmm','taxid','run')]
matrix_df <- d1 %>% 
  spread(key = run, value = rpmm)
matrix_df <- replace(matrix_df, is.na(matrix_df), 0)
rownames(matrix_df)=matrix_df$taxid
matrix_df$taxid<-NULL 

count_table <- matrix_df
ta=as.character(rownames(count_table))
count_table=data.frame(t(count_table))
colnames(count_table)=ta
count_table=count_table[sam,]
rows_table <- rownames(count_table)
count_matrix <- data.matrix(count_table)



micro=data.frame()
for (i in sam){
  path=paste('result/step1/',i,'.kraken.report.txt',sep='')
  data=read.csv(path,header = F,sep='\t')
  unclass=data[which(data$V7 == 0),][1,2]
  class=data[which(data$V7 == 1),][1,2]
  human=data[which(data$V7 == 9606),][1,3]
  root=data[which(data$V7 == 1),][1,3]
  cellular=data[which(data$V7 == 131567),][1,3]
  bacviru=sum(data[which(data$V7 %in% c(2,10239)),]$V2)
  reads=unclass+human+root+cellular+bacviru
  count=as.numeric(c(unclass,human,root,cellular,bacviru,reads))
  micro=rbind(micro,count)
}
colnames(micro)=c('unclassified','human','root','cellular','bacviru','reads')
micro$sample=sam
#micro$unclssified_percent=micro$unclassified/micro$reads
#micro$human_percent=micro$human/micro$reads
#micro$root_cellular_bacviru=rowSums(data.frame(micro$bacviru,micro$root,micro$cellular))
#micro$root_cellular_bacviru_percent=micro$root_cellular_bacviru/micro$reads

concentrations <- micro
concentrations=concentrations[,c('sample','bacviru')]
colnames(concentrations)=c('V1','V2')
concentrations_list <- concentrations[ (concentrations[, "V1"] %in% rows_table) , "V2"]
decontam_out <- isContaminant(count_matrix, normalize=FALSE, conc=concentrations_list, method="frequency", threshold=args$threshold)

write.table(data.frame(taxid=rownames(decontam_out),decontam_out), file=args$resout, sep="\t", quote=FALSE,row.names=F)

