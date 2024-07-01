library(dplyr)
library(tidyverse)
argv=commandArgs(T)
data=read.csv(argv[1],header = T,sep='\t',check.names = F)

for (i in c('k','p','c','o','f','g','s'))
{
data_s=data[which(data$rank %in% i),]
d1=data_s[,c("barcode",'taxid','counts')]
matrix_df <- d1 %>% 
  spread(key = barcode, value = counts)
matrix_df <- replace(matrix_df, is.na(matrix_df), 0)
rownames(matrix_df)=matrix_df$taxid
matrix_df$taxid=NULL
path=paste(argv[3],'/',argv[2],'_exprss_matrix_',i,'.csv',sep='')
write.table(data.frame(taxid=rownames(matrix_df),matrix_df),path,sep = '\t',row.names = F,col.names = T,quote = F)
}
