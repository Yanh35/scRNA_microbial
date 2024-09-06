argv=commandArgs(T)

library(dplyr)
library(tidyverse)
library(readr)




q_df=read.delim('cell_line_quantile.txt', header = T,sep='\t') %>% tibble()
q_df=q_df[,c('name','rank','threshold_0.90_exclude0')]
colnames(q_df)=c('name','rank','CLrpmm')  
kr=read.csv(argv[2],header=T,sep='\t') 
kr=kr %>% subset(rank %in% c('G','S'))
sprintf("kraken rpmm: %d",nrow(kr))
bh=read.csv(argv[1],header=T,sep='\t') 
sprintf("Barcode level signal denoising_get_taxid: %d", nrow(bh)) 
#merge_bn_qdf=left_join(bh,q_df, by = c('name' )) %>% left_join(kr %>% select(name,rpmm),by = c('name')) 
#merge_bn_qdf$CLrpmm[is.na(merge_bn_qdf$CLrpmm)]= 0 
#hits_df <- merge_bn_qdf # %>% 
#filter(rpmm > CLrpmm)

merge_bn_qdf=left_join(bh,q_df, by = c('name' )) 
merge_bn_qdf_kr=left_join(kr %>% select(name,taxid,rpmm),merge_bn_qdf,by = c('taxid'))

merge_bn_qdf_kr$CLrpmm[is.na(merge_bn_qdf_kr$CLrpmm)]= 0
hits_df <- merge_bn_qdf_kr

sample_level=read.csv(argv[6],header=T,sep='\t')
#hits_df=hits_df[ which(hits_df$taxid %in% sample_level$taxid), ]
s1=left_join(hits_df,sample_level, by = c('taxid' ))
s1$name.x<-NULL
s1$rank.x<-NULL
s1$name=s1$name.y
s1$rank=s1$rank.y
s1$name.y<-NULL
s1$rank.y<-NULL


kr_report=read.csv(argv[7],header=F,sep='\t')
kr_1=unique(data.frame(kr_report$V6,kr_report$V7,kr_report$V8))
colnames(kr_1)=c('rank','taxid','name')
s1$rank=ifelse(s1$taxid %in% kr_1$taxid, kr_1[match(s1$taxid,kr_1$taxid),]$rank,NA)
s1$name=ifelse(s1$taxid %in% kr_1$taxid, trimws(kr_1[match(s1$taxid,kr_1$taxid),]$name),NA)
s1=s1 %>% subset(taxid != 9606 & taxid !=9605)
s1$mreads=ifelse(s1$taxid %in% kr$taxid, kr[match(s1$taxid,kr$taxid),]$reads,NA)

sprintf("cell line quantile test(exclude zero counts,0.90): %d",nrow(hits_df))
write.table(s1,argv[4],row.names=F,col.names=T,quote=F,sep='\t')
write.table(s1$taxid,argv[5],row.names=F,col.names=F,quote=F,sep='\t')
