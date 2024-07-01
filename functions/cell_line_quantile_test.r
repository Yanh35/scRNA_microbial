argv=commandArgs(T)

library(dplyr)
library(tidyverse)
library(readr)


#cell.lines = read.delim('/hdd/yanhuan/micro/software/SAHMI/cell.lines.txt', header = T,sep=' ') %>% tibble() 
#cell.lines=cell.lines %>% 
#    filter (reads>0 & uniq >0)
#qtile = 0.90
#q_df = cell.lines %>%
#group_by(name,rank)%>%
#summarize(CLrpmm = 10^quantile(log10(rpmm),qtile, na.rm = T),
#.groups = 'keep')


q_df=read.delim('/hdd/yanhuan/micro/software/SAHMI/cell_line_quantile.txt', header = T,sep='\t') %>% tibble()
q_df=q_df[,c('name','rank','threshold_0.90_exclude0')]
colnames(q_df)=c('name','rank','CLrpmm')  
kr=read.csv(argv[2],header=T,sep='\t') 
sprintf("kraken rpmm: %d",nrow(kr))
bh=read.csv(argv[1],header=T,sep='\t') 
sprintf("Barcode level signal denoising_get_taxid: %d", nrow(bh)) 
merge_bn_qdf=left_join(bh,q_df, by = c('name' )) %>% left_join(kr %>% select(name,rpmm),by = c('name')) 
merge_bn_qdf$CLrpmm[is.na(merge_bn_qdf$CLrpmm)]= 0 
hits_df <- merge_bn_qdf %>% 
filter(rpmm > CLrpmm)
sprintf("cell line quantile test(exclude zero counts,0.90): %d",nrow(hits_df))
write.table(hits_df,argv[4],row.names=F,col.names=T,quote=F,sep='\t')
write.table(hits_df$taxid,argv[5],row.names=F,col.names=F,quote=F,sep='\t')
