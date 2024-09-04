library(dplyr)
library(tidyverse)
argv=commandArgs(T)
sam=unlist(strsplit(argv[1], ","))


all=data.frame()
for ( i in sam){
  path <- paste0('result/step5_v5/',i,'.kraken.report.rpmm.tsv')
  step5=read.csv(path,header=T,sep='\t')
  step5$study=rep('',dim(step5)[1])
  step5$sample=rep(i,dim(step5)[1])
  if (dim(all)[1] != 0) { all=rbind(all,step5) }
  if (dim(all)[1] == 0) { all=step5 }    
}
###############################
kr=all
kr = kr %>%
  group_by(taxid) %>%
  mutate(nn = n()) %>%
  subset(nn > 2) #%>%   select(-nn)
dim(all)
dim(kr)

c2 = kr %>%
  subset(rank %in% c('G', 'S')) %>% 
  group_by(taxid) %>%
  summarize(r1 = cor(min,uniq,method='spearman'),
            r2 = cor(min,reads,method='spearman'),
            r3 = cor(reads,uniq,method='spearman'),
            p1 = cor.test(min,uniq,method='spearman')$p.value,
            p2 = cor.test(min,reads,method='spearman')$p.value,
            p3 = cor.test(reads,uniq,method='spearman')$p.value
  )

######################################
kr=all
kr1 = kr %>%
  group_by(taxid) %>%
  mutate(nn = n()) %>%
  subset(nn < 3)
  
c2_1 = kr1 %>%
    subset(rank %in% c('G', 'S')) %>% 
    group_by(taxid) %>%
    summarize(r1 = NA,
              r2 = NA,
              r3 = NA,
              p1 = NA,
              p2 = NA,
              p3 = NA
    )
c2=data.frame(rbind(c2,c2_1))
################################

kr=all
kr = kr %>%
  group_by(taxid) %>%
  mutate(nn = n()) 

c3 = c2 %>% 
  left_join(select(kr, rank, taxid,name,nn) %>% distinct()) %>% 
  subset(rank %in% c('G', 'S'))

path=paste0(argv[2],'/sample_denoising.txt')
write.table(c3,path,row.names=F,col.names=T,sep='\t',quote=F)

