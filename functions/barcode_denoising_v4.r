library(dplyr)
library(tidyverse)
argv=args=commandArgs(T)

script=argv[1]
kraken_report=argv[2]
sckmer=argv[3]
barcode_hit=argv[4]
barcode_tax=argv[5]
source(script)

# report = read.delim('output/kraken2/P1/SCAF2963_3_Live.kraken.report.txt', header = F)
report = read.delim(argv[2], header = F)
report$V8 = trimws(report$V8)
report[report$V8 %in% c('Homo sapiens', 'Bacteria', 'Fungi', 'Viruses'), ]

# sckmer data
kmer_data = read.table(argv[3], header = T,sep=' ')
head(kmer_data)

c = kmer_data %>% 
  subset(kmer > 1) %>%
  group_by(taxid) %>%
  mutate(nn = n()) %>%
  subset(nn > 3) %>% 
  group_by(taxid) %>%
  summarize(r =  cor.test(kmer, uniq, method = 'spearman')$estimate,
            p =   cor.test(kmer, uniq, method = 'spearman')$p.value,
            .groups='keep') %>%
  mutate(p = p.adjust(p))

print(c)
c$name = report$V8[match(c$taxid, report$V7)] # add taxa names 
#c=c[which(c$p < 0.05),]

k1 <- read_kraken_reports(c(argv[2]))
write.table(c, file=argv[4], sep="\t")
write.table(k1, file=argv[5], sep="\t")
