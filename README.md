# 1 Kraken2 Database
Since [SAHMI](https://github.com/sjdlabgroup/SAHMI) is optimized to run with Kraken2Uniq, please refer to the Kraken2 manual for more details on the installation and usage of [Kraken2/KrakenUniq](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown). This analysis follows the steps to construct the Kraken2 database, primarily including human, bacterial, and viral.
## Human
We used the [CHM13 human genome](https://github.com/marbl/CHM13) as the reference, and the file [chm13v2.0.fa.gz](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz) was downloaded.
## Bacteria
`./kraken2/kraken2-build --download-library bacteria --db ./DATABASE`
## Viral
We obtained 15,037 RefSeq genome IDs from the [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Genome&SourceDB_s=RefSeq) and ultimately downloaded 18,524 sequences.
##
`kraken2/kraken2-build --download-taxonomy  --db ./DATABASE`


The downloaded human and viral sequences need to have their FASTA headers changed to the string kraken:taxid|XXX, with XXX replaced by the desired taxon ID. 
Next, use `kraken2-build --build --db ./DATABASE` to construct the database.



# 2 Identifying contaminants and false positives (cell line quantile test)
[SAHMI](https://github.com/sjdlabgroup/SAHMI) is designed to identify contaminant taxa and spurious false positive assignments by providing a negative control resource. We downloaded the [raw cell lines microbiome data](https://www.dropbox.com/s/r6xvw1589lqyqts/cell.lines.txt?dl=0) it provides and calculated the 90th percentile data. The 90th percentile values are calculated only from samples in which a taxon was reported with >0 reads and >0 unique k-mers. For more detailed descriptions, refer to the [SAHMI](https://github.com/sjdlabgroup/SAHMI).
```
library(dplyr)
library(readr)
cell.lines = read.delim('cell_lines.txt', header = T,sep=' ') %>% tibble()
cell.lines=cell.lines %>%
    filter (reads>0 & uniq >0)
qtile = 0.90
q_df = cell.lines %>%
group_by(name,rank)%>%
summarize(CLrpmm = 10^quantile(log10(rpmm),qtile, na.rm = T),
.groups = 'keep')
```
