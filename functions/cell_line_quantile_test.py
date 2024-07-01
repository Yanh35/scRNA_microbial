import pandas as pd
import sys
kmer_hit=sys.argv[1]
kmer_rpmn=sys.argv[2]
cell=sys.argv[3]
cell_hit=sys.argv[4]
cell_hit_taxa=sys.argv[5]
# kmer_df = pd.read_csv("output/SAHMI/P1/SCAF2963_3_Live.barcode.kmer.hits.tsv", sep="\t")
kmer_df = pd.read_csv(sys.argv[1], sep="\t")
# report_df = pd.read_csv("output/SAHMI/P1/SCAF2963_3_Live.kraken.report.rpmm.tsv", sep="\t")
report_df = pd.read_csv(sys.argv[2], sep="\t")
# cell_lines_df = pd.read_excel("../../SAHMI/Table S4.xlsx")
cell_lines_df = pd.read_excel(sys.argv[3])
merged_df = kmer_df.merge(report_df.drop(columns=["study", "sample"]), on=["taxid", "name"])
merged_df = merged_df.loc[merged_df.p.notna()]
cell_lines_df = cell_lines_df[["name", "rank", "taxid", 0.90]]
merged_df = merged_df.merge(cell_lines_df, on=["name", "rank", "taxid"],how='outer')
merged_df=merged_df.fillna(value=0)
#merged_df.to_csv(sys.argv[6],sep='\t')
hits_df = merged_df.loc[merged_df.rpmm > merged_df[0.90]]
hits_df.to_csv(sys.argv[4], sep="\t")
hits_df["taxid"].to_csv(sys.argv[5], sep="\t", header=False, index=False)
