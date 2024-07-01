import os
import pandas as pd
from itertools import combinations
import re
from os.path import join
import subprocess

workdir: config['Output_dir']
config_path = config["config_path"]
Output_dir = config['Output_dir']
SAHMI = config['SAHMI']
fq_path = config['fq_path']
ncbi_blast = config['ncbi_blast']
Kraken2Uniq_path = config['Kraken2Uniq_path']
kraken_database_path = config['kraken_database_path']
sample = config['sample'].strip().split(',')
umi = config['umi']
cb = config['cb']
r1=config['R1']
r2=config['R2']

rule all:
    input: 
#        expand('result/step1/{sample}.kraken.report.txt', sample=sample),
#        expand('result/step1/{sample}.kraken.report.mpa.txt', sample=sample),
#        expand('result/step2/{sample}_1.fa', sample=sample),
#        expand('result/step2/{sample}_2.fa', sample=sample),
#        expand('result/step3/{sample}.microbiome.output.txt', sample=sample),
#        expand('result/step4/{sample}.sckmer.txt', sample=sample),
#        expand('result/step5/{sample}.barcode.kmer.hits.tsv', sample=sample),
#        expand('result/step5/{sample}.kraken.report.rpmm.tsv', sample=sample),
#        expand('result/step6/{sample}.cell_line_quantile_hits.tsv', sample=sample),
#        expand('result/step6/{sample}.cell_line_quantile_hits_taxa.tsv', sample=sample),
#       expand('result/step7/{sample}.counts.txt', sample=sample),
        expand('result/cellranger/{sample}.done',sample=sample)

rule run_kraken:
    input: 
#        join(fq_path, '{sample}_{r1}.fastq.gz'),
#        join(fq_path, '{sample}_{r2}.fastq.gz')
        fq1 = lambda wildcards: join(fq_path, f"{wildcards.sample}_{r1}.fastq.gz"),
        fq2 = lambda wildcards: join(fq_path, f"{wildcards.sample}_{r2}.fastq.gz")
    output: 
        'result/step1/{sample}.kraken.report.txt',
        'result/step1/{sample}.kraken.report.mpa.txt',
        'result/step1/{sample}.kraken.output.txt',
        'result/step1/{sample}_1.fq',
        'result/step1/{sample}_2.fq'
    params:
        sam='{sample}',
        outdir='result/step1/'
    threads: 4
    shell:
        '''
        mkdir -p {params.outdir}
        Rscript {SAHMI}/functions/run_kraken.r \
        --sample {params.sam} \
        --fq1 {input[0]} \
        --fq2 {input[1]} \
        --out_path {params.outdir} \
        --ncbi_blast_path {ncbi_blast} \
        --Kraken2Uniq_path {Kraken2Uniq_path} \
        --kraken_database_path {kraken_database_path} \
        --kreport2mpa_path {SAHMI}/functions/kreport2mpa.py \
        --paired T
        '''

rule extract:
    input: 
        'result/step1/{sample}.kraken.report.txt',
        'result/step1/{sample}.kraken.report.mpa.txt',
        'result/step1/{sample}_1.fq',
        'result/step1/{sample}_2.fq'
    output: 
        'result/step2/{sample}_1.fa',
        'result/step2/{sample}_2.fa'
    params:
        sam='{sample}',
        outdir='result/step2/'
    threads: 4
    shell:
        '''
        Rscript {SAHMI}/functions/extract_microbiome_reads.r \
        --sample_name {params.sam}_1 \
        --fq {input[2]} \
        --kraken_report {input[0]} \
        --mpa_report {input[1]} \
        --out_path {params.outdir}


        Rscript {SAHMI}/functions/extract_microbiome_reads.r \
        --sample_name {params.sam}_2 \
        --fq {input[3]} \
        --kraken_report {input[0]} \
        --mpa_report {input[1]} \
        --out_path {params.outdir}
        '''

rule extract_microbiome_output:
    input: 
        'result/step1/{sample}.kraken.report.txt',
        'result/step1/{sample}.kraken.report.mpa.txt',
        'result/step1/{sample}.kraken.output.txt'
    output:
        'result/step3/{sample}.microbiome.output.txt'
    params:
        sam='{sample}',
        outdir='result/step3/'
    threads: 4
    shell:
        '''
        Rscript {SAHMI}/functions/extract_microbiome_output.r --sample_name {params.sam} \
        --output_file {input[2]} \
        --kraken_report {input[0]} \
        --mpa_report {input[1]} \
        --out_path {params.outdir}
        '''

rule sckmer:
    input: 
        'result/step2/{sample}_1.fa',
        'result/step2/{sample}_2.fa',
        'result/step1/{sample}.kraken.report.txt',
        'result/step1/{sample}.kraken.report.mpa.txt',
        'result/step3/{sample}.microbiome.output.txt'
    output: 
        'result/step4/{sample}.sckmer.txt'
    params:
        sam='{sample}',
        outdir='result/step4/'
    threads: 4
    shell:
        '''
        Rscript {SAHMI}/functions/sckmer.r --sample_name={params.sam} \
        --fa1={input[0]} \
        --fa2={input[1]} \
        --microbiome_output_file={input[4]} \
        --kraken_report={input[2]} \
        --mpa_report={input[3]} \
        --out_path={params.outdir} \
        --cb_len={cb} \
        --umi_len={umi} \
        --ranks=c\\('G','S'\\)  \
        --host=9606 \
        '''

rule barcode_denoising:
    input: 
        'result/step1/{sample}.kraken.report.txt',
        'result/step4/{sample}.sckmer.txt'
    output: 
        'result/step5/{sample}.barcode.kmer.hits.tsv',
        'result/step5/{sample}.kraken.report.rpmm.tsv'
    params: 
        outdir='result/step5/'
    threads: 4
    shell:
        '''
        Rscript {SAHMI}/functions/barcode_denoising.r {SAHMI}/functions/read_kraken_reports.r {input[0]} {input[1]} {output[0]} {output[1]}
        '''

rule cell_line_quantile_test:
    input: 
        'result/step5/{sample}.barcode.kmer.hits.tsv',
        'result/step5/{sample}.kraken.report.rpmm.tsv'
    output: 
        'result/step6/{sample}.cell_line_quantile_hits.tsv',
        'result/step6/{sample}.cell_line_quantile_hits_taxa.tsv'
    threads: 4
    shell:
        '''
        Rscript {SAHMI}/functions/cell_line_quantile_test.r {input[0]} {input[1]} {SAHMI}/cell_line_quantile.txt {output[0]} {output[1]}
        '''

rule taxa_counts:
    input: 
        'result/step2/{sample}_1.fa',
        'result/step2/{sample}_2.fa',
        'result/step6/{sample}.cell_line_quantile_hits_taxa.tsv',
        'result/step1/{sample}.kraken.report.txt',
        'result/step1/{sample}.kraken.report.mpa.txt'
    output: 
        'result/step7/{sample}.counts.txt'
    params: 
        sam='{sample}',
        outdir='result/step7/'
    threads: 4
    shell:
        '''
        Rscript {SAHMI}/functions/taxa_counts.r \
        --sample_name {params.sam} \
        --fa1 {input[0]} \
        --fa2 {input[1]} \
        --taxa {input[2]} \
        --kraken_report {input[3]} \
        --mpa_report {input[4]} \
        --out_path {params.outdir} \
        --cb_len {cb} \
        --umi_len {umi}
        '''
rule cellranger:
        input: join(fq_path,'{sample}_S1_L001_R1_001.fastq.gz'),
            join(fq_path,'{sample}_S1_L001_R2_001.fastq.gz')
        output: touch('result/cellranger/{sample}.done')
        params: sam='{sample}',
                outdir='result/cellranger/{sample}'
        threads:  4
        shell:
                '''
/hdd/cheng/projects/MPXV/software/cellranger-7.2.0/cellranger count --id={params.sam} --transcriptome=/hdd/data/reference/Ensembl-GRCh38-2020A/refdata-gex-GRCh38-2020-A/ \
--fastqs={fq_path}/ \
--sample={params.sam} \
--output-dir={params.outdir}
                '''
