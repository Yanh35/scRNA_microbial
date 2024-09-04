import os
import pandas as pd
from itertools import combinations
import re
from os.path import join
import subprocess

workdir: config['Output_dir']
config_path = config["config_path"]
Output_dir = config['Output_dir']
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
        expand('result/step1/{sample}.kraken.report.txt', sample=sample),
        expand('result/step1/{sample}.kraken.report.mpa.txt', sample=sample),
        expand('result/step2/{sample}_1.fa', sample=sample),
        expand('result/step2/{sample}_2.fa', sample=sample),
        expand('result/step3/{sample}.microbiome.output.txt', sample=sample),
        expand('result/step4/{sample}.sckmer.txt', sample=sample), 
        expand('result/step5_v4/{sample}.barcode.kmer.hits.tsv', sample=sample),
        expand('result/step5_v4/{sample}.kraken.report.rpmm.tsv', sample=sample),
        'result/sample_denoising_v5/sample_denoising.txt',
        expand('result/step6_sample_v5/{sample}.cell_line_quantile_hits.tsv', sample=sample),
        expand('result/step6_sample_v5/{sample}.cell_line_quantile_hits_taxa.tsv', sample=sample),
        expand('result/step7_sample_v5/{sample}.counts.txt', sample=sample),
	'result/decontamn_sample_v5/decontamn.txt',
        expand('result/cellranger/{sample}.done',sample=sample),
        expand('result/micro_cellranger/{sample}_sahmi_cellranger.txt',sample=sample),

rule run_kraken:
    input: 
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
        Rscript functions/run_kraken.r \
        --sample {params.sam} \
        --fq1 {input[0]} \
        --fq2 {input[1]} \
        --out_path {params.outdir} \
        --ncbi_blast_path {ncbi_blast} \
        --Kraken2Uniq_path {Kraken2Uniq_path} \
        --kraken_database_path {kraken_database_path} \
        --kreport2mpa_path functions/kreport2mpa.py \
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
        Rscript functions/extract_microbiome_reads.r \
        --sample_name {params.sam}_1 \
        --fq {input[2]} \
        --kraken_report {input[0]} \
        --mpa_report {input[1]} \
        --out_path {params.outdir}


        Rscript functions/extract_microbiome_reads.r \
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
        Rscript functions/extract_microbiome_output.r --sample_name {params.sam} \
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
        Rscript functions/sckmer.r --sample_name={params.sam} \
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
        'result/step5_v4/{sample}.barcode.kmer.hits.tsv',
        'result/step5_v4/{sample}.kraken.report.rpmm.tsv'
    params: 
        outdir='result/step5_v4/'
    threads: 4
    shell:
        '''
        Rscript functions/barcode_denoising_v4.r functions/read_kraken_reports.r {input[0]} {input[1]} {output[0]} {output[1]}
        '''

rule sample_denoising:
	input: expand('result/step5_v4/{sample}.kraken.report.rpmm.tsv',sample=sample) 
	output: 'result/sample_denoising_v5/sample_denoising.txt'
	params:
		outdir='result/sample_denoising_v5',
		sam=config['sample']
	threads: 4
	shell:
		'''
Rscript functions/sample_denoising_v5.r {params.sam} {params.outdir}
		'''



rule cell_line_quantile_test:
    input: 
        'result/step5_v4/{sample}.barcode.kmer.hits.tsv',
        'result/step5_v4/{sample}.kraken.report.rpmm.tsv',
	'result/sample_denoising_v5/sample_denoising.txt',
	'result/step1/{sample}.kraken.report.txt'
    output: 
        'result/step6_sample_v5/{sample}.cell_line_quantile_hits.tsv',
        'result/step6_sample_v5/{sample}.cell_line_quantile_hits_taxa.tsv'
    threads: 4
    shell:
        '''
        Rscript functions/cell_line_quantile_test_v5.r {input[0]} {input[1]} '' {output[0]} {output[1]} {input[2]} {input[3]}
        '''

rule taxa_counts:
    input: 
        'result/step2/{sample}_1.fa',
        'result/step2/{sample}_2.fa',
        'result/step6_sample_v5/{sample}.cell_line_quantile_hits_taxa.tsv',
        'result/step1/{sample}.kraken.report.txt',
        'result/step1/{sample}.kraken.report.mpa.txt'
    output: 
        'result/step7_sample_v5/{sample}.counts.txt'
    params: 
        sam='{sample}',
        outdir='result/step7_sample_v5/'
    threads: 4
    shell:
        '''
        Rscript functions/taxa_counts.r \
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

rule decontamn:
    input: expand('result/step6_sample_v5/{sample}.cell_line_quantile_hits.tsv',sample=sample)
    output: 'result/decontamn_sample_v5/decontamn.txt'
    params:
        outdir='result/decotamn_sample_v5',
        sam=config['sample']
    threads: 4
    shell:
        '''
Rscript functions/decontamn.r {params.sam} {output} 0.1
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
/hdd/cheng/projects/MPXV/software/cellranger-7.2.0/cellranger count \
--id={params.sam} \
--transcriptome=/hdd/data/reference/Ensembl-GRCh38-2020A/refdata-gex-GRCh38-2020-A/ \
--fastqs={fq_path}/ \
--sample={params.sam} \
--output-dir={params.outdir}
                '''


rule micro_cellranger:
    input: 'result/decontamn_sample_v5/decontamn.txt',
        'result/cellranger/{sample}.done',
        'result/step7_sample_v5/{sample}.counts.txt',
        'result/step6_sample_v5/{sample}.cell_line_quantile_hits.tsv',
        'result/cellranger/{sample}/outs/filtered_feature_bc_matrix'
    output: 'result/micro_cellranger/{sample}_sahmi_cellranger.txt'
    threads: 2
    shell:
        '''
Rscript functions/micro_cellranger.r  {input[0]} {input[3]} {input[4]} {input[2]} {output}
        '''

