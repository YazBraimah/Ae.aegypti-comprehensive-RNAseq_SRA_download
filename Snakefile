"""
Author: Y. Ahmed-Braimah
--- Snakemake workflow to process public Ae. aegypti RNAseq data
"""

import json
import pandas as pd
from os.path import join, basename, dirname
from os import getcwd
from subprocess import check_output
import os

##--------------------------------------------------------------------------------------##
## Functions
##--------------------------------------------------------------------------------------##

# To print process messages
def message(x):
  print()

# To remove suffix from a string
def rstrip(text, suffix):
    if not text.endswith(suffix):
        return text
    return text[:len(text)-len(suffix)]

##--------------------------------------------------------------------------------------##
## Global config parameters: 
##--------------------------------------------------------------------------------------##

configfile: 'config.yml'

# Load SRA runs file
RUN_INFO = pd.read_csv(config['SRA_INFO'])

# define SRA sample list
SAMPLES = RUN_INFO["Sample_Name"].tolist()

# split SRA sample list by library type (pe vs. se)
seSAMPLES = list(RUN_INFO[RUN_INFO["LibraryLayout"] == "SINGLE"]["Sample_Name"])
peSAMPLES = list(RUN_INFO[RUN_INFO["LibraryLayout"] == "PAIRED"]["Sample_Name"]) 

# define run and ID
RUNS = RUN_INFO["Run"].tolist()
LIBL = RUN_INFO["LibraryLayout"].tolist()

# Full path to a folder where final output files will be deposited.
OUT_DIR = config['OUT_DIR']
WORK_DIR = config['WORK_DIR']
HOME_DIR = config['HOME_DIR']  # the "launch_snakemake.sh" and "config.yml" files should be here

## set the usr and job environments for each job (specific for CBSU qsub jobs)
USER = os.environ.get('USER')
JOB_ID = os.environ.get('JOB_ID')


## Create the final output directory if it doesn't already exist
if not os.path.exists(OUT_DIR):
            os.makedirs(OUT_DIR)


##--------------------------------------------------------------------------------------##
## RULES
##--------------------------------------------------------------------------------------##


## Final expected output(s)
rule all: 
    input: 
    	expand(join(OUT_DIR, 'Reads', 'fastq', '{lib}', '{sample}', 'merged', '{sample}.R1.fastq.gz'), zip, sample=SAMPLES, run=RUNS, lib=LIBL),
        expand(join(OUT_DIR, 'Reads', 'fastq', 'PAIRED', '{sample}', 'merged', '{sample}.R1.fastq.gz'), zip, sample=peSAMPLES, run=RUNS, lib=LIBL),
        expand(join(OUT_DIR, 'Reads', 'fastq', '{lib}', '{sample}', '{run}', '{run}_1.fastq.gz'), zip, sample=peSAMPLES, run=RUNS, lib=LIBL)
        
##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to fetch raw SRA data
rule get_sra_data:
    output: 
        raw_sra = join(OUT_DIR, 'Reads', 'sra', '{lib}', '{sample}', '{run}', '{run}.sra')
    log:
        join(OUT_DIR, 'Reads', 'sra', '{lib}', '{sample}', '{run}', '{run}.log')
    benchmark:
        join(OUT_DIR, 'Reads', 'sra', '{lib}', '{sample}', '{run}', '{run}.benchmark.tsv')
    message: 
        """--- Downloading raw SRA record for "{wildcards.sample}" """
    params: 
        run_prefix=lambda wildcards: wildcards.run[:6], 
        sra_prefix=lambda wildcards: wildcards.run[:3]
    run:
        # shell('/programs/bin/labutils/mount_server cbsufsrv5 /data1')
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && wget -O ' + join(WORK_DIR, USER, JOB_ID, '{wildcards.run}.sra') + 
                ' ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/{params.sra_prefix}/{params.run_prefix}/{wildcards.run}/{wildcards.run}.sra')
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, '{wildcards.run}.sra') + ' ' + join(OUT_DIR, 'Reads', 'sra', '{wildcards.lib}', '{wildcards.sample}', '{wildcards.run}'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to convert SRA format raw file to fastq
rule get_fastq_files_from_sra_file:
    input: 
        raw_sra = rules.get_sra_data.output.raw_sra
    output: 
        raw_fastq = join(OUT_DIR, 'Reads', 'fastq', '{lib}', '{sample}', '{run}', '{run}_1.fastq.gz')
    log:
        join(OUT_DIR, 'Reads', 'fastq', '{lib}', '{sample}', '{run}', '{run}.log')
    benchmark:
        join(OUT_DIR, 'Reads', 'fastq', '{lib}', '{sample}', '{run}', '{run}.benchmark.tsv')
    message: 
        """--- running fastq-dump for "{wildcards.sample}" """
    run:
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.raw_sra} ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && fastq-dump --defline-seq \'@[$ac_]$sn/$ri\' --defline-qual \'+\' --split-files --gzip --outdir {wildcards.run}.output {wildcards.run}.sra')
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, '{wildcards.run}.output', '{wildcards.run}') + '*fastq.gz ' + join(OUT_DIR, 'Reads', 'fastq', '{wildcards.lib}', '{wildcards.sample}', '{wildcards.run}'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to merge multiple se fastq files from teh same sample
rule merge_se_fastq_files:
    # input:
    #     join(OUT_DIR, 'Reads', 'fastq', 'SINGLE', '{sample}', '{run}', '{run}_1.fastq.gz')
    output:
        join(OUT_DIR, 'Reads', 'fastq', 'SINGLE', '{sample}', 'merged', '{sample}.R1.fastq.gz')
    message: 
        """--- Merging SE fastq files for sample "{wildcards.sample}" """
    run:
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && cp ' + join(OUT_DIR, 'Reads', 'fastq', 'SINGLE', '{wildacrds.sample}', '*', '*_1.fastq.gz') + ' .' +        
                ' && zcat *_1.fastq.gz | gzip - > {wildcards.sample}.R1.fastq.gz')
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, '{wildcards.sample}.R1.fastq.gz') + ' ' + join(OUT_DIR, 'Reads', 'fastq', 'SINGLE', '{wildcards.sample}', 'merged', '{wildcards.sample}.R1.fastq.gz'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))


##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to merge multiple se fastq files from teh same sample
rule merge_left_pe_fastq_files:
    # input:
    #     join(OUT_DIR, 'Reads', 'fastq', 'PAIRED', '{sample}', '{run}', '{run}_1.fastq.gz')
    output:
        join(OUT_DIR, 'Reads', 'fastq', 'PAIRED', '{sample}', 'merged', '{sample}.R1.fastq.gz')
    message: 
        """--- Merging left fastq PE reads for sample "{wildcards.sample}" """
    run:
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && cp ' + join(OUT_DIR, 'Reads', 'fastq', 'PAIRED', '{wildacrds.sample}', '*', '*_1.fastq.gz') + ' .' +        
                ' && zcat *_1.fastq.gz | gzip - > {wildcards.sample}.R1.fastq.gz')
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, '{wildcards.sample}.R1.fastq.gz') + ' ' + join(OUT_DIR, 'Reads', 'fastq', 'PAIRED', '{wildcards.sample}', 'merged', '{wildcards.sample}.R1.fastq.gz'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))


##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to merge multiple se fastq files from teh same sample
rule merge_right_pe_fastq_files:
    # input:
    #     join(OUT_DIR, 'Reads', 'fastq', 'PAIRED', '{sample}', '{run}', '{run}_2.fastq.gz')
    output:
        join(OUT_DIR, 'Reads', 'fastq', 'PAIRED', '{sample}', 'merged', '{sample}.R2.fastq.gz')
    message: 
        """--- Merging right fastq PE reads for sample "{wildcards.sample}" """
    run:
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && cp ' + join(OUT_DIR, 'Reads', 'fastq', 'PAIRED', '{wildacrds.sample}', '*', '*_2.fastq.gz') + ' .' +        
                ' && zcat *_2.fastq.gz | gzip - > {wildcards.sample}.R2.fastq.gz')
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, '{wildcards.sample}.R2.fastq.gz') + ' ' + join(OUT_DIR, 'Reads', 'fastq', 'PAIRED', '{wildcards.sample}', 'merged', '{wildcards.sample}.R2.fastq.gz'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))
