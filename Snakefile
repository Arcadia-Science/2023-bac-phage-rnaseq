import pandas as pd

SUFFIXES = ['genomic.gff', 'genomic.gtf', 'genomic.fna']
SAMPLES = ['J1']
ORGS = ['t4', 'spo1', 'ecoli', 'bsub']

rule all:
    input:
        expand("inputs/genomes_combined/all_{suffix}.gz", suffix = SUFFIXES),
        expand("outputs/bwa_align/{sample}.flagstat", sample = SAMPLES),
        expand("outputs/bwa_align/{sample}.bam.bai", sample = SAMPLES),
        expand("outputs/bwa_align/{sample}.depth", sample = SAMPLES),
        expand("outputs/bwa_align/{sample}.idxstats", sample = SAMPLES),
        expand("outputs/featurecounts/{sample}_featurecounts_counts.txt", sample = SAMPLES)

###########################################################################
## Download reference genomes and "transcriptomes" for mapping and counting
###########################################################################

rule download_references:
    input: "inputs/genomes.csv"
    output: "inputs/genomes_raw/{org}-{suffix}.gz"
    params: out_prefix = lambda wildcards: "inputs/genomes_raw/" + wildcards.org + "-"
    run:
        genomes = pd.read_csv(str(input[0]))
        row = genomes.loc[genomes['genome'] == wildcards.org]
        prefix = row['url_prefix'].values[0]
        url = prefix + wildcards.suffix + ".gz"
        shell("curl -JLo {params.out_prefix}{wildcards.suffix}.gz {url}")


rule combine_references:
    input: expand("inputs/genomes_raw/{org}-{{suffix}}.gz", org = ORGS)
    output: "inputs/genomes_combined/all_{suffix}.gz"
    shell:'''
    cat {input} > {output}
    '''

###########################################################################
## Map reads to genome to create BAM file
###########################################################################

rule gunzip_ref:
    input: "inputs/genomes_combined/all_genomic.fna.gz"
    output: "inputs/genomes_combined/all_genomic.fna"
    shell:'''
    gunzip -c {input} > {output}
    '''

rule bwa_index:
    input: "inputs/genomes_combined/all_genomic.fna"
    output: "inputs/genomes_combined/all_genomic.fna.bwt"
    conda: "envs/bwa.yml"
    shell:'''
    bwa index {input}
    ''' 

rule map_nucleotide_reads_against_assembly:
    input: 
        ref = "inputs/genomes_combined/all_genomic.fna",
        ref_bwt = "inputs/genomes_combined/all_genomic.fna.bwt",
        r1 = "inputs/raw/{sample}_1.fq.gz",
        r2 = "inputs/raw/{sample}_2.fq.gz",
    output: "outputs/bwa_align/{sample}.bam"
    conda: "envs/bwa.yml"
    shell:'''
    bwa mem {input.ref} {input.r1} {input.r2} | samtools sort -o {output} -
    '''

rule samtools_index:
    input: "outputs/bwa_align/{sample}.bam"
    output: "outputs/bwa_align/{sample}.bam.bai"
    conda: "envs/bwa.yml"
    shell:'''
    samtools index {input}
    '''
   
rule samtools_flagstat:
    input: "outputs/bwa_align/{sample}.bam"
    output: "outputs/bwa_align/{sample}.flagstat"
    conda: "envs/bwa.yml"
    shell:'''
    samtools flagstat {input} > {output}
    '''

rule samtools_depth:
    input: "outputs/bwa_align/{sample}.bam"
    output: "outputs/bwa_align/{sample}.depth"
    conda: "envs/bwa.yml"
    shell:'''
    samtools depth {input} > {output}
    '''

rule samtools_idxstats:
    input: "outputs/bwa_align/{sample}.bam"
    output: "outputs/bwa_align/{sample}.idxstats"
    conda: "envs/bwa.yml"
    shell:'''
    samtools idxstats {input} > {output}
    '''

rule featurecounts:
    input:
        bam = "outputs/bwa_align/{sample}.bam",
        gtf = "inputs/genomes_combined/all_genomic.gtf.gz"
    output: 
        counts="outputs/featurecounts/{sample}_featurecounts_counts.txt",
        stat="outputs/featurecounts/{sample}_featurecounts_stats.txt"
    conda: "envs/rsubread.yml"
    script: "scripts/snakemake_featurecounts.R"
