import pandas as pd

SUFFIXES = ['genomic.gff', 'genomic.gtf', 'cds_from_genomic.fna', 'genomic.fna']
SAMPLES = ['J1']
ORGS = ['t4', 'spo1', 'ecoli', 'bsub']

rule all:
    input:
        expand("inputs/genomes_combined/all_{suffix}.gz", suffix = SUFFIXES),
        expand("outputs/salmon_quant/{sample}_quant/quant.sf", sample = SAMPLES),
        expand("outputs/bwa_align/{sample}.flagstat", sample = SAMPLES),
        expand("outputs/bwa_align/{sample}.bam.bai", sample = SAMPLES)


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
## Count transcripts
###########################################################################

rule salmon_index:
    input: "inputs/genomes_combined/all_cds_from_genomic.fna.gz"
    output: "outputs/salmon_indx/info.json"
    params: indx_dir = "outputs/salmon_indx"
    conda: "envs/salmon.yml"
    shell:'''
    salmon index -t {input} -i {params.indx_dir} -k 31
    '''

rule salmon_quant:
    input: 
        r1 = "inputs/raw/{sample}_1.fq.gz",
        r2 = "inputs/raw/{sample}_2.fq.gz",
        indx = "outputs/salmon_indx/info.json"
    output: "outputs/salmon_quant/{sample}_quant/quant.sf"
    params: 
        indx_dir = "outputs/salmon_indx",
        out_dir = lambda wildcards: "outputs/salmon_quant/" + wildcards.sample + "_quant"
    conda: "envs/salmon.yml"
    shell:'''
    salmon quant -i {params.indx_dir} -l A -1 {input.r1} -2 {input.r2} -o {params.out_dir} --validateMappings
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
   
rule flagstat:
    input: "outputs/bwa_align/{sample}.bam"
    output: "outputs/bwa_align/{sample}.flagstat"
    conda: "envs/bwa.yml"
    shell:'''
    samtools flagstat {input} > {output}
    '''
