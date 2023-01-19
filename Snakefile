SAMPLES = ['J1']
SUFFIXES = ['genomic.gff', 'genomic.gtf', 'cds_from_genomic.fna', 'genomic.fna']
ORGS = ['t4', 'spo1', 'ecoli', 'bsub']

rule all:
    input:
        expand("inputs/genomes/all_{suffix}.gz", suffix = SUFFIXES),
        expand("outputs/salmon_quant/{sample}_quant/quant.sf", sample = SAMPLES),
        expand("outputs/bwa_align/{sample}.flagstat", sample = SAMPLES),
        expand("outputs/bwa_align/{sample}.bam.bai", sample = SAMPLES)


###########################################################################
## Download reference genomes and "transcriptomes" for mapping and counting
###########################################################################

rule download_phage_T4:
    output: "inputs/genomes/t4_{suffix}.gz"
    shell:'''
    curl -JLo {output} https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/836/945/GCF_000836945.1_ViralProj14044/GCF_000836945.1_ViralProj14044_{wildcards.suffix}.gz
    '''

rule download_phage_spo1:
    output: "inputs/genomes/spo1_{suffix}.gz"
    shell:'''
    curl -JLo {output} https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/881/675/GCF_000881675.1_ViralProj32379/GCF_000881675.1_ViralProj32379_{wildcards.suffix}.gz
    '''

rule download_ecoli:
    output: "inputs/genomes/ecoli_{suffix}.gz"
    shell:'''
    curl -JLo {output} https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/559/635/GCF_001559635.1_ASM155963v1/GCF_001559635.1_ASM155963v1_{wildcards.suffix}.gz
    '''

rule download_bsub:
    output: "inputs/genomes/bsub_{suffix}.gz"
    shell:'''
    curl -JLo {output} https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1/GCF_000009045.1_ASM904v1_{wildcards.suffix}.gz
    '''

rule combine_references:
    input: expand("inputs/genomes/{org}_{{suffix}}.gz", org = ORGS)
    output: "inputs/genomes/all_{suffix}.gz"
    shell:'''
    cat {input} > {output}
    '''

###########################################################################
## Count transcripts
###########################################################################

rule salmon_index:
    input: "inputs/genomes/all_cds_from_genomic.fna.gz"
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
    input: "inputs/genomes/all_genomic.fna.gz"
    output: "inputs/genomes/all_genomic.fna"
    shell:'''
    gunzip -c {input} > {output}
    '''

rule bwa_index:
    input: "inputs/genomes/all_genomic.fna"
    output: "inputs/genomes/all_genomic.fna.bwt"
    conda: "envs/bwa.yml"
    shell:'''
    bwa index {input}
    ''' 

rule map_nucleotide_reads_against_assembly:
    input: 
        ref = "inputs/genomes/all_genomic.fna",
        ref_bwt = "inputs/genomes/all_genomic.fna.bwt",
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
