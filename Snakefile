#!/usr/bin/env snakemake

# Load configuration file
configfile: "./config.yaml"


# Reading a file, that contains SRA Run Ids
with open(config["file_with_genomes"], "r") as f:
    analysed_genomes = dict(line.strip().split("\t") for line in f)


rule all:
    input:
        expand("rrn_genome_variants/{fas_name}.vcf", fas_name=analysed_genomes.keys()),


# Using fasterq-dump to get fastq files
rule get_fastq:
    output:
        temp("raw_reads/{fas_name}_1.fastq"),
        temp("raw_reads/{fas_name}_2.fastq"),
    threads: config["get_fastq_threads"]
    params:
        sample=lambda wildcards: analysed_genomes[wildcards.fas_name],
    shell:
        """
        fasterq-dump --split-files --threads {threads} -O raw_reads {params.sample}
        mv raw_reads/{params.sample}_1.fastq {output[0]}
        mv raw_reads/{params.sample}_2.fastq {output[1]}
        """


# Using bwa for mapping fastq files to a reference and creating a bam file.
rule bwa_map_short:
    input:
        fasta="ref_short_genome/revc_who_f_2024.fa",
        fq1="raw_reads/{fas_name}_1.fastq",
        fq2="raw_reads/{fas_name}_2.fastq",
    output:
        bam=temp("bam_files/{fas_name}.bam"),
        bai=temp("bam_files/{fas_name}.bam.bai"),
    threads: config["bwa_map_short_threads"]
    params:
        rg=lambda wildcards: rf"@RG\tID:{analysed_genomes[wildcards.fas_name]}\tSM:{analysed_genomes[wildcards.fas_name]}",
    shell:
        "bwa mem -t {threads} -R '{params.rg}' {input.fasta} {input.fq1} {input.fq2} "
        "| samtools sort --threads {threads} > {output.bam} && samtools index {output.bam} --threads {threads} -o {output.bai}"


# Creating a vcf file with freebayes for rrn operon
rule call_variants_short:
    input:
        fa="ref_short_genome/revc_who_f_2024.fa",
        bam="bam_files/{fas_name}.bam",
        bai="bam_files/{fas_name}.bam.bai",
    output:
        temp("pre_normalized_vcf/{fas_name}.vcf"),
    threads: config["normalize_short_variants_threads"]
    shell:
        "freebayes -f {input.fa} -p 4 --region CP145052.1:1058983-1064385 {input.bam} > {output}"


# Normalization of a vcf file
rule normalize_short_variants:
    input:
        fa="ref_short_genome/revc_who_f_2024.fa",
        vcf="pre_normalized_vcf/{fas_name}.vcf",
    output:
        "rrn_genome_variants/{fas_name}.vcf",
    threads: 1
    shell:
        "bcftools norm -f {input.fa} --threads {threads} {input.vcf} -O v > {output}"
