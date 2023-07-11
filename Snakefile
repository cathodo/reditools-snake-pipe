
##########################################
########### reditools pipe ###############
##########################################

### which R
SIDE = "R2"

import os
from glob import glob

configfile: "config.yaml"
master_threads = 48
GENFILEDIR=os.path.dirname(config["inputs"]["reference"])

glob_list = glob(os.path.join(config["inputs"]["fastq_dir"], "*"+SIDE+"*.fastq*"))
SAMPLE_DICT = dict()
for s in glob_list:
    sample_name = os.path.basename(s).split(".")[0]
    SAMPLE_DICT[sample_name] = s

rule all:
    input:
        expand("output/redifolders/{sample}", sample=SAMPLE_DICT.keys())

rule bwa_mem:
    log:
        o = "logs/bwa_mem/{sample}.stderr1",
        e = "logs/bwa_mem/{sample}.stderr2",
    input:
        fq = lambda wildcards: SAMPLE_DICT[wildcards.sample],
        ref = config["inputs"]["reference"],
    output:
        sam = "output/bwa_mem/{sample}.sam",
    threads: master_threads
    shell:
        "if [ ! -e {input.ref}.bwt ]; then bwa index {input.ref} 2> {log.o}; fi "
        "&& bwa mem -t {threads} {input.ref} {input.fq} > {output.sam} 2> {log.e} "


rule bwa_aln:
    log:
        o = "logs/bwa_aln/{sample}.stderr1",
        e = "logs/bwa_aln/{sample}.stderr2",
    input:
        fq = lambda wildcards: SAMPLE_DICT[wildcards.sample],
        ref = config["inputs"]["reference"]
    output:
        sai = "output/bwa_aln/{sample}.sai",
    threads: master_threads
    shell:
        "if [ ! -e {input.ref}.bwt ]; then bwa index {input.ref} 2> {log.o}; fi "
        "&& bwa aln -t {threads} {input.ref} {input.fq} > {output.sai} 2> {log.e} "

rule bwa_samse:
    log:
        o = "logs/bwa_samse/{sample}.stdout",
        e = "logs/bwa_samse/{sample}.stderr",
    input:
        sai = "output/bwa_aln/{sample}.sai",
        ref = config["inputs"]["reference"],
        fq = lambda wildcards: SAMPLE_DICT[wildcards.sample]
    output:
        sam = "output/bwa_samse/{sample}.sam",
    threads: master_threads
    shell:
        "bwa samse {input.ref} {input.sai} {input.fq} > {output.sam} 2> {log.e} "

rule samtools_view:
    log:
        o = "logs/samtools_view/{sample}.stdout",
        e = "logs/samtools_view/{sample}.stderr",
    input:
        sam = "output/bwa_mem/{sample}.sam",
        ref = config["inputs"]["reference"],
    output:
        bam = "output/samtools_view/{sample}.bam",
    threads: master_threads
    shell:
        "samtools view -@ {threads} -bS -T {input.ref} {input.sam} > {output.bam} 2>> {log.e}"

rule separate_sense:
    log:
        o = "logs/separate_sense/{sample}.stderr1",
        e = "logs/separate_sense/{sample}.stderr2",
    input: "output/samtools_view/{sample}.bam"
    output:
        plus = "output/separate_sense/{sample}.plus.bam",
        minus = "output/separate_sense/{sample}.minus.bam",
    threads: master_threads
    shell:
        "samtools view -@ {threads} -b -F 16 {input} > {output.plus} 2> {log.o} "
        "&& samtools view -@ {threads} -b -f 16 {input} > {output.minus} 2> {log.e} "


rule index_bam:
    log:
        o = "logs/index_bam/{sample}.stdout",
        e = "logs/index_bam/{sample}.stderr",
    input: "output/separate_sense/{sample}.plus.bam"
    output: "output/index_bam/{sample}.sorted.bam"
    threads: master_threads
    shell:
        "samtools sort -@ {threads} -o {output} {input} 1> {log.o} 2> {log.e} "
        "&& samtools index -@ {threads} {output} 1>> {log.o} 2>> {log.e}"

rule reditools:
    resources:
        time = "48:00:00",
        mem_gb = 32,
    singularity: "singularity/reditools.sif"
    log:
        o = "logs/reditools/{sample}.stdout",
        e = "logs/reditools/{sample}.stderr",
    input:
        bam = "output/index_bam/{sample}.sorted.bam",
        ref = config["inputs"]["reference"],
    output: "output/redifolders/{sample}"
    threads: master_threads
    shell:
        "'python REDItoolDenovo.py -i {input.bam} -f "
        "{input.ref} -o {output} -t {threads} "
        "1> {log.o} 2> {log.e}' "
