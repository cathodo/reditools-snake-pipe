
##########################################
########### reditools pipe ###############
##########################################

### which R
SIDE = "R2"

import os
from glob import glob 

configfile: "config.yaml"
master_threads = 48

glob_list = glob(os.path.join(config["inputs"]["fastq_dir"], "*"+SIDE+"*.fastq*"))
SAMPLE_DICT = dict()
for s in glob_list:
    sample_name = os.path.basename(s).split(".")[0]
    SAMPLE_DICT[sample_name] = s

rule all:
    input:
        expand("output/redifolders/{sample}", sample=SAMPLE_DICT.keys())

rule bwa_aln:
    log:
        o = "logs/bwa_aln/{sample}.stderr1",
        e = "logs/bwa_aln/{sample}.stderr2",
    input: 
        fq = lambda wildcards: SAMPLE_DICT[wildcards.sample],
        ref = config["inputs"]["reference"]
    output: "output/bwa_aln/{sample}.sai"
    threads: master_threads
    shell:
        "if [ ! -e {input.ref}.bwt ]; then bwa index {input.ref} 2> {log.o}; fi "
        "&& bwa aln -n 18 -t {threads} {input.ref} {input.fq} > {output} "
        "2> {log.e} "

rule bwa_samse:
    log:
        o = "logs/bwa_samse/{sample}.stdout",
        e = "logs/bwa_samse/{sample}.stderr",
    input: 
        fq = lambda wildcards: SAMPLE_DICT[wildcards.sample],
        ref = config["inputs"]["reference"],
	sai = "output/bwa_aln/{sample}.sai",
    output: "output/bwa_samse/{sample}.sam"
    threads: master_threads
    shell:
        "bwa samse {input.ref} {input.sai} {input.fq} > {output} 2> {log.e}"

rule sam_to_bam:
    log:
        o = "logs/sam_to_bam/{sample}.stdout",
        e = "logs/sam_to_bam/{sample}.stderr",
    input: 
        sam = "output/bwa_samse/{sample}.sam",
        ref = config["inputs"]["reference"],
    output: "output/sam_to_bam/{sample}.bam"
    threads: master_threads
    shell:
        "samtools view -@ {threads} -bS -T {input.ref} {input.sam} "
        "> {output} 2> {log.e} "

rule separate_sense:
    log:   
        o = "logs/separate_sense/{sample}.stderr1",
        e = "logs/separate_sense/{sample}.stderr2",
    input: "output/sam_to_bam/{sample}.bam"
    output: 
        plus = "output/separate_sense/{sample}.plus.bam",
        minus = "output/separate_sense/{sample}.minus.bam",
    threads: master_threads
    shell:
        "samtools view -@ {threads} -b -F 16 {input} > {output.plus} 2> {log.o} "
        "samtools view -A {threads} -b -f 16 {input} > {output.minus} 2> {log.e} "

rule index_bam:
    log:
        o = "logs/index_bam/{sample}.stdout",
        e = "logs/index_bam/{sample}.stderr",
    input: "output/separate_sense/{sample}.plus.bam"
    output: "output/index_bam/{sample}.sorted.bam"
    threads: master_threads
    shell:
        "samtools sort -A {threads} -o {output} {input} 1> {log.o} 2> {log.e} "
        "&& samtools index -A {threads} {output} 1>> {log.o} 2>> {log.e}"

rule reditools:
    singularity:
        "singularity_files/reditools.sif"
    log:
        o = "logs/reditools/{sample}.stdout",
        e = "logs/reditools/{sample}.stderr",
    input: 
        bam = "output/index_bam/{sample}.sorted.bam",
        ref = config["inputs"]["reference"],
    output: "output/redifolders/{sample}"
    threads: master_threads
    shell:
        "docker run -i --cpus {threads} --memory 120G reditools -c "
        "'python /REDItools/main/REDItoolDenovo.py -I {input.bam} -f "
        "{input.ref} -o {output} -t {threads}' 1> {log.o} 2> {log.e}"

