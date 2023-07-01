# reditools-snake-pipe
snakemake pipeline to run reditools on illumina data


## steps
bwa index reference genome
bwa aln map reads
split sense from antisense reads
sort & index bams
run reditools
merge reditools tables
