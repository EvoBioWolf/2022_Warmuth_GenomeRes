#!/bin/bash -l
#SBATCH -A snic2018-3-655 
#SBATCH -p core
#SBATCH -n 1 
#SBATCH -t 48:00:00
#SBATCH -J map 

# In this script, long reads of a male hooded crow are aligned to the female
# de novo assembly and the average depth per contig is determined. Since the
# W chromosome is absent in male individuals, contigs with much lower than
# average read depth can be considered as W-linked.


REF=/proj/uppstore2017156/b2010059_nobackup/S_Up_H60_assembly/flye_assembly_S_Up_H60_all_reads/Corvus_cornix__S_Up_H60__v1.0_flye.10kbFILTER.fasta

ID=$(echo $IND | awk -F "__" '{print $18}')
proj=$(echo $IND | awk -F "__" '{print $21}')
label=$(echo ${proj}.${ID})
module load bioinfo-tools samtools


samtools depth S_Up_H37.vs.S_Up_H60_flye_v1.0.sorted.bam > S_Up_H37.vs.S_Up_H60_flye_v1.0.sorted.bam_coverage

while IFS= read -r line; do
   echo $line
	samtools depth S_Up_H37.vs.S_Up_H60_flye_v1.0.sorted.bam  -r "${line}" \
	| awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' |\
	paste <(echo $line) - >> coverage_depth_per_contig 
done<contig_list
