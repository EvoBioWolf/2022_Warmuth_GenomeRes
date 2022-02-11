#!/bin/bash -l                                                                                                                                                                                                                                                             
#SBATCH -A snic2018-3-655 
#SBATCH -p core                                                                                                                                                                                                                                                             
#SBATCH -n 2
#SBATCH -t 240:00:00

# In this script, the de novo assembly of a female hooded crow is aligned
# to the assembled W chromosome of the New Caledonian crow (Corvus moneduloides)
# to identify syntenic (i.e., W-chromosome linked) contigs in the hooded crow.


module load bioinfo-tools lastz

lastz  Cmoneduloides_W_chr.fasta \
   /proj/uppstore2017156/b2010059_nobackup/S_Up_H60_assembly/flye_assembly_S_Up_H60_all_reads/Corvus_cornix__S_Up_H60__v1.0_flye.10kbFILTER.fasta \
   M=254 K=4500 L=3000 Y=15000 C=2 T=2 --matchcount=10000 --format=general:name1,start1,end1,length1,name2,start2,end2,strand2 > S_Up_H60_v1.0_1kbFILTER.vs.bCorMon1_W_chrom
