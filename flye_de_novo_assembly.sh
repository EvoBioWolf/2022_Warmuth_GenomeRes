#!/bin/bash -l
#SBATCH -A snic2018-3-655 
#SBATCH -p node
#SBATCH -C mem256GB
#SBATCH -n 20
#SBATCH -t 186:00:00
#SBATCH -J flye_genome_assembly 
#SBATCH --mail-type=ALL

# In this script, long-read sequencing data (Oxford Nanopore Technologies)
# is used to generate a de-novo assembly of the genome of a female hooded
# crow (Corvus cornix). 

flye --pacbio-raw /proj/uppstore2017156/b2010059/INBOX/Matthias_22_OCT_2019/blum.galaxy.lafuga.genzentrum.lmu.de/data/sequencers/by_group/Wolf_J/20190416_DNA_WolfJ-SUpH60_kraehe/fastq/20190416_DNA_WolfJ-SUpH60_kraehe-guppy-3.0.3-hac.porechop.fastq \
   /proj/uppstore2017156/b2010059/INBOX/Matthias_22_OCT_2019/blum.galaxy.lafuga.genzentrum.lmu.de/data/sequencers/by_group/Wolf_J/20190605_DNA_WolfJ_S_UP_60/fastq/20190605_DNA_WolfJ_S_UP_60-guppy-3.1.5-prom-hac.porechop.fastq \
   --out-dir flye_assembly_S_Up_H60_all_reads \
   --genome-size 1200m \
   --threads 20 \
   --resume
