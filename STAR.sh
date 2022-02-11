### indexing and read mapping using STAR mapper ###

# STAR indexing
STAR --runThreadN 20 --runMode genomeGenerate --genomeDir index --genomeFastaFiles myFasta.fa --sjdbGTFfile myGTF.gtf --sjdbOverhang 49 --genomeSAindexNbases 13

path=/path/to/fastq

# STAR mapping using 'random mode'

for mate1 in $path/*fq.gz
do
prefix=$(echo ${mate1} | cut -f11 -d'/' | cut -f1-4 -d'_')

STAR --runThreadN 8 --runMode alignReads --genomeDir index --readFilesCommand zcat --readFilesIn $mate1  --sjdbGTFfile myGTF.gtf  --sjdbOverhang 49 --outFileNamePrefix $prefix --outSAMtype BAM SortedByCoordinate \
--outFilterMultimapNmax 1000 --outSAMmultNmax 1 --outFilterMismatchNmax 3 --outMultimapperOrder Random --winAnchorMultimapNmax 1000 --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 350

done
