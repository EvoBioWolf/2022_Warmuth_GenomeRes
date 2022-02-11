### run NGM (v. 0.5.4) ###



for mate1 in *.fq
do

prefix=$(echo ${mate1} | cut -f22-24 -d'_')
ngm -q $mate1 --ref myGenome.fa --very-sensitive -i 0.88 -b -t 4 --output ./$prefix.ngm.bam

done
