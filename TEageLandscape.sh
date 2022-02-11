### create repeat landscape from repeatmasker output (.align and masked genome)
#source: https://github.com/4ureliek/Parsing-RepeatMasker-Outputs

# whole genome
perl parseRM.pl -i myGenome.fa.align -f myGenome.fa.masked -n -p -l 50,1 -v

# autosomes only
perl parseRM.pl -i autosomes.fa.align -f autosomes.fa.masked -n -p -l 50,1 -v

# Z chromosome only
perl parseRM.pl -i chrZ.fa.align -f chrZ.fa.masked -n -p -l 50,1 -v

# W chromosome only
perl parseRM.pl -i Wscaffs.fa.align -f Wscaffs.fa.masked -n -p -l 50,1 -v




 