# Quantifying abundances of transcripts using Kallisto
by Francesco Delongu
## Prepare transcripts database

### Go to the folder where the bins are
```
for i in *.fa; do j=`echo $i | awk -F"." '{print $2}'`; awk -v var=$j '/>/{print substr($1, 2, length($1))"\t"var}' $i >> bin_map.txt; done
awk '{a++; if(a>1) print $1}' annotations.txt | awk -F"_" '{a[$0"\t"$1"_"$2"_"$3"_"$4"_"$5"_"$6]=1}END{for(i in a)print i}' > ORF_contig.txt
awk -F"\t" 'NR==NRF{a[$1]=$2; next} {print $0"\t"a[$2]}' bin_map.txt ORF_contig.txt > ORF_contig_bin.txt

awk '{if($1~/^>/){print header"\n"seq; header=$1; seq=""} else seq=seq""$1}END{print header"\n"seq}' nostops.ffn | awk '{a++; if(a>2)print $0}' > nostops_2lines.ffn 
awk 'NR==FNR{a[">"$1]=1; next} {if(a[$1]==1){print $1; getline; print $1}}' ORF_contig_bin.txt nostops_2lines.ffn > genes.fna

module add kallisto
kallisto index -i genes.fna.kallisto_index --make-unique genes.fna

```
### Quantify transcripts
```
cd ~
mkdir kallisto
```
### Build sbatch file (e.g. quantifier.sh):

```
	#!/bin/bash

        #SBATCH -N 1
        #SBATCH -J kallisto
        #SBATCH -c 8
        #SBATCH -a 1-21

        path1='/mnt/project/Courses/BIN420-2019/Day3_MTandNetAnalysis/data/'
        path0='/mnt/users/student54/BIN420/Wednesday'
        path2='/mnt/users/student54/kallisto'

        module add kallisto

        i=$(( SLURM_ARRAY_TASK_ID - 1 ))

        i1=$(( i * 2 ))
        i2=$(( i1 + 1 ))

        a=(`ls $path1/*.fastq`)

        r1=${a[@]:$i1:1}
        r2=${a[@]:$i2:1}
        r3=`echo $r2 | awk -F"/|_" '{print $9}'`

        # Remember to build the index first!

        mkdir $path2/$r3

        kallisto quant --index=$path0/genes.fna.kallisto_index --output=$path2/$r3 --plaintext $r1 $r2
```
### submit the job
```
sbatch quantifier.sh

ls -lh kallisto/t2A # Check presence and size (if it is 0 something went wrong) of an abundance table.
```
### Merge TPM tables into one
```
cd kallisto
ls -l */abundance.tsv | perl -ne 'chomp $_; if ($_ =~ /(\S+)\/abundance\.tsv/){print "\t$1"}' | perl -ne 'print "target_id\tlength$_\n"' > header.tsv
paste */abundance.tsv | cut -f 1,2,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105 > transcript_tpms_all_samples.tsv
cat header.tsv transcript_tpms_all_samples.tsv | grep -v "tpm" > transcript_tpms_all_samples.tsv2
mv transcript_tpms_all_samples.tsv2 kallisto_tpms_7x3.tsv
rm -f header.tsv
rm transcript_tpms_all_samples.tsv
```
