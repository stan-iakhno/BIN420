#!/bin/bash
#SBATCH --ntasks=16             # Number of cores (CPUs, threads)
#SBATCH --job-name=align.stan         # Sensible name for the job
#SBATCH --nodes=1                # We *always* use 1 node
#SBATCH --partition=BIN420       # We *always* use this partition

## Software must be loaded before you start using it
module purge
module load bowtie2
module load samtools

## Then you enter the command line(s) of your choice...


data_path='/mnt/users/student54/BIN420/data/DNA/results/'
assembly_path='/mnt/users/student54/BIN420/data/DNA/results/assembly'
out_path='/mnt/users/student54/BIN420/data/DNA/results/assembly'

mkdir $assembly_path/indices

bowtie2-build -f $assembly_path/contigs.fasta $assembly_path/indices/contigs

bowtie2 -p 16 -x $assembly_path/indices/contigs -1 $data_path/D1B_R1.fastq.00.0_0.cor.fastq.gz -2 $data_path/D1B_R2.fastq.00.0_0.cor.fastq.gz --very-sensitive -X 400 -I 180 | samtools view -bS > $out_path/assembly_D1B.bam
samtools sort -l 9 -O bam -o  $out_path/assembly_D1B.sorted_bam $out_path/assembly_D1B.bam

bowtie2 -p 16 -x $assembly_path/indices/contigs -1 $data_path/D2B_R1.fastq.00.0_1.cor.fastq.gz -2 $data_path/D2B_R2.fastq.00.0_1.cor.fastq.gz --very-sensitive -X 400 -I 180 | samtools view -bS > $out_path/assembly_D2B.bam
samtools sort -l 9 -O bam -o  $out_path/assembly_D2B.sorted_bam $out_path/assembly_D2B.bam