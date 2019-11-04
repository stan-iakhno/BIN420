#!/bin/bash
#SBATCH --ntasks=1               # Number of cores (CPUs, threads)
#SBATCH --job-name=first.stan         # Sensible name for the job
#SBATCH --nodes=1                # We *always* use 1 node
#SBATCH --partition=BIN420       # We *always* use this partition

## Software must be loaded before you start using it
module load seqtk

## Then you enter the command line(s) of your choice...
INDIR=/mnt/project/Courses/BIN420-2019/Day1_AssemblyAndBinning/data/trim  ## data folder
OUTDIR=$HOME                                                              ## output here
IN_FASTQ=D1B_U1.fastq.gz     ## input fastq file
OUT_FASTA=small.fasta        ## create this fasta file
OUT_TXT=counts.txt           ## create this text file


zcat $INDIR/$IN_FASTQ | head -n400 | seqtk seq -A > $OUTDIR/$OUT_FASTA
seqtk comp $OUTDIR/$OUT_FASTA > $OUTDIR/$OUT_TXT