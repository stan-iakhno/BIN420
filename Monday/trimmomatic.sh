#!/bin/bash
#SBATCH --ntasks=8               # Number of cores (CPUs, threads)
#SBATCH --job-name=trim.stan         # Sensible name for the job
#SBATCH --nodes=1                # We *always* use 1 node
#SBATCH --partition=BIN420       # We *always* use this partition

## Software must be loaded before you start using it
module purge
module load trimmomatic

## Then you enter the command line(s) of your choice...
INDIR=/mnt/users/student54/BIN420/data/DNA  ## data folder
                                                              ## output here

cd $INDIR/D1B/
trimmomatic PE -threads 8 SRR6820513_1.fastq.gz SRR6820513_2.fastq.gz trim/D1B_R1.fastq.gz trim/D1B_U1.fastq.gz trim/D1B_R2.fastq.gz trim/D1B_U2.fastq.gz SLIDINGWINDOW:10:28 MINLEN:100

cd $INDIR/D2B/
trimmomatic PE -threads 8 SRR6820512_1.fastq.gz SRR6820512_2.fastq.gz trim/D2B_R1.fastq.gz trim/D2B_U1.fastq.gz trim/D2B_R2.fastq.gz trim/D2B_U2.fastq.gz SLIDINGWINDOW:10:28 MINLEN:100