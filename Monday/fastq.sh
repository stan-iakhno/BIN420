#!/bin/bash
#SBATCH --ntasks=8               # Number of cores (CPUs, threads)
#SBATCH --job-name=fastq.stan         # Sensible name for the job
#SBATCH --nodes=1                # We *always* use 1 node
#SBATCH --partition=BIN420       # We *always* use this partition

## Software must be loaded before you start using it
module purge
module load fastqc

## Then you enter the command line(s) of your choice...
INDIR=/mnt/users/student54/BIN420/data/  ## data folder
                                                              ## output here

cd $INDIR
fastqc /DNA/*/*.fastq.gz