#!/bin/bash
#SBATCH --ntasks=8               # Number of cores (CPUs, threads)
#SBATCH --job-name=spades.stan         # Sensible name for the job
#SBATCH --nodes=1                # We *always* use 1 node
#SBATCH --partition=BIN420       # We *always* use this partition

## Software must be loaded before you start using it
module purge
module load spades

## Then you enter the command line(s) of your choice...
INDIR=/mnt/users/student54/BIN420/data/DNA  ## data folder
                                                              ## output here

cd $INDIR/

spades.pspades.py --meta -t 8 -m 70 -1 $INDIR/D1B/trim/D1B_R1.fastq.gz \
	-2 $INDIR/D1B/trim/D1B_R2.fastq.gz \
	-1 $INDIR/D2B/trim/D2B_R1.fastq.gz \
	-2 $INDIR/D2B/trim/D2B_R2.fastq.gz \
	-s $INDIR/D1B/trim/D1B_U1.fastq.gz \
	-s $INDIR/D1B/trim/D1B_U2.fastq.gz \
	-s $INDIR/D2B/trim/D2B_U1.fastq.gz \
	-s $INDIR/D2B/trim/D2B_U2.fastq.gz \
	-o results/assembly