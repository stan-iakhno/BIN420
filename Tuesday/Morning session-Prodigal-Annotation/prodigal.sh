#!/bin/bash
#SBATCH --ntasks=1               # Number of cores (CPUs, threads)
#SBATCH --job-name=prodigal      # Sensible name for the job
#SBATCH --nodes=1                # We *always* use 1 node
#SBATCH --partition=BIN420       # We *always* use this partition

## Software must be loaded before you start using it
module load prodigal

## Then you enter the command line(s) of your choice...
INDIR=/mnt/users/student54/BIN420/Tuesday
OUTDIR=$HOME/BIN420/Tuesday/prodigal
GENOME_FILE=contigs.fasta

prodigal -a $OUTDIR/proteins.faa \
  -d $OUTDIR/coding.ffn \
  -f gff -o $OUTDIR/table.gff \
  -p meta -i $INDIR/$GENOME_FILE
