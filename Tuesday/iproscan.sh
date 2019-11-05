#!/bin/bash
#SBATCH --ntasks=8               # Number of cores (CPUs, threads)
#SBATCH --job-name=iprscan      # Sensible name for the job
#SBATCH --nodes=1                # We *always* use 1 node
#SBATCH --partition=BIN420       # We *always* use this partition

## Software must be loaded before you start using it
module load interproscan

## Then you enter the command line(s) of your choice...
TMPDIR=$HOME/BIN420/Tuesday/prodigal/tmp
INFILE=$HOME/BIN420/Tuesday/prodigal/proteins_nostops_mini.faa
OUT_PREFIX=$HOME/BIN420/Tuesday/iprscan/prodigal&proteins_nostops_mini

interproscan.sh --input $INFILE --disable-precalc --applications CDD,HAMAP,Pfam,TIGRFAM --formats GFF3,TSV --iprlookup --goterms --pathways --cpu 8  --tempdir $TMPDIR --output-file-base $OUT_PREFIX
