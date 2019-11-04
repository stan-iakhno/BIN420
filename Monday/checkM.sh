#!/bin/sh
#SBATCH -N 1
#SBATCH -J checkm
#SBATCH -n 4
#SBATCH --partition=BIN420
	
module load checkm

path1='/mnt/users/student54/BIN420/Monday/data/DNA/results/bins'
	
checkm lineage_wf -t 4 -r -x fa $path1 $path1/checkm
checkm qa -o 2 $path1/checkm/lineages.ms $path1 $path1/checkm_table
checkm bin_qa_plot -x fa $path1/checkm $path1 $path1/checkm_plot