#!/bin/bash                                                                                                                        
#SBATCH -N 1
#SBATCH -J binning
#SBATCH -n 4
#SBATCH --partition=BIN420
	
module load metabat 

path1='/mnt/users/student54/BIN420/Monday/data/DNA/results/assembly/'

runMetaBat.sh -t 4 --verysensitive --unbinned $path1/contigs.fasta $path1/assembly_D1B.sorted_bam $path1/assembly_D2B.sorted_bam
