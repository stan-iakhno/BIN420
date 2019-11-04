
git commit -m "first"
git remote add origin https://github.com/stan-iakhno/BIN420.git
git push -u origin master
## The folder /mnt/SCRATCH/defr is a copy of what your folder should look like at the end of the practical session. There you can find sbatch files, data and results if you need them.


## Download and prepare the sra-toolkit

# Connect to the Orion cluster
cd ~	# Go to your home directory if you were not there yet
mkdir scripts	# Create a new folder called "scripts"
cd scripts	# Enter folder "scripts"
wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.0/sratoolkit.2.10.0-ubuntu64.tar.gz	# Download the sra-toolkit
tar -xvzf sratoolkit.2.10.0-ubuntu64.tar.gz	# Uncompress the directory
cd sratoolkit.2.10.0-ubuntu64
cd bin
pwd	# Print current directory
# Copy the printed path
vi ~/.bashrc	# Open the configuration bash file with the text editor vi
# Enter writing mode pressing the key "i"
# Add the line "export PATH=$PATH:" + the copied path
# Enter the writing mode pressing the key "ESC"
# Enter the command mode pressing the key ":"
# Save and quit typing "wq" and pressing "Enter"
# Close the connection with the cluster using "Ctrl+D"
# Open the connection with the cluster again
cd ~	# Go to your home directory if you were not there yet
module load perl	# Load perl for your user
fastq-dump --help	# Check if the command works
mkdir data
cd data

mkdir DNA
cd DNA

fastq-dump -X 5 -Z SRX3777359	# Print the first 5 mate pairs of reads form sample D1B on screen

## If we can download the sequences:

fastq-dump -X 10000000 --split-3 SRX3777359 -O D1B --gzip	# Download the first 10000000 sequences of sample D1B and split them according to direction
fastq-dump -X 10000000 --split-3 SRX3777360 -O D2B --gzip	# Download the first 10000000 sequences of sample D2B and split them according to direction

## Else, copy the seqeunces from the cluster with a sbatch script (e.g. copier.sh):

cd ~

	#!/bin/bash

	#SBATCH -N 1
	#SBATCH -J copy_reads
	#SBATCH -n 1

	cd ~/DNA

	cp -r /mnt/SCRATCH/defr/DNA/D1B ./
	cp -r /mnt/SCRATCH/defr/DNA/D2B ./


sbatch copier.sh

## Inspect sequence quality

cd ~	# Go to your home directory if you were not there yet
module load fastqc
cd data
cd D1B
fastqc *.fastq	# Analyze D1B
cd ..
cd D2B
fastqc *.fastq	# Analyze D2B
cd ..
pwd
# Copy printed path
# Open a shell on your machine (Ctrl+C on Ubuntu)
scp <YOUR USER>@orion.nmbu.no:<COPIED PATH>D1B/*html ./Desktop/
# Enter your password
scp <YOUR USER>@orion.nmbu.no:<COPIED PATH>D2B/*html ./Desktop/
# Enter your password


#####
## Clean sequences
#####

# Build sbatch file (e.g. trimmer.sh):

	#!/bin/bash

	#SBATCH -N 1
	#SBATCH -J assembly
	#SBATCH -n 8

	dc ~/DNA
	module load trimmomatic

	trimmomatic -threads 8 PE D1B/SRR6820513_1.fastq.gz D1B/SRR6820513_2.fastq.gz trim/D1B_R1.fastq.gz trim/D1B_U1.fastq.gz trim/D1B_R2.fastq.gz trim/D1B_U2.fastq.gz SLIDINGWINDOW:10:28 MINLEN:100
	trimmomatic -threads 8 PE D2B/SRR6820512_1.fastq.gz D2B/SRR6820512_2.fastq.gz trim/D2B_R1.fastq.gz trim/D2B_U1.fastq.gz trim/D2B_R2.fastq.gz trim/D2B_U2.fastq.gz SLIDINGWINDOW:10:28 MINLEN:100


#####
## Assemble metagenome
#####

# Build sbatch file (e.g. assembler.sh):

	#!/bin/bash                                                                                                                        

	#SBATCH -N 1
	#SBATCH -J assembly
	#SBATCH -n 8

	module load spades

	data_path=$HOME'/DNA/trim'

spades.py --meta -t 8 -m 70 -1 $INDIR/D1B/trim/D1B_R1.fastq.gz 
	-2 $INDIR/D1B/trim/D1B_R2.fastq.gz 
	-1 $INDIR/D2B/trim/D2B_R1.fastq.gz 
	-2 $INDIR/D2B/trim/D2B_R2.fastq.gz 
	-s $INDIR/D1B/trim/D1B_U1.fastq.gz 
	-s $INDIR/D1B/trim/D1B_U2.fastq.gz 
	-s $INDIR/D2B/trim/D2B_U1.fastq.gz 
	-s $INDIR/D2B/trim/D2B_U2.fastq.gz 
	-o results/assembly

sbatch assembler.sh


#####
## Compute coverage
#####

# Build sbatch file:

	#!/bin/bash                                                                                                                        

	#SBATCH -N 1
	#SBATCH -J align
	#SBATCH -n 8

	module load bowtie2
	module load samtools

	data_path='~/assembly/corrected'
	assembly_path='~/assembly'
	out_path='~/assembly'

	mkdir $assembly_path/indices

	bowtie2-build -f $assembly_path/contigs.fasta $assembly_path/indices/contigs

	bowtie2 -p 16 -x $assembly_path/indices/contigs -1 $data_path/D1B_R1.fastq.00.0_0.cor.fastq.gz -2 $data_path/D1B_R2.fastq.00.0_0.cor.fastq.gz --very-sensitive -X 400 -I 180 | samtools view -bS > $out_path/assembly_D1B.bam
	samtools sort -l 9 -O bam -o  $out_path/assembly_D1B.sorted_bam $out_path/assembly_D1B.bam

	bowtie2 -p 16 -x $assembly_path/indices/contigs -1 $data_path/D2B_R1.fastq.00.0_1.cor.fastq.gz -2 $data_path/D2B_R2.fastq.00.0_1.cor.fastq.gz --very-sensitive -X 400 -I 180 | samtools view -bS > $out_path/assembly_D2B.bam
	samtools sort -l 9 -O bam -o  $out_path/assembly_D2B.sorted_bam $out_path/assembly_D2B.bam

sbatch aligner.sh

