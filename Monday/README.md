# Metagenome assemby 
by Francesco  Delogu

The folder /mnt/SCRATCH/defr is a copy of what your folder should look like at the end of the practical session. There you can find sbatch files, data and results if you need them.


 Download and prepare the sra-toolkit

Connect to the Orion cluster
```
cd ~	# Go to your home directory if you were not there yet
mkdir scripts	# Create a new folder called "scripts"
cd scripts	# Enter folder "scripts"
wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.0/sratoolkit.2.10.0-ubuntu64.tar.gz	# Download the sra-toolkit
tar -xvzf sratoolkit.2.10.0-ubuntu64.tar.gz	# Uncompress the directory
cd sratoolkit.2.10.0-ubuntu64
cd bin
pwd	# Print current directory
```

Copy the printed path
```
vi ~/.bashrc	# Open the configuration bash file with the text editor vi
```
 Enter writing mode pressing the key "i"
 Add the line "export PATH=$PATH:" + the copied path
 Enter the writing mode pressing the key "ESC"
 Enter the command mode pressing the key ":"
 Save and quit typing "wq" and pressing "Enter"
 Close the connection with the cluster using "Ctrl+D"
 Open the connection with the cluster again
```
cd ~	# Go to your home directory if you were not there yet
module load perl	# Load perl for your user
fastq-dump --help	# Check if the command works
mkdir data
cd data

mkdir DNA
cd DNA

fastq-dump -X 5 -Z SRX3777359	# Print the first 5 mate pairs of reads form sample D1B on screen
```
If we can download the sequences:
```
fastq-dump -X 10000000 --split-3 SRX3777359 -O D1B --gzip	# Download the first 10000000 sequences of sample D1B and split them according to direction
fastq-dump -X 10000000 --split-3 SRX3777360 -O D2B --gzip	# Download the first 10000000 sequences of sample D2B and split them according to direction
```
 Else, copy the seqeunces from the cluster with a sbatch script (e.g. copier.sh):

```
cd ~

	#!/bin/bash

	#SBATCH -N 1
	#SBATCH -J copy_reads
	#SBATCH -n 1

	cd ~/DNA

	cp -r /mnt/SCRATCH/defr/DNA/D1B ./
	cp -r /mnt/SCRATCH/defr/DNA/D2B ./


sbatch copier.sh
```
Inspect sequence quality
```
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
```
Copy printed path
Open a shell on your machine (Ctrl+C on Ubuntu)
```
scp <YOUR USER>@orion.nmbu.no:<COPIED PATH>D1B/*html ./Desktop/
```
Enter your password
```
scp <YOUR USER>@orion.nmbu.no:<COPIED PATH>D2B/*html ./Desktop/
```
Enter your password



## Clean sequences
Build sbatch file (e.g. trimmer.sh):
```
	#!/bin/bash

	#SBATCH -N 1
	#SBATCH -J assembly
	#SBATCH -n 8

	dc ~/DNA
	module load trimmomatic

	trimmomatic -threads 8 PE D1B/SRR6820513_1.fastq.gz D1B/SRR6820513_2.fastq.gz trim/D1B_R1.fastq.gz trim/D1B_U1.fastq.gz trim/D1B_R2.fastq.gz trim/D1B_U2.fastq.gz SLIDINGWINDOW:10:28 MINLEN:100
	trimmomatic -threads 8 PE D2B/SRR6820512_1.fastq.gz D2B/SRR6820512_2.fastq.gz trim/D2B_R1.fastq.gz trim/D2B_U1.fastq.gz trim/D2B_R2.fastq.gz trim/D2B_U2.fastq.gz SLIDINGWINDOW:10:28 MINLEN:100

```

# Assemble metagenome

 Build sbatch file (e.g. assembler.sh):
```
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

```

# Compute coverage
Build sbatch file:
```
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
```
```
sbatch aligner.sh
```

# Metagenome Binning 
by Live Heldal Hagen Bin420 – November 2019


Summary
In this exercise, we will reconstruct genomes from metagenome assemblies (metagenome-assembled
genomes/MAG). We use the contigs from the metaSPAdes assembly an do a de-novo metagenome binning.  Since the results of different binning algorithms varies between samples, several programs (such as MetaBAT, MetaBAT2, MaxBin2 and CONCOCT) should optimally be tested and evaluated for each metagenome project. The results of multiple binning algorithms can furthermore be combined using DAS Tool to generate more accurate MAGs. Due to time restrictions, we will only be using MetaBAT in this exercise, which clusters the contigs using information such as nucleotide composition and read coverage. This clustering gives different genome bins ranging in size from a few thousand bp to almost complete genomes, depending on the sample. The MAGs, basically a list of FASTA files, need to be checked for completeness and contamination. We will do so by examining specific sets of marker genes for the corresponding clades using CheckM. An ultimate quality check maps the reconstructed MAGs onto the Spades assembly graph using Bandage (a demonstration will be given). There are several different ways to assign taxonomy to the metagenome, either to reads, to the contigs or to the MAGs. In this exercise, we will use MiGA Web (http://microbial-genomes.org/) to taxonomically classify the MAGs that we generated using MetaBAT. MiGA Web is a user-friendly online tool to make taxonomic inferences about genomes from isolates as well as from metagenomes. It uses a combination of ANI (Average Nucleotide Identity) and AAI (Average Amino-acid Identity) to taxonomically classify the genome sequence against sequences in its reference databases.  


1 PREPARATION
We will run both the binning and the quality check on the orion cluster: 
Connect to the Orion cluster using the terminal in RStudio (or via PuTTY) and go to your home directory (if you were not there yet) and create a new folder for the binning analysis in the results directory: 
>mkdir results/binning/


2 DE-NOVO BINNING WITH METABAT (Cluster)
De-novo contig grouping using Metabat requires only the contigs FASTA and BAM files, which
contain a mapping of reads onto the contigs. The BAM files are available in $PATH/TO/USER/results/assembly/
and were pre-calculated using the following commands: 

First, lets load and run MetaBat to see the available command line arguments: 
module load metabat
metabat -h

The version we are running is v0.26.3

MetaBAT requires the contig fasta file and sorted BAM files, which contain a mapping of reads onto the contigs. The mapping was done in the previous assembly section (using bowtie2 and samtools), and the output are called assembly_D1B.sorted_bam and assembly_D2B.sorted_bam


Make a new folder in the results directory:
cd /mnt/PATH/TO/USER/results/
mkdir binning/
cd binning/

Build sbatch file: 
>vi metabat_binning.sh

	#!/bin/bash                                                                                                                        
	#SBATCH -N 1
	#SBATCH -J binning
	#SBATCH -n 4
	#SBATCH --partition=BIN420
	
	module load metabat 

	path1='/mnt/PATH/TO/USER/results/assembly'

	runMetaBat.sh -t 4 --verysensitive --unbinned $path1/contigs.fasta $path1/assembly_D1B.sorted_bam $path1/assembly_D2B.sorted_bam

#end of sbatch

Run sbatch file: 
>sbatch metabat_binning.sh

This will generate a FASTA file for each bin (5 bins and one file of unbinned), with a rather long and complicated name. To make the downstream analysis easier, you could change this to a simpler name (such as out.1.fa, out.2.fa etc), and move them to a output folder (Please note that on newer version of metaBAT, you will be able to both specify output folder and a simpler name. Makes life easier!). 

Here is a step-by-step (the long way) around this. 
First we rename files using the "mv" command: 
>mv contigs.fasta.metabat-bins-_-t_4_--verysensitive_--unbinned.1.fa out.1.fa	#do this for all bins

Then we make the file called "bins", and move the binned fasta files to this folder:
>mkdir bins/
>mv out.* bins/  #this moves all files you have named out.somethingsomething to the bin folder, and we will use this as a input in the next step.

3 MAG QUALITY CHECKING WITH CHECKM (Cluster)
Remember, although reconstruction of MAGs can greatly improve data interpretation, binning can also be an important source of error. To check the quality of the MAGs, we will run CheckM. First, lets check out the version and the options: 

>module load checkm
>checkm -h

We will use the lineage_wf workflow with the following options:
threads: $threads; number of cpu cores
reduced-tree; for systems with less main memory
extension: fa; FASTA file suffix used by MetaBAT
We will then assess quality (contamination and completeness) through a table and a plot. We can choose between 9 different outputs, where the 9th is the most comprehensive. In this run, I use -o 2, which gives us a table. Other output formats can be found here: https://github.com/Ecogenomics/CheckM/wiki/Genome-Quality-Commands 

The input files are the bins in fasta format.
 

Build sbatch file: 
>vi checkM.sh
 
	#!/bin/sh

	#SBATCH -N 1
	#SBATCH -J checkm
	#SBATCH -n 4
	#SBATCH --partition=BIN420
	
	module load checkm

	path1='/mnt/PATH/TO/USER/results/binning/bins'
	
	checkm lineage_wf -t 4 -r -x fa $path1 $path1/checkm
	checkm qa -o 2 $path1/checkm/lineages.ms $path1 $path1/checkm_table
	checkm bin_qa_plot -x fa $path1/checkm $path1 $path1/checkm_plot

#end of sbatch 

Run sbatch: 
>sbatch checkM.sh

The output is a plot and a table with quality scores for each MAG. Take a look at the completeness and contamination columns. Are the MAGs of good quality?   

4 TAXONOMIC CLASSIFICATION (Web)
Download the 5 bins (FASTA formate) to your local computer (e.g. using the scp command). 
In your favorite web browser, go to http://microbial-genomes.org/ - the Microbial Genomes Atlas (MiGA) oneline, and create a user account (this is easy, and free!). After creating your account, you log in to MiGA and you first page should be the Dashboard. Here, the Query datasets will take you to your results, once your analysis are submitted. To submit your MAGs, go to Home, chose the right analysis (we will use NCBI Prok), click Upload genome and upload the FASTA file of each MAG. Remember to choose the right type of dataset, “popgenome”. Click “Upload new dataset” and wait for your results (but don’t hold your breath, this can take hours - lets check the results tomorrow)!

Video tutorials for an introduction to MiGA Oneline can be found here: https://manual.microbial-genomes.org/part3/web

Note: Although MiGA is easy to use, it might not be appropriate for datasets consisting of hundreds of MAGs. For this purpose, the GTDB-Tk software https://github.com/Ecogenomics/GTDBTk) is a promising alternative. 
