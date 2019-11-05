###Metagenome Binning Exercise
Bin420 – November 2019


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
