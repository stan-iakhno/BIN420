The following link contains the the walkthrough of the morning session:
```
file:///C:/Users/stia/OneDrive%20-%20Norwegian%20University%20of%20Life%20Sciences/PhD_courses/BIN420/BIN420_tuesday_morning.html

```

  BIN420 - Tuesday morning


        /Lars Snipen/


  1 Predicting genes with |prodigal|

From the Monday session you should have some contigs. You will use these
as input to this session. In the code below I use some pre-computed
contigs, and you may also use these if you for some reason did not (yet)
have the results from the assembly part.

First, start by creating a folder named |tuesday| in your |$HOME| folder. You can do this in RStudio, or from the Terminal by

|mkdir tuesday|

All scripts and results from today we put into this |$HOME/tuesday| folder.


    1.1 The shell script

We will use the |prodigal| software to predict genes from our contigs.

A simple shell script for running |prodigal| is as follows:
```
|#!/bin/bash
#SBATCH --ntasks=1               # Number of cores (CPUs, threads)
#SBATCH --job-name=prodigal      # Sensible name for the job
#SBATCH --nodes=1                # We *always* use 1 node
#SBATCH --partition=BIN420       # We *always* use this partition

## Software must be loaded before you start using it
module load prodigal

## Then you enter the command line(s) of your choice...
INDIR=/mnt/project/Courses/BIN420-2019/Day1_AssemblyAndBinning/assembly
OUTDIR=$HOME/tuesday/prodigal
GENOME_FILE=contigs.fasta

prodigal -a $OUTDIR/proteins.faa -d $OUTDIR/coding.ffn -f gff -o $OUTDIR/table.gff -p meta -i $INDIR/$GENOME_FILE|
```
Open a new text file in RStudio (*File - New File - Text File*) and copy
the code above into it. Save the file as |prodigal.sh| in your |$HOME/tuesday| folder.

You may need to edit the path in |INDIR| to match where you output the |spades| results from Monday. You may also want to edit the |OUTDIR| folder. Note that this folder must exist /before/ you run |prodigal|, i.e. we need to create this folder as well.

Once you are happy with your script, submit it to SLURM. In your
Terminal window, you may need to jump to the |tuesday| folder first:
```
cd $HOME/tuesday
```
To verify you are actually in the correct folder, run a listing of its
content:
```
ls -als
```
and verify the file |prodigal.sh| is actually listed. If not, you are not either not in the |$HOME/tuesday| folder, or you saved the shell script to some other folder. Next,
create the |prodigal| folder inside |$HOME/tuesday|.

From the |$HOME/tuesday| folder you |sbatch|:

|sbatch prodigal.sh|

To verify your job is running, inspect the queue:

|squeue -u username|

where you enter your actual username instead of |username|. The |prodgal| is superfast, and should be done within a couple of minutes.


    1.2 The output

We will now use R to process some of the output from |prodigal|. By the options we used above, we get three output files:

  * The |.ffn| file. This is a fasta file with the predicted coding genes,
    i.e. it contains only DNA sequences.
  * The |.faa| file. This is a fasta file with the protein sequences we get by
    translating the genes in the |.ffn| file. This is the file we are interested in when we start
    annotating shortly.
  * The |.gff| file. This is a table with one row for each gene, containing some
    information in the GFF-format.

The |prodigal| software will sometimes predict genes with a stop codon inside (when
using the |-p meta| option). This is most easily seen in the protein sequences, since the
symbol for a stop codon after translation is |*|. Proteins should only have a |*| at the end. Later we will use the |interproscan| software to annotate the genes/proteins, and this software does not
tolerate |*| in the input at all. Thus, we need to

  * Remove all |*| from the end of the protein sequences
  * Discard all genes with an |*| inside

Let us make a small R script for doing this:

|library(tidyverse)
library(microseq)
setwd("~/tuesday")

## Reading files
ffn.tbl <- readFasta("prodigal/coding.ffn")                  ## reads coding file
gff.tbl <- readGFF("prodigal/table.gff")                     ## reads GFF file
faa.tbl <- readFasta("prodigal/proteins.faa") %>%            ## reads protein file
  mutate(Sequence = str_remove(Sequence, "\\*$"))            ## remove * at the ends

## Filtering out genes with stop codon inside (if any)
has.stop.inside <- str_detect(faa.tbl$Sequence, "\\*")       ## detects proteins with *
faa.tbl %>% 
  filter(!has.stop.inside) %>%                               ## discards proteins with *
  writeFasta(out.file = "prodigal/proteins_nostops.faa")     ## write to new fasta file
ffn.tbl %>% 
  filter(!has.stop.inside) %>%                               ## discard the same genes as above
  writeFasta(out.file = "prodigal/coding_nostops.ffn")       ## writes to new fasta file
gff.tbl %>% 
  filter(!has.stop.inside) %>%                               ## same filtering as above
  writeGFF(out.file = "prodigal/table_nostops.gff")|

|## [1] "gff.table written to prodigal/table_nostops.gff"|

Make a new R script (*File - New File - R Script*) and copy the code
into it. Save it in your |$HOME/tuesday| folder.

Small R scripts we run directly in RStudio. In the upper right corner of
the Source (editor) window in RStudio, you find a Source button. Press
this to run the script.

We may inspect |has.stop.inside| to see if any proteins had a stop codon inside, but the new files we
created are guaranteed to be without such problems.

There is one last piece of code that we want to add to the script above.
Let us output also a small subset of 100 proteins. This we need below,
to test the software |interproscan| that we will use for annotation. The annotations by |interproscan| takes some time, and it is nice to have a mini data set for testing
first. Add these lines to the script above, and re-run:

|faa.tbl %>% 
  filter(!has.stop.inside) %>% 
  slice(1:100) %>%
  writeFasta(out.file = "prodigal/proteins_nostops_mini.faa")
ffn.tbl %>% 
  filter(!has.stop.inside) %>% 
  slice(1:100) %>%
  writeFasta(out.file = "prodigal/coding_nostops_mini.ffn")
gff.tbl %>% 
  filter(!has.stop.inside) %>% 
  slice(1:100) %>%
  writeGFF(out.file = "prodigal/table_nostops_mini.gff")|





  2 Annotation with |interproscan|

The |prodigal| predicts genes, but says nothing about what these genes are. For that
purpose we now use |interproscan|. The input is a fasta file of /protein sequences/ (the |.faa| file), and these must not contain any stop codons (i.e. no |*|).

You will hear more about annotations after lunch, but since the |interproscan| takes some time to complete, we need to start it as soon as possible.


    2.1 The shell script

Here is a small shell script for running |interproscan|, using our mini protein file from above:

|#!/bin/bash
#SBATCH --ntasks=1               # Number of cores (CPUs, threads)
#SBATCH --job-name=iprscan      # Sensible name for the job
#SBATCH --nodes=1                # We *always* use 1 node
#SBATCH --partition=BIN420       # We *always* use this partition

## Software must be loaded before you start using it
module load interproscan

## Then you enter the command line(s) of your choice...
TMPDIR=$HOME/tuesday/tmp
INFILE=$HOME/tuesday/prodigal/proteins_nostops_mini.faa
OUT_PREFIX=$HOME/tuesday/iprscan/proteins_nostops_mini

interproscan.sh --input $INFILE --disable-precalc --applications CDD,HAMAP,Pfam,TIGRFAM --formats GFF3,TSV --iprlookup --goterms --pathways --cpu 8  --tempdir $TMPDIR --output-file-base $OUT_PREFIX
|

Make a new shell script, and copy this code into it, and save it into your |$HOME/tuesday| folder as |iprscan.sh|. Let us have a brief look at the command line.

The input and output is straightforward, and note that we also supply a
folder for |interproscan| to store temporary results (the |tmp| folder). Notice also we output to a subfolder named |$HOME/tuesday/iprscan| here, create this folder right away.

When annotating proteins, many proteins are similar across virtually all
organisms, and it seems like a waste of resources to re-compute their
annotation each time. For this reason the |interproscan| has a look-up service, where some results are downloaded to speed it
up. However, this service frequently fails, and by the option |--disable-precalc| we turn it off and do all calculation here and now.

The annotations are obtained by comparing the proteins to various
databases of already annotated proteins. the option |--applications CDD,HAMAP,Pfam,TIGRFAM| specifies which databases to look into. There are many more databases
you may specify, see the |interproscan| wiki pages <https://github.com/ebi-pf-team/interproscan/wiki> for more
details on this. The four databases we selected here is a compromise
between speed and comprehensive annotations. By using all databases |interproscan| will be /very slow/.

The |--formats| indicate we get results in two formats, both are text files with
results in columns (tables).

The options |--iprlookup --goterms --pathways| are important. This turns on the listing of gene onthology (GO) terms
and KEGG pathways information in the output. This is something we would
like to have when we start the other omics-studies later.

Again, we run this by |sbatch|ing to SLURM:

|sbatch iprscan.sh|

Since we use the mini file with only 25 proteins, this will take a short
time to complete (few minutes). Later we will use the full protein file,
and then it will take a long time.


    2.2 The output

In the |iprscan| folder you should find two result files, with extensions |.gff3| and |.tsv|. Both files are text files with results in a table format. This time
we focus on the |.gff3| file only. Both tables have a variable number of rows for each gene,
since the number of database hits varies. Also, much of the information
we seek is ‘buried’ inside longer texts We would like to have a table where

  * Each row corresponds to one gene
  * Only columns of some interest are displayed (GO-terms, KEGG-terms etc)

Here is some R code for converting the |.gff3| file into such a table:

|library(tidyverse)
library(microseq)
setwd("~/tuesday")

## Edit filenames here
infile <- "iprscan/proteins_nostops_mini.gff3"
outfile <- "annotations_mini.txt"

## Reading GFF file, and extracting information
readGFF(infile) %>% 
  select(Seqid, Database = Source, Attributes) %>% 
  filter(Database != ".") %>% 
  mutate(Database.accession = str_remove_all(str_extract(Attributes, "Name=.+?;"), "Name=|;")) %>% 
  mutate(InterPro.accession = str_extract(Attributes, "IPR[0-9]+")) %>% 
  mutate(Description = str_remove_all(str_extract(Attributes, "signature_desc=.+?;"), "signature_desc=|;")) %>% 
  mutate(GO.terms = str_extract_all(Attributes, "GO:[0-9]+")) %>% 
  mutate(GO.terms = lapply(GO.terms, function(x){if(length(x) == 0) return(NA) else return(x)})) %>% 
  mutate(KEGG.terms = str_extract_all(Attributes, "KEGG:.+?\"")) %>% 
  mutate(KEGG.terms = lapply(KEGG.terms, str_remove, "\"")) %>% 
  mutate(KEGG.terms = lapply(KEGG.terms, function(x){if(length(x) == 0) return(NA) else return(x)})) %>%
  group_by(Seqid) %>% 
  summarize(Database = str_c(str_replace_na(Database), collapse = ","),
            Database.accession = str_c(str_replace_na(Database.accession), collapse = ","),
            InterPro.accession   = str_c(str_replace_na(unique(InterPro.accession)), collapse = ","),
            Description = str_c(str_replace_na(Description), collapse = ","),
            GO.terms   = str_c(str_replace_na(unique(unlist(GO.terms))), collapse = ","),
            KEGG.terms = str_c(str_replace_na(unique(unlist(KEGG.terms))), collapse = ",")) %>% 
  mutate_all(str_remove_all, ",NA|NA,") -> annotation.tbl
write_delim(annotation.tbl, path = outfile, delim = "\t")|

Make a new R script, copy this code, and save it as |gff2annotations.R| in your |$HOME/tuesday| folder. Run the script. A small file named |annotations_mini.txt| should appear in your |$HOME/tuesday| folder. We may also inspect the table directly in RStudio.

Once you have re-run the |iprscan.sh|, using the full set of proteins as input, you use this R script again
to convert the corresponding |.gff3| file to a similar table, useful for the downstream analysis.



  3 Extra - inspecting results in R

This part is just something we look into if we should have some extra
time before lunch.


    3.1 Predicted genes not annotated

Let us have a look at the results from the mini data set. We read the
annotations into R again, as if we started a new R session:

|library(tidyverse)
library(microseq)

## Read annotations_mini.txt
annot.tbl <- suppressMessages(read_delim("~/tuesday/annotations_mini.txt", delim = "\t"))|

We notice there are 84 rows, i.e. 16 out of the 100 genes in the mini
data set were not annotated. Let us locate these genes in the |prodigal| output. We read the GFF table produced by |prodigal|:

|readGFF("~/tuesday/prodigal/table_nostops_mini.gff") %>% 
  mutate(ID = word(Attributes, 1, 1, sep = ";")) %>% 
  mutate(Tag = str_c(Seqid, str_remove(ID, "ID=[0-9]+"))) %>% 
  select(-ID) %>% 
  arrange(desc(Score)) -> prdgl.tbl|

Here we also made two new columns, the |ID| and the |Tag|, then we de-selected the former. The |Tag| column is needed shortly. We also sorted the genes by their |Score|. This is the score given by |prodigal|, and the higher score value, the more confident we should be in this
really being a gene (according to the |prodigal| algorithm).

Let us now add another column, indicating if the gene was annotated.
This is where we need the |Tag| information, because we need to match this to the |Seqid| column of the |annot.tbl| (which is /not/ the same as the |Seqid| column in the |prdgl.tbl|):

|prdgl.tbl %>% 
  mutate(Annotated = Tag %in% annot.tbl$Seqid) %>% 
  mutate(Length = abs(Start - End) + 1) -> prdgl.tbl|

We also added the column |Length|, which is simply the length of the gene (in bases, and including the
stop codon). From this table, we make a plot

|p1 <- ggplot(prdgl.tbl) +
  geom_point(aes(x = Score, y = Length, color = Annotated), size = 3, alpha = 0.5) +
  labs(x = "Prodigal score", y = "Gene length")
print(p1)|

We clearly see the un-annotated genes tend to have lower scores. Also,
this score is highly correlated with gene length.


    3.2 All predicted genes

The |prodigal| also outputs a |.gff| file, which is a table with one row for each predicted gene (protein).
Let us read this into R, the full table this time, and have a look at
its content:

|library(tidyverse)
library(microseq)

readGFF("~/tuesday/prodigal/table_nostops.gff") %>% 
  mutate(Contig.length = as.numeric(word(Seqid, 4, 4, sep = "_"))) %>% 
  mutate(Contig.coverage = round(as.numeric(word(Seqid, 6, 6, sep = "_")))) %>% 
  separate(Attributes, into = c("ID", "Partial", "Start.codon", "RBS.motif",
                                "RBS.spacer", "GC", "Confidence", "Total.score",
                                "C.score", "S.score", "rscore", "uscore", "tscore","empty"),
           sep = ";") %>% 
  mutate(Tag = str_c(Seqid, str_remove(ID, "ID=[0-9]+"))) %>% 
  select(-c(ID, Partial, RBS.spacer, rscore, uscore, tscore, empty)) %>% 
  mutate_at(9:15, str_remove, "^.+=") %>% 
  mutate_at(11:15, as.numeric) -> prdgl.tbl|

After reading the file, we extracted the length and K-mer coverage for
each contig, which is always part of the contig headers from |spades|. Next, we split the |Attributes| column (last column) of the GFF format into new columns by the
semicolon separator. We also make a new column |Tag| like we did in the previous exercise. Then we discard some of the
columns, and modify some of the remaining (removing som text and
converting to numeric). Open the table in the RStudio viewer for inspection.

We have more than 22 000 predicted genes, from a long range of contigs.
Let us plot the |Total.score| for each gene against the |Contig.length|, i.e. the length of the contig where the gene is found:

|p2 <- ggplot(prdgl.tbl, aes(x = Contig.length, y = Total.score)) +
  geom_point(alpha = 0.3) +
  geom_smooth() +
  scale_x_log10() + scale_y_log10()
print(p2)|

|## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'|

Note that both axes are log-transformed. We clearly see that genes from
shorter contigs tend to have a lower |prodigal| score, which means they are less reliable. There are also a large
number of short contigs, as usual from Illumina data:

|prdgl.tbl %>% 
  distinct(Seqid, .keep_all = T) %>% 
  ggplot() +
  geom_histogram(aes(x = Contig.length), bins = 50) +
  scale_x_log10() +
  labs(x = "Contig length", y = "Number of contigs") -> p3
print(p3)|

Once we have the annotation table for all genes (from |interproscan|) we may again mark which genes were annotated, and then add some color
to the plot from above:

|annot.tbl <- suppressMessages(read_delim("~/tuesday/annotations.txt", delim = "\t")) ## once this
                                                                                     ## is available...
prdgl.tbl %>% 
  mutate(Annotated = Tag %in% annot.tbl$Seqid) -> prdgl.tbl
p4 <- ggplot(prdgl.tbl, aes(x = Contig.length, y = Total.score)) +
  geom_point(aes(color = Annotated), alpha = 0.2) +
  geom_smooth() +
  scale_x_log10() + scale_y_log10() +
  facet_wrap(~Annotated, nrow = 1)
print(p4)|

|## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'|

The main trend is as expected, shorter contigs produce less reliable
gene predictions, and a large number of un-annotated genes, even if
there are a good number of exceptions. Having longer contigs would most
likely have improved the results from gene prediction and annotations.

