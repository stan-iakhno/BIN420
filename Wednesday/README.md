# Quantifying abundances of transcripts using Kallisto
by Francesco Delongu
## Prepare transcripts database

### Go to the folder where the bins are
```
for i in *.fa; do j=`echo $i | awk -F"." '{print $2}'`; awk -v var=$j '/>/{print substr($1, 2, length($1))"\t"var}' $i >> bin_map.txt; done
awk '{a++; if(a>1) print $1}' annotations.txt | awk -F"_" '{a[$0"\t"$1"_"$2"_"$3"_"$4"_"$5"_"$6]=1}END{for(i in a)print i}' > ORF_contig.txt
awk -F"\t" 'NR==NRF{a[$1]=$2; next} {print $0"\t"a[$2]}' bin_map.txt ORF_contig.txt > ORF_contig_bin.txt

awk '{if($1~/^>/){print header"\n"seq; header=$1; seq=""} else seq=seq""$1}END{print header"\n"seq}' nostops.ffn | awk '{a++; if(a>2)print $0}' > nostops_2lines.ffn 
awk 'NR==FNR{a[">"$1]=1; next} {if(a[$1]==1){print $1; getline; print $1}}' ORF_contig_bin.txt nostops_2lines.ffn > genes.fna

module add kallisto
kallisto index -i genes.fna.kallisto_index --make-unique genes.fna

```
### Quantify transcripts
```
cd ~
mkdir kallisto
```
### Build sbatch file (e.g. quantifier.sh):

```
	#!/bin/bash

        #SBATCH -N 1
        #SBATCH -J kallisto
        #SBATCH -c 8
        #SBATCH -a 1-21

        path1='/mnt/project/Courses/BIN420-2019/Day3_MTandNetAnalysis/data/'
        path0='/mnt/users/student54/BIN420/Wednesday'
        path2='/mnt/users/student54/kallisto'

        module add kallisto

        i=$(( SLURM_ARRAY_TASK_ID - 1 ))

        i1=$(( i * 2 ))
        i2=$(( i1 + 1 ))

        a=(`ls $path1/*.fastq`)

        r1=${a[@]:$i1:1}
        r2=${a[@]:$i2:1}
        r3=`echo $r2 | awk -F"/|_" '{print $9}'`

        # Remember to build the index first!

        mkdir $path2/$r3

        kallisto quant --index=$path0/genes.fna.kallisto_index --output=$path2/$r3 --plaintext $r1 $r2
```
### submit the job
```
sbatch quantifier.sh

ls -lh kallisto/t2A # Check presence and size (if it is 0 something went wrong) of an abundance table.
```
### Merge TPM tables into one
```
cd kallisto
ls -l */abundance.tsv | perl -ne 'chomp $_; if ($_ =~ /(\S+)\/abundance\.tsv/){print "\t$1"}' | perl -ne 'print "target_id\tlength$_\n"' > header.tsv
paste */abundance.tsv | cut -f 1,2,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105 > transcript_tpms_all_samples.tsv
cat header.tsv transcript_tpms_all_samples.tsv | grep -v "tpm" > transcript_tpms_all_samples.tsv2
mv transcript_tpms_all_samples.tsv2 kallisto_tpms_7x3.tsv
rm -f header.tsv
rm transcript_tpms_all_samples.tsv
```





---
title: "MT_analysis"
author: "Francesco Delogu"
date: "September 6, 2017"
output: html_document
---

# 0 Getting started

RNA-seq expression data from a simplistic microbial community. 7 time points with 3 replicates each.

## 0.1 Loading libraries and data

Load all the libraries and data required. If the libraries are not installed you can install them with the command:
> install.packages("YOUR_PACKAGE")
or via RStudio using the button "Tools" -> "Install Packages..." and then typing the package in the box.

If the package is from Bioconductor you need to install it with the following commands:
> source("http://bioconductor.org/biocLite.R")
> biocLite("YOUR_PACKAGE") 

If you have a doubt about a function or a package use the command:
> ?YOUR_FUNCTION
if no information is provided, follow the advice from R and use the command:
> ??YOUR_FUNCTION
if still no information has been retrieved it's time to ask a human (or google).

```{r load libraries, message=F, warning=F, results='hide'}

library(Biobase)
library(WGCNA)
library(matrixStats)
library(tidyverse)
library(gplots)
library(corrplot)
library(GOstats)
library(GSEABase)
library(AnnotationDbi)
library(BiasedUrn)

options(stringsAsFactors = FALSE)

#setwd("/mnt/SCRATCH/guestN") # Put your path here!

replicates <- read.csv('Day3_MTandNetAnalysis/kallisto_tpms_7x3.tsv', header = T, row.names = 1, sep = '\t') # Change it with your path and file
datRep1 <- as.data.frame(t(replicates[-1]))

ann.df <- read_tsv('Day2_GenesAndAnnotations/annotations.txt') # Change it with your path and file

bin.df <- read_tsv('Day3_MTandNetAnalysis/ORF_contig_bin.txt', col_names=c("ORF", "contig", "bin")) %>% # Change it with your path and file
  filter(bin!="")

```

## 0.2 - Raw data exploration:

How many genes have an expression signal?

```{r number of genes, eval=FALSE}

dim(datRep1)

```

What does an expression profile from this dataset look like?

```{r first gene plot, out.width=1000, eval=FALSE}

plot(datRep1[,10000])

```

Which is the general structure of the dataset?

```{r structure, results='hide'}

str(datRep1)

```

# 1 Filtering and exploration of expression data

## 1.1 Signal filter

RNA-seq data are ususally noisy and we want to remove genes with no signal and the ones with outlier profiles. Therefore we use the goodSamplesGenes function (how does it work?) to filter out all the not desired genes.

```{r filtering, results='hide'}

?goodSamplesGenes

datRep1[datRep1 == 0] <- NA
gsg = goodSamplesGenes(datRep1, verbose = 3, minNSamples = 2)
if (!gsg$allOK){
  if (sum(!gsg$goodGenes)>0)
    #printFlush(paste("Removing genes:", paste(names(datRep1)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datRep1)[!gsg$goodSamples], collapse = ", ")))
  datExpr0 = datRep1[gsg$goodSamples, gsg$goodGenes]
} else{
  datExpr0 = datRep1
}

datExpr0[is.na(datExpr0)] <- 0

```

Try to change the parameter(s) in the "goodSamplesGenes" function. Do you have some first consideration about the expression the signal?

## 1.2 Visualize the samples

### 1.2.1 Data transformation

Before proceding to the main analysis transform the data with the log2 function. Of course don't try to compute log(0).
Log transformation is a "must have" in RNA data analysis: why and what biological feature(s) is(are) behind this reason?

```{r data to log2}

datExpr1 <- log2(datExpr0+1)

```

### 1.2.2 Replicate correlation

How strong is the variability among replicates from the same time point? We can have a first look computing the correlation among the replicates. What kind of correlation does the function "cor" compute by default? Is this a sound method to show variablity?

```{r replicate correlation, out.width=1000, results='hide'}

?stats::cor # "stats::" forces R to use the cor function from the library "stats"
?corrplot

for(i in (0:6)){
  j=i*3+1
  corr_mat <- cor(t(datExpr1[j:(j+2),]), t(datExpr1[j:(j+2),])) # Correlation matrix
  corrplot(corr_mat, type = "upper", method = "number", tl.pos = "n") # Correlation plot
}

```

### 1.2.3 Sample clustering

In order to visualize better the relationship among samples compute the distance among them (e.g. euclidean), perform a hierarchical clustering, and plot the results as a tree.

```{r hc samples}

sampleTree = hclust(dist(datExpr1), method = "average")
plot(sampleTree, xlab = "Samples")

```

Does the result change with different distance metrics or clustering methods?

## 1.3 Time series filter

A time series usually inquires something changing over time. Therefore time should account for the greatest amount of variability in the data.

### 1.3.1 Sample PCA

Compute the principal components, plot the variance explained by the principal components and the 2D combination of the first 3 ones.

```{r, PCA samples}

expr.pca <- prcomp((datExpr1), scale. = TRUE)

pcr_var <- expr.pca$sdev**2
plot(pcr_var/sum(pcr_var)*100, type = 'b', ylab = "Variance explained [%]", xlab = "Principa Component")

ggplot(data = as.data.frame(expr.pca$x), aes(x=PC1, y=PC2)) + geom_text(aes(col=rownames(datExpr1), label=rownames(datExpr1))) + ggtitle("Principal component analysis - PC1 vs PC2")
ggplot(data = as.data.frame(expr.pca$x), aes(x=PC1, y=PC3)) + geom_text(aes(col=rownames(datExpr1), label=rownames(datExpr1))) + ggtitle("Principal component analysis - PC1 vs PC3")
ggplot(data = as.data.frame(expr.pca$x), aes(x=PC2, y=PC3)) + geom_text(aes(col=rownames(datExpr1), label=rownames(datExpr1))) + ggtitle("Principal component analysis - PC2 vs PC3")

```

How much variance is explained from the first 3 components? What can you infer from the PC vs PC plots? Given the information here and in 1.2.3 is there any outlier in the samples? If you think so, remove it (them) for the downstream analysis.

```{r, outlier removal}

datExpr2 <- datExpr1[rownames(datExpr1)!="t7C",]
rownames(datExpr2)
dim(datExpr2)

```

# 2 Network analysis

Let's move the focus from the samples to the genes. 

### 2.0.1 Density filter

We are interested in genes that reach a good level of expression. Print the density of the expression matrix obteined so far, chose a threshold and filter out all the genes that don't reach it even once.

```{r, density}

plot(density(as.matrix(datExpr2)), main = "Density plot", xlab = "Sample")
abline(v=2, col='red') # TPM 4

datExpr3 <- datExpr2[,colMaxs(as.matrix(datExpr2))>3]

```

### 2.0.2 Heatmap visualization

Produce a heatmap for the data and the normalized data. Is there any pattern emerging? What does it mean when two scaled profiles are similar?

```{r, heatmap}

datExpr3.scal <- scale(datExpr3, center = TRUE, scale = TRUE)

heatmap(as.matrix(datExpr3), Rowv = NA, Colv = TRUE, scale = "none", main = ("Gene expression heatmap"))
heatmap(as.matrix(datExpr3.scal), Rowv = NA, Colv = TRUE, scale = "none", main = ("Gene expression change heatmap"))

```

## 2.1 Network inference

Enable the multithreading for the library "WGCNA", compute the scale-free topology fit index and the mean connectivity for an array of condidate powers and plot the results. Take a power from the biginning of the plateau of the scale-free topolgy index. Ideally it should be above 0.90, but it is not always possible to reach such a high score.

```{r, network inference}

enableWGCNAThreads()

powers = (1:20)

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr3, corFnc = 'cor', corOptions = list(use = 'p', method = "pearson") , powerVector = powers, verbose = 5)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower = 6

```

Which power fits the data? Is the network scale-free? Why is scale-freedom considered important (in biology)?

## 2.2 Module inference

### 2.2.1 Adjacency, TOM, dissimilarity

Compute the adjacency matrix using the selected power from 2.1, then compute the Topological Overlap Matrix (TOM) and make a dissimilarity matrix out of it (1-TOM). 


```{r, adjacency, results='hide'}

adjacency = adjacency(datExpr3, power = softPower, type = "signed")

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

```

### 2.2.2 Dissimilarity clustering

Cluster the genes using the dissimilaity matrix using the hierarchical clustering and plot it as a tree. Chose a minimum number of genes per module, then perform a clustering through the "cutreeDynamic" function. Plot the tree and the modules together.

```{r, module inference, results='hide'}

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 10;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")

```

Which minimum number of genes seems to be the best one? How many modules have you retrieved?

### 2.2.3 Eigengene clustering

Compute the eigengenes for the modules, correlate dem and compite the dissimilarity (1-cor). Then compute the hierarchical clustering, plot the results. Choose a threshold and merge the modules according to it.

```{r eigengene clustering}

# Calculate eigengenes
MEList = moduleEigengenes(datExpr3, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")

MEDissThres = 0.2
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr3, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

```

What threshold fits the best? How many clusters have you retrieved with this process?

Plot the results for the new cut.

```{r, merged clusters}

sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

```

## 2.3 Module analysis

## 2.3.1 Epression plots

For each module print the most informative representation of the expression, such as:
- eigengene
- heatmap (expression and standardized change)
- matplot (expression and standardized change)
- boxplot for the time points (expression and standardized change)

```{r, single module analysis}

arrk <- sort(unique(mergedColors))

for (i in (1:length(arrk))){
  
    smallmat <- datExpr3[,mergedColors==arrk[i]]
    smallmat.scal <- datExpr3.scal[,mergedColors==arrk[i]]
  
  somePDFPath = paste("auto_report/", toString(arrk[i]), ".pdf", sep="")
  pdf(file=somePDFPath, width = 10, height = 6)
  
  barplot(MEList$eigengenes[,i], xlab = "Sample", ylab = "Eigengene value")
  heatmap.2(as.matrix(smallmat), colsep = FALSE, Rowv = NA, scale = "none",
            dendrogram = 'none', trace = 'none', lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(0.25, 4, 0.25),
            lwid=c(0.25, 4), key=FALSE, keysize=1.0, symkey=FALSE, density.info='none', labCol = NA,
            xlab = "Gene", ylab = "Sample")
  heatmap.2(as.matrix(smallmat.scal), colsep = FALSE, Rowv = NA, scale = "none",
            dendrogram = 'none', trace = 'none', lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(0.25, 4, 0.25),
            lwid=c(0.25, 4), key=FALSE, keysize=1.0, symkey=FALSE, density.info='none', labCol = NA,
            xlab = "Gene", ylab = "Sample")
  
  matplot(as.matrix(smallmat),type="l", xlab="Sample",ylab="log2 expression")
  matplot(as.matrix(smallmat.scal),type="l", xlab="Sample",ylab="standardized change")
  
  time_replicates <- as.factor(rep(c(c(2,2), rep(3,3), rep(4,3), rep(5,3) ,rep(6,3), c(7,7), rep(8,3)), length(colnames(smallmat))))
  
  small.plot <- as.data.frame(smallmat) %>%
    rownames_to_column(var="sample") %>%
    gather(ORF, value, -sample) %>%
    separate(sample, into=c("time", "replicate"), 2) %>%
    ggplot(aes(x=time, y=value)) +
    geom_boxplot() +
    labs(y="Unscaled expression", x="Time point")
  plot(small.plot)
  
  smallscaled.plot <- as.data.frame(smallmat.scal) %>%
    rownames_to_column(var="sample") %>%
    gather(ORF, value, -sample) %>%
    separate(sample, into=c("time", "replicate"), 2) %>%
    ggplot(aes(x=time, y=value)) +
    geom_boxplot() +
    labs(y="Scaled expression", x="Time point")
  plot(smallscaled.plot)
  
  dev.off()
}

```


## 2.3.1.1 Module-module relationships

In order to understand the relationships among modules it is possible to visualize the scatterplots of their eigengenes and compute again the correlation (as in 2.2.3). What can this visualization say about the interplay of the modules?

```{r, global topological features analysis}

datME=moduleEigengenes(datExpr3,mergedColors)$eigengenes
sizeGrWindow(8,9)
plotMEpairs(datME)

```

## 2.3.2 Topological features analysis

### 2.3.2.1 Global features

Compute some topological features (degree) that can be meaninful to answer the biological questions about the community and formulate new hypothesis. Since the network we are analysing is a weigthed one, does it affect the computation of such features?

The function "intramodularConnectivity" from WGCNA computes both the connectivity (aka the degree) inside each module and the total one.

```{r, global topological features analysis 1}

Alldegrees1 <- intramodularConnectivity(adjacency, mergedColors)
head(Alldegrees1)

```

Compare the first ranking genes for the total degree and the degree computed within each module. Do you obtain the same results?

```{r, global topological features analysis 2}

head(Alldegrees1[order(Alldegrees1$kTotal, decreasing = TRUE),])
head(Alldegrees1[order(Alldegrees1$kWithin, decreasing = TRUE),])

```

In order to have a more clear perception of the difference, make a scatter plot to compare the total and the within measures.

```{r, global topological features analysis 3}

plot(Alldegrees1$kTotal, Alldegrees1$kWithin, col=(mergedColors), main="Within and Total Degree comparison", xlab="Total Degree", ylab="Within Degree")

```


### 2.3.2.2 Local features

Compute the Module Membership (correlation between a gene and its module eigengene) for all the genes.

```{r, local topological features analysis 1}

datKME=signedKME(datExpr3, datME, outputColumnName="MM.")
head(datKME)

```

Compare the Within degree with the Module Membership for each module. Do those metrics give different results?

```{r, local topological features analysis 2}

for (i in (1:length(arrk))){

  which.color=arrk[i];
  restrictGenes=mergedColors==which.color
  verboseScatterplot(Alldegrees1$kWithin[ restrictGenes],
  (datKME[restrictGenes, paste("MM.", which.color, sep="")])^6,
  col=which.color,
  xlab="Intramodular degree",
  ylab="Module membership")

}

```

Retrieve the five highest ranking genes according at least to one of the two intramodular metrics.

```{r, local topological features analysis 3}

for (i in (1:dim(datKME)[2])){
  print(colnames(datKME)[i])
  print(rownames(datKME[order(datKME[,i], decreasing = TRUE),])[1:5])
}
  
```

## 2.3.3 Taxonomic enrichment

Load the taxonomic information and perform an enrichment analysis on the modules. Can you see some species specific modules?

```{r, taxonomy enrichment}

#tax_file <- read.table('genes_tax.txt', header = FALSE, row.names = 1, sep = '\t') # Change it with your path and file
bins <- table(bin.df$bin)

for (i in (1:length(arrk))){

  which.color=arrk[i];
  sub_genes <- colnames(datExpr3)[mergedColors==which.color]
  sub_tax <- (bin.df %>% filter(ORF%in%sub_genes) %>% dplyr::select(bin))$bin
  #sub_tax[is.na(sub_tax)]="other"
  tab_tax <- table(sub_tax)
  mod_counts <- tab_tax[names(bins)]
  names(mod_counts) <- names(bins)
  mod_counts[is.na(mod_counts)]=0
  enrich_tax <- dMWNCHypergeo(mod_counts, bins, sum(mod_counts), rep(1,length(bins)), precision = 1E-7)
  #print(list(color=which.color, probability=enrich_tax, bin_counts=mod_counts))
  
  
  scores <- rbind(mod_counts, mod_counts)
  rownames(scores) <- c("depleted", "enriched")
  
  for (j in (names(bins))){
    
    x <- mod_counts[names(mod_counts)==j]
    X <- bins[names(bins)==j]
    Y <- sum(bins[names(bins)!=j])
    tot <- sum(mod_counts)
    scores[1,j]<-phyper(x, X, Y, tot, lower.tail = TRUE)
    scores[2,j]<-phyper(x, X, Y, tot, lower.tail = FALSE)
    
  }
  print(which.color)
  print(scores)
}

```

## 2.3.4 Functional enrichment

Load the GO term annotation and make the dataframe according to the AnnotationDbi criteria (structure: dataframe, 1st column=GO, 2nd column=evidence, 3rd column=gene name). Then build the GeneSetCollection object using GSEABase package.

```{r, GO term loading}

go_file <- ann.df %>% mutate(ORF=Seqid) %>% dplyr::select(ORF, GO.terms) %>% filter(!is.na(GO.terms)) %>% separate_rows(GO.terms, convert=T, sep=",")
nd_rep <- rep("ND", nrow(go_file))
frameData <- data.frame(cbind(go_file$GO.terms, nd_rep, go_file$ORF))
head(frameData)

frame=GOFrame(frameData,organism="Metagenome")
allFrame=GOAllFrame(frame)

gsc <- GeneSetCollection(allFrame, setType = GOCollection())

```

For each module compute the enrichment analysis using the hypergeometric distribution. Is there any functional characterization of the modules? GO term can be extremely vague, but try to put a label on the modules. You could also inspect the results for a different ontology ("BP", "CC", "MF").

```{r, GO enrichment}

for (i in (1:length(arrk))){

  which.color=arrk[i];
  sub_genes <- colnames(datExpr3)[mergedColors==which.color]
  params <- GSEAGOHyperGParams(name="GSEA Params", geneSetCollection=gsc, geneIds = sub_genes, universeGeneIds = colnames(datExpr3), ontology = "BP", pvalueCutoff = 0.05, conditional = FALSE, testDirection = "over")
  Over <- hyperGTest(params)
  t <- summary(Over)
  print(t)

}

```

# 3 Network visuaization

If you managed to reach this point: congratulation, you can install Cytoscape on your machine! It's pretty strightforward from the site (http://www.cytoscape.org/). So far we used a network, an object with nodes and edges, without seeing its classical graphical representation. Therefore export the graph in 2 files that can be used imported in Cytoscape.

```{r, network visualization}

cyt = exportNetworkToCytoscape(TOM, edgeFile = "CytoscapeInput-edges.txt", nodeFile = "CytoscapeInput-nodes.txt", weighted = TRUE,
threshold = 0.2, nodeNames = colnames(datExpr3), nodeAttr = mergedColors)

```

Open Cytoscape and close the central window, then load the edges file using the third icon on the top left of the screen (a blue three-nodes network with an orange arrow). After asking to inport the edges file Cytoscape will ask you to select source and target. Click on the small black arrow on the header of the first column, then on the source symbol (green plain dot), do the same on the second but select the target icon (red concentric dots). Leave all the other fields as they are and import the network. Probably a message will ask you if you want to change to big network layout: say yes. Now import the nodes as a table using the fifth icon (blue grid and orange arrow), asking to merge the information using the first column (the merging column has a small key picture on it), but Cytoscspe should try to match the first columns, which is correct.

Now compute some topological features on the network using the buttons "Tools" -> "NetworkAnalyzer" -> "Network Analysis" -> "Analyze Network" and wait for the results.

Open the "Style" and "Node" panel on the left column. Now change the color clicking the right box in the "Fill Color" filed, select the column that stores the modules' colors using the box "Column" and select "Discrete Mapping" for the "Mappyng Type" box. Assign different color for the different modules.

Change the general layout of the graph clicking the button "Layout" from the top of the screen, then click "Perfuse Force Directed Layout" and select "weight" to encode the edges' weight.


