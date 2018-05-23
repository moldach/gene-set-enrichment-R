# install packages: edgeR, statmod

# Data from Tuch, B.B. et al. (2010). Tumor transcriptome sequencing
# reveals allelic expression imbalances associated with copy number 
# alterations. PLoS ONE 5, e9317
rawdata <- read.csv("rnaseq_example.csv", stringsAsFactors=FALSE)

library(edgeR)
library(statmod)
y <- DGEList(counts=rawdata[,4:9], genes=rawdata[,1:3])  # Change the counts to where the columns for counts occur
# column one here is ID, 2 is the gene name and column 3 is the number of exons in the gene

## filter the data to keep genes above background
## rule of thumb is roughly 5-10 counts in 20,0000 (0.05 counts per million)
## this doesn't really change anything, just demonstration
keep <- rowSums(cpm(y) > 0.5) > 1
dim(y)
y <- y[keep, ] # This doesn't do anything for this data but will for realy datasets
dim(y)


## install annotation internal R-package from Entrez and different databases
## our authors provided refseq and we need to convert to GO IDs
# source("https://bioconductor.org/biocLite.R")
# biocLite("org.Hs.eg.db")
library(org.Hs.eg.db)
eg2refseq <- toTable(org.Hs.egREFSEQ)
y <- y[ y$genes$idRefSeq %in% eg2refseq$accession, ] # Keep only the genes we can find with refseq identifiers (example pseudogenes or non-coding genes)
y$genes$EntrezID <- eg2refseq$gene_id[ match(y$genes$idRefSeq, eg2refseq$accession) ] # add entrez identifier to genes
# head(y$genes)

# remove gene duplicates (Entrez genes) because we matched to transcripts but we do GSEA on genes and not transcripts - in this case we are dealing with genes and not differential transcript expression
## take the dominating transcript (that with the largest expression)
head(which(duplicated(y$genes$EntrezID)))
ord <- order(rowSums(y$counts), decreasing=TRUE)  # When I uncommented this line the script works
y <- y[ord, ]
y <- y[ !duplicated(y$genes$EntrezID), ]

# recalculate library size and calculate normalization factors
# equivalent to apply(y$counts, 2, sum)
y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y)
plotMDS(y)
# We see normal samples group to right and tumors group to the left
# MDS is somewhat like PCA

## also, a PCA. note that we need (1) log counts, (2) remove 
## counts with zero variance

x <- y$counts
vars <- apply(y$counts, 1, var)
keep <- vars > 0 # Keep only those with variance above zero
x <- x[keep,]
x <- cpm(x, log=T, prior.count=3)
## with x, heatmaps can be made (see script for day 2)
## or pca
pca <- prcomp(t(x), scale.=TRUE)
plot(pca$x[,1], pca$x[,2], pch=19, col=rep(c(1,2), 3))
## we actually see a separation between tumour (red points)
## and healthy tissue.

## prepare the model design
Patient <- factor(paste0("P.", c(8,8,33,33,51,51)))
Tissue <- factor(c("N","T","N","T","N","T"))
design <- model.matrix(~Patient+Tissue)
# View(design)
## TissueT column is what interests us - compare tissue tumor to normal

## this is what takes the most time
y <- estimateDisp(y, design, robust=TRUE)

## fit the model
fit <- glmQLFit(y, design)

## limma makes eBayes for every coefficient at once
## calculate the results for coefficient of interest, otherwise it will take average of all the coefficients (ones your not interested in)
lrt <- glmQLFTest(fit, coef="TissueT")
topTags(lrt)
# which genes are significant
de <- decideTestsDGE(lrt, p.value=0.01, lfc=1)
table(de)

## install GO.db
# biocLite("GO.db")
library(GO.db)

go <- goana(lrt)
topGO(go, ont="BP", sort="Up", n=30)
# We see a large number of terms that are absolutely not interesting at all
# fairly vague 

## alternatively, we could use the GOrilla online tool
## at http://cbl-gorilla.cs.technion.ac.il/
## for this, we need to prepare a list of genes ordered by p-value
tt <- topTags(lrt, n=Inf)
write.csv(tt$table$nameOfGene, file="output.txt", row.names = FALSE, quote = FALSE)

## now, go to the GOrilla web site and upload "output.txt"
