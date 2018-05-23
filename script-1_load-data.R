## create a new directory and change the workspace to there
library(here)
here()

## load packages
library(GEOquery)
geo <- getGEO(GEO="GSE28623") # this may take a while
# if this throws up an error from using a locally cached file you need to specify a new folder

## alternatively, using local copy:
# load("gse28623.rda")
# geo <- GSE28623
# rm(GSE28623)

## what is this?
geo
class(geo)
str(geo)
names(geo)

## select the first (and only) member of the geo list
## you can select multiple GEO sets in an object but since we only have we are calling the first and only member
geo <- geo[[1]]

## we now convert the "ExpressionSet" object into 
## an "EList" object from limma
## (if you have two-color arrays, watch out!)
ph <- as(phenoData(geo), "data.frame")
E.mat  <- exprs(geo)
genes <- as(featureData(geo), "data.frame")

## Take a brief look at these objects with View()
# View(ph)

## ph now contains the phenotypic data
## much of that is redundant. We only choose 
## these columns which we need
## When working with your own projects, you 
## will have to figure out which columns
## are interesting
head(ph)
ph <- ph[ , c("geo_accession", "source_name_ch1", "characteristics_ch1.1")]
colnames(ph) <- c("ID", "group", "sex")
ph$group <- gsub(".* ", "", ph$group)  # Remove all that comes before a . and replace with nothing
ph$sex <- gsub(".* ", "", ph$sex) # Remove all that comes before a . and replace with nothing
head(ph)

## E.mat contains the actual expression data
str(E.mat)
head(E.mat[,1:10])

## genes is meta information on the features
## again, we will choose some of the columns
## when working with your own projects, you 
## will have to figure out which columns
## are interesting
## to do this run:
# colnames(genes)
## and select which ones you want to keep
genes <- genes[ , c("NAME", "CONTROL_TYPE", "REFSEQ", "GENE_SYMBOL", "DESCRIPTION")]
head(genes)
tail(genes)

## do the objects fit?
dim(E.mat)
dim(genes)
dim(ph)

## we now load limma and create an EList
library(limma)
E <- new("EListRaw", list(E=E.mat, genes=genes, targets=ph))
str(E)
class(E)
E[1:10,1:2]
E[1:10,1:4]

dim(E)
