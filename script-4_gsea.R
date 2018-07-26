## This example uses the data from script 1 and script 2

## find out the Entrez gene ids, this time using biomaRt
# source("https://bioconductor.org/biocLite.R")
# biocLite("biomaRt")
library(biomaRt)

mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
# Look at the avilable datasets in DB
# foo <- listDatasets(mart)


## this show us what keys we can use to get the data, and what atributes we
## can get from the mart
f <- listFilters(mart)
# View(f)
# If you don't have default annotations you can get the homolgs
# getBM(attributes = c("entrezgene", "description"), filters="hgnc_symbol", values = "IL10", mart=mart)
## or
# d <- getBM(attributes = c("entrezgene", "description", "refseq_mrna"), filters="refseq_mrna", values=E$genes$REFSEQ, mart=mart)

a <- listAttributes(mart)
# View(a)



## we match the existing REFSEQ ids to the mrna refseq
d <- getBM(attributes=c("entrezgene", "description", "refseq_mrna"), filters="refseq_mrna", mart=mart, values=E$genes$REFSEQ)
En$genes$EntrezID <- d$entrezgene[ match(En$genes$REFSEQ, d$refseq_mrna) ]

# At this stage you could remove all those without an Entrez id and then re-run the differential expression analysis - this is not done here but you should

## we repeat the fit to get the EntrezIDs in the fit object
d <- model.matrix(~ 0 + group, data=En$targets)
colnames(d) <- levels(En$targets$group)
c <- makeContrasts(TBvsNID="TB-NID", LTBIvsNID="LTBI-NID", TBvsLTBI="TB-LTBI", levels=d)
fit2 <- lmFit(En, d)
fit2 <- eBayes(contrasts.fit(fit2, c))

## now that we defined fit2 (which contains the interesting coefficents),
## we run basic gene set enrichment analysis using tools at hand 
## GO analysis with topGO from limma
res <- goana(fit2, geneid=En$genes$EntrezID, coef="TBvsNID")
restable <- topGO(res, ontology="BP")
View(restable)

## similar analysis can be done for KEGG pathways
res <- kegga(fit2, geneid=En$genes$EntrezID, coef="TBvsNID")
topKEGG(res)

## using topGO
# source("https://bioconductor.org/biocLite.R")
# biocLite("topGO")
library(topGO)

## first, create a new results table with all genes
## and get rid of these without an EntrezID
tt <- topTable(fit2, coef="TBvsNID", number=Inf, genelist=En$genes)  # Get table with results, extract list of genes with p-value
tt <- tt[ !is.na(tt$EntrezID), ]   # Removes genes with no p-value
tt <- tt[!duplicated(tt$EntrezID), ]  # Remove duplicated entries
rownames(tt) <- tt$EntrezID  


## topGO wants a vector with p-values where the names of the vector
## are Entrez IDs
universe <- tt$adj.P.Val   # contains gene ontology
names(universe) <- tt$EntrezID
data("geneList")

## first topology is calculated
go <- new("topGOdata", ontology="BP", # You can also put in MF and CC if you want but these are often not wanted
  allGenes=universe,
  geneSel=topDiffGenes, nodeSize=10, 
  annotationFun=annFUN.org, # There are many annotations you can choose, as well as your own GO terms (internal IDs to GO terms)
  mapping="org.Hs.eg.db")

## second, we can run the tests
## a couple of these may be of interest
resF    <- runTest(go, algorithm="classic", statistic="fisher" ) # The `classic` algorithm tests each category at a time
resKS   <- runTest(go, algorithm="classic", statistic="ks" )   # 
resF.e  <- runTest(go, algorithm="elim", statistic="fisher" )  # The `elim` algorithm means that once we have tested a low level function we are removing terms from parent node (for example if interferon is significant and immune response is just interferon category just driving the immune response significance? This will remove the interferon category from immune response to see if there is still an immune response)
resKS.e <- runTest(go, algorithm="elim", statistic="ks" )  # The `elim` method is more conservative 
GenTable(go, classicFisher=resF, classicKS=resKS, elimFisher=resF.e,
             elimKS=resKS.e)

## using GOrilla: only online interface
## now, save the genes as a list
write.csv(tt$GENE_SYMBOL, row.names=FALSE, quote=FALSE, file="export.csv")
## go to http://cbl-gorilla.cs.technion.ac.il/ and upload the file
# filter based on 
# ignore anything not orange or red

# Notice that for some GO categories even though there is small p-value there is higher enrichment (this is like logFC vs p-value in differential expression)


## tmod and the
## build-in gene expression modules from tmod

## first, a "manual" analysis to show how tmod works
tt <- topTable(fit2, coef="TBvsNID", number=Inf, sort.by="p") # no MSD sorting yet
gnames <- tt$GENE_SYMBOL
# install.packages("tmod")
library(tmod)
res <- tmodCERNOtest(gnames)
head(res)

# N is the number of genes, 

## show the fragment of the results which corresponds to the genes from a selected module 
## Say we are interested in which genes are in ID=LI.M7.0
# View(res)
# View(tt)
# showModule(tt$REFSEQ, x$GENE_SYMBOL, module = "LI.M7.0") # Extract just the REFSEQ IDs
# foo <- showModule(tt, x$GENE_SYMBOL, module = "LI.M7.0")  # find a subtable of genes in this module
# View(foo)

## Look into the structure of the module set - in thise case,
## the default module set built into tmod
data(tmod)
names(tmod)
# View(tmod$MODULES)
tmod$MODULES2GENES$LI.M7.0


## inspect res now
head(res)

## inspect single modules -- check for evidence of enrichment
## (generate a ROC curve for the first three modules in res)
evidencePlot(tt$GENE_SYMBOL, m = res$ID[1:3], col=2:4, legend = "right")  # res#ID looks at the first three IDs from res
# evidencePlot(gnames, m = res$ID[1:3], col=2:4, legend = "right")  
evidencePlot(gnames, m = "LI.M7.0", col=2:4, legend = "right", gene.labels = TRUE) # res#ID looks at the first three IDs from res
## you see blue module has low p-value just because there is so many genes (n=320), the logFC isnt that much

## CD96 looks interesting
showGene(En$E, En$targets$group) ### Getting error here
tt[ tt$GENE_SYMBOL] ### Getting error here
# boxplot(En$E[ "A_23_P44155", ] ~ En$targets$group) # Not as good as the beeswarm plot in the next line
showGene(En$E[ "A_23_P44155", ], En$targets$group)



## create a tag cloud
library(tagcloud)
tagcloud(res$Title, weights=-log10(res$P.Value), col=smoothPalette(res$AUC))
## size of tags=P value, color=AUC


## tmod has a facility to make running the results on limma coefficients easier
res <- tmodLimmaTest(fit2, fit2$genes$GENE_SYMBOL)
# res <- tmodLimmaTest(fit2, fit2$genes$GENE_SYMBOL, mset=go.bp)
tmodPanelPlot(res, text.cex=0.6, grid="between", filter.rows.pval = 1e-9)

## res is now a list of three tables, one for each contrast

## tmod panel plot
tmodPanelPlot(res, text.cex=0.7, filter.rows.pval = 1e-7)

## but which direction? For that, we need to inspect the results and decide,
## which genes go significantly up, and which go down
pie <- tmodLimmaDecideTests(fit2, genes=fit2$genes$GENE_SYMBOL, 
                            pval.thr = 0.01, lfc.thr = 0.5)
tmodPanelPlot(res, text.cex=1.1, filter.rows.pval = 1e-7, pie=pie)
## a list with an entry for each coefficient 
## a sligthly different representation
tmodPanelPlot(res, text.cex=0.6, filter.rows.pval = 1e-7, pie=pie,
              pie.style = "r", grid = "between")

## using MSigDB database
## For this, downloading the XML MSigDB file is necessary
## this cannot be automatically distributed with tmod
# Go to: http://software.broadinstitute.org/gsea/downloads.jsp 
## Select the most recent data for "All gene sets" under current MSigDB xml file
msig <- tmodImportMSigDB("msigdb_v6.1.xml")
# Dont worry about it complaining about duplicated ones, this is ok

## msig is HUGE. We can, however, test for several different types of
## modules
## E.g. KEGG pathways
kegg <- msig[ msig$MODULES$Subcategory == "CP:REACTOME" ]  # You can replace `REACTOME` with `KEGG`
# kegg <- msig[ msig$MODULES$Subcategory == "CP:KEGG" ]  # You can replace `REACTOME` with `KEGG`

# Every gene set has been derived in a different way, REACTOME was based on chemical ontologies, modules come from co-expression analysis
# Depending on how these were derived there will be a bias

res <- tmodLimmaTest(fit2, fit2$genes$GENE_SYMBOL, mset=kegg)

hallmark <- msig[ msig$MODULES$Category == "H" ]
hallmark$MODULES$Title <- gsub("Hallmark", "", hallmark$MODULES$Title)
res <- tmodLimmaTest(fit2, fit2$genes$GENE_SYMBOL, mset=hallmark)
# C7 looks interesting
res <- tmodLimmaTest(fit2, fit2$genes$GENE_SYMBOL, mset=msig$MODULES$Category == "C7") ### Getting an error here

## what other sets are there?
with(msig$MODULES, table(Category, Subcategory))

go.bp <- msig[ msig$MODULES$Subcategory == "BP" ]
res <- tmodLimmaTest(fit2, fit2$genes$GENE_SYMBOL, mset=go.bp)

## use directly KEGG instead of MSigDB
library(KEGGREST)  # Uses entrez identifiers with an HSA infront HSA:#
pathways <- as.matrix(keggLink("pathway", "hsa"))
pathways <- keggLink("pathway", "hsa")
paths    <- sapply(unique(pathways), function(p) keggGet(p)[[1]]$NAME) # For each path way calls keggGet to get description
paths <- gsub(" - Homo sapiens.*", "", paths)

## here, we construct our own tmod data set
## for non-model organisms you should construct your own data set, 
## this can be done by BLASTing your data set against KEGG, GO, etc.
m2g <- sapply(unique(pathways), function(p) names(pathways)[pathways == p ], simplify=F)
head(m2g)
m <- data.frame(ID=unique(pathways), Title=paths)
kegg2 <- makeTmod(modules=m, modules2genes=m2g)
res <- tmodLimmaTest(fit2, paste0("hsa:", fit2$genes$EntrezID), mset=kegg2)

## we can use KEGGREST to show us the pathway
png <- keggGet("path:hsa05322", option="image")
grid::grid.raster(png)

## we can color the genes in the figure by the log fold change
# source("https://bioconductor.org/biocLite.R")
# biocLite("pathview")
library(pathview)
install.packages("png")
library(png)
genes <- tt$logFC
names(genes) <- tt$EntrezID
foo <- pathview(genes, species="hsa", pathway.id="05322")
image <- readPNG("hsa05322.pathview.png")
grid::grid.raster(image)



## the influence of the sample size on gene set enrichment
sizes <- seq(3, 25, by=3)
tbs <- which(En$targets$group == "TB")
nids <- which(En$targets$group == "NID")
sel <- sapply(sizes, function(s) {
  c(sample(tbs, s), sample(nids, s))
}, simplify=FALSE)
names(sel) <- paste0("S.", sizes) 

 
## for each sample size, run the whole analysis pipeline
res <- sapply(sel, function(s) {
  ee <- En[ , s]
  ee$targets$group <- factor(ee$targets$group)
  d <- model.matrix(~ group, data=ee$targets)
  f <- eBayes(lmFit(ee, d))
  tt.temp <- topTable(f, number=Inf, p.value=0.05)
  print("Number of significant genes:")
  print(nrow(tt.temp))
  res <- tmodLimmaTest(f, f$genes$GENE_SYMBOL, coef=2)
  res <- res[[1]]
}, simplify=FALSE)

tmodPanelPlot(res, filter.rows.pval =1e-7)
