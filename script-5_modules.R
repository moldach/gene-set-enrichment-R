load("kaforou.rda")

## When using read.csv() make sure to use stringsAsfactor=FALSE
## to assure you the character strings are not 
## google global .Rprofile to change this so you will never accidently forget
## to turn character strings into factors

## some sanitizing
kaforou$genes <- kaforou$genes[ , c("ID", "Symbol", "Definition", "Entrez_Gene_ID", "RefSeq_ID") ]
kaforou$genes$Symbol <- as.character(kaforou$genes$Symbol)
kaforou$genes$ID <- as.character(kaforou$genes$ID)
kaforou$genes$Entrez_Gene_ID <- as.character(kaforou$genes$Entrez_Gene_ID) # R converted characters into factors with load() so lets turn them back into characters
kaforou$targets$cg <- gsub("kaforou\\.", "", kaforou$targets$cg)

## basic differential analysis - just to look what is there
d <- model.matrix(~ 0 + cg, data=kaforou$targets)
colnames(d) <- gsub("cg", "", colnames(d))
c <- makeContrasts(Malawi="(Malawi.TB-Malawi.LTB)", SA="(SA.TB-SA.LTB)", 
  TBvsLTBI="((Malawi.TB-Malawi.LTB)+(SA.TB-SA.LTB))/2", levels=d)
fit <- eBayes(contrasts.fit(lmFit(kaforou, d), c))
topTable(fit, coef="TBvsLTBI")

library(tmod)
res <- tmodLimmaTest(fit, fit$genes$Symbol)
tmodPanelPlot(res, filter.rows.pval = 1e-5)
# We see type 1 interferon response which we expect, B cells as well


## first, a small example
## we manually inspect correlations for some of the genes
i <- "ILMN_1799848"
#Ctrl + Shift + 2
kaforou$genes[i,]

## absolute correlation coefficients 
cc <- abs(cor(kaforou$E[i, ], t(kaforou$E))) # Making a correlation between this one gene and all the other genes
dim(cc) # 1 row, many columns
cc <- t(cc)[,1]  # Transpose this very long row into one column
head(sort(cc, decreasing=TRUE))  # Sort from largest to smallest

## what are these genes?
kaforou$genes[ names(head(sort(cc, decreasing=T))), ]

## plot showing correlation between two genes
j <- "ILMN_2261600"
plot(kaforou$E[i, ], kaforou$E[j,], pch=19, col=factor(kaforou$targets$group), 
     bty="n", xlab=kaforou$genes[i, "Symbol"], ylab=kaforou$genes[j, "Symbol"])
legend("topleft", levels(factor(kaforou$targets$group)), pch=19, col=1:2, bty="n")

# plot(ks$E[i, ], ks$E[j,], pch=19, col=factor(ks$targets$group), bty="n", xlab=ks$genes[i, "Symbol"], ylab=ks$genes[j, "Symbol"])
# legend("topleft", levels(factor(ks$targets$group)), pch=19, col=1:2, bty="n")


## 1. primitive clustering

## A rudimentary clustering can be directly produced when creating a heatmap.
## for a heatmap, we want to select genes that are significantly regulated

tt <- topTable(fit, coef="TBvsLTBI", number=Inf, sort.by="p")
X <- kaforou$E[tt$ID[1:150], ]
head(X[,1:5])

## select a color scheme
## the fu is a *function* that can be used to generate colors. 
## for example, try fu(5) and fu(15)
## also note how R can understand colors: using the HTML hexadecimal notation style, "#320A4B" 
library(RColorBrewer)
fu <- colorRampPalette(c("purple", "black", "cyan"))
fu(10) # If you use an uneven number you get black in the middle

## heatmap
## trace: no idea who uses this option, ever
## scale: the expression of a gene will be converted to a z-score.
##        alternative: calculate z-scores for a column
## dendrogram: row, col or both [default]. 
library(gplots)
group.col <- c("cyan", "purple")[ factor(kaforou$targets$group)]
colnames(X) <- paste0(kaforou$targets$group, ".", colnames(X))
par(bg="#eeeeee")
heatmap.2(X, trace="n", scale="r", col=fu, 
          labCol=kaforou$targets$group,
          ColSideColors = group.col)
# We see how many modules we would be able to identify (vertically) 7 bars


## a second version of the same figure. Here, we do not want to reorder the samples and
## calculate a dendrogram for the columns, but simply show the samples grouped
ordo <- order(kaforou$targets$group)
heatmap.2(X[,ordo], trace="n", scale="r", dendrogram="row", col=fu, Colv=NULL,
          labCol=kaforou$targets$group[ordo],
          ColSideColors = group.col[ordo])

# Correlation clusters take ^2 more compute resources than genes (e.g. 40,000genes * 40,000)
plotDensities(kaforou,legend=F)
axis(1)
axis(side=1,at=1:15)

## 2. Simple clustering. 

## First, select a number of genes based on their absolute expression
## level. We select only genes for which the upper boundary of IQR
## is higher than 7 (arbitrary threshold)
 
uq <- apply(kaforou$E, 1, function(x) quantile(x, 0.75)) 
sel <- uq > 7 # Select upper quartile greater than 7
ks <- kaforou[sel, ] 

## from what is left, we select 20% of the genes with the largest IQR
iqrs <- apply(ks$E, 1, IQR)
sel  <- iqrs > quantile(iqrs, 0.8) 
sum(sel)
ks <- ks[sel, ]
X <- ks$E

## Clustering with hclust
## We start with clustering the samples
## First, you need to calculate the distances. 
dd <- dist(t(X))
hc <- hclust(dd)
plot(hc, labels=ks$targets$group)
# what if we choose three clusters
# abline(h=101)
clusts <- cutree(hc, h = 101)

clusts <- cutree(hc, h = 180)
# we started with clustering the samples and not the genes because we have fewer and it will be easier to see

## clusts is a vector; names are the sample IDs, 
## we can ask how these clusters relate to our known phenotype
table(clusts, ks$targets$group)  # We see that first cluster is almost all TB, second and third clusters are predominantly LTB
table(clusts, ks$targets$Cohort) # We see that clusters also seem to correspond to the different populations; couldn't decide what is more important (TB or cohort) we didn't include it but it would likely cluster by sex as well

## We can do a similar thing for genes
## first, with euclidean distances
dd <- dist(X)
hc <- hclust(dd)
plot(hc, labels=FALSE)
# Much harder to see if there is anything meaninful here
# Decision to cut tree is rather arbitrary then
abline(h=40)
clusts <- cutree(hc, h = 40)
uc <- unique(clusts)
uc # number from 1:11
head(clusts) # which ones associate to clusters
m2g <- sapply(uc, function(i) unique(ks$genes[ names(clusts[ clusts == i ]), "Symbol" ]), simplify=FALSE)
m <- data.frame(ID=paste0("M.", uc), Title=uc, N=sapply(uc, function(i) sum(clusts == i)))
names(m2g) <- m$ID

library(tmod) # Let's do a hypergeometric to test against the background of all genes in this reduced dataset (not the full dataset)
res <- sapply(m2g, function(cl) tmodHGtest(cl, ks$genes$Symbol), simplify=F) # 
res # we see nulls where we couldn't find any genes, we also see there aren't that many genes in each of the modules
## Go back up the script and instead of choosing 40 choose closer to the bottom of the tree and then
## sort(sapply(m2g,length))
## res shows N=gene universe, n=genes in modules, B=number of ?, b=number of genes in this module
res <- res[ !sapply(res, is.null) ]
res <- res[ sapply(res, nrow) > 0 ]
res <- sapply(res, function(x) { x$AUC <- log(x$E) ; x }, simplify=F) # log 
names(res) <- paste0("M.", seq_along(res)) # Create module names
tmodPanelPlot(res)
tmodPanelPlot(res, pval.thr = 0.05)

## now, we will use correlation as distance metrics
cc <- cor(t(ks$E))
cc <- 1 - abs(cc)
cd <- as.dist(cc)

hc <- hclust(cd)
plot(hc, labels=ks$targets$group)
clusts <- cutree(hc, h = 0.5)

uc <- unique(clusts)
m2g <- sapply(uc, function(i) ks$genes[ names(clusts[ clusts == i ]), "Symbol" ], simplify=FALSE)
m <- data.frame(ID=paste0("M.", uc), Title=uc, N=sapply(uc, function(i) sum(clusts == i)))
names(m2g) <- m$ID

sel <- m$N > 30
sum(sel)
m2g <- m2g[sel]
m <- m[sel,]

res <- sapply(m2g, function(cl) tmodHGtest(cl, ks$genes$Symbol), simplify=F)
res <- res[ !sapply(res, is.null) ]
res <- res[ sapply(res, nrow) > 0 ]
res <- sapply(res, function(x) { x$AUC <- log(x$E) ; x }, simplify=F)
tmodPanelPlot(res)


library(WGCNA)
## uncomment the following in terminal R, but not in Rstudio
# allowWGCNAThreads()

## we test a range of the "power" parameter
powers <- c(c(1:5), seq(from = 6, to=30, by=3))
s <- pickSoftThreshold(t(ks$E), powerVector=powers, verbose=5)

with(s$fitIndices, plot(Power, -sign(slope) * SFT.R.sq, type="n",
  xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit,signed R^2",
  main = "Scale independence"))

with(s$fitIndices, text(Power, -sign(slope) * SFT.R.sq, Power, col="red"))
abline(h=0.90,col="red")

s <- 16
adjacency <- adjacency(t(ks$E), power = s)
TOM <- TOMsimilarity(adjacency)
d.t <- 1 - TOM
geneTree <- hclust(as.dist(d.t), method = "average")

plot(geneTree, labels=FALSE, hang=0.4)
m <- cutreeDynamic(dendro=geneTree, distM=dissTOM, deepSplit=2, minClusterSize=30, pamRespectsDendro=FALSE)
col <- labels2colors(m)
plotDendroAndColors(geneTree, col)
length(m)

library(tmod)
head(ks$genes)
 
res <- sapply(unique(names(m)), function(i) tmodHGtest(fg=ks$genes$Symbol[m == i], bg=ks$genes$Symbol), simplify=F)
names(res) <- paste0("M", names(res))
res <- res[ sapply(res, nrow) > 0 ]
sapply(res, function(x) { colnames(x)[7] <- "AUC" ; x }, simplify=F)
res <- sapply(res, function(x) { colnames(x)[7] <- "AUC" ; x }, simplify=F)
tmodPanelPlot(res)
