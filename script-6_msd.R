## First, we calculate a principal component analysis (PCA) transformation
## This will allow us to visualize the highly multivariate data on a two- 
## or three-dimensional plot

pca <- prcomp(t(En$E), scale.=TRUE)

## using pca3d, we can visualize the PC's. Each point on the graph
## corresponds to one sample, and we can color these samples by the group 
## that they belong to.

## First, we use a group in which we replace "NID" as well as "LTBI" with "CTRL",
## such that we only have two categories. NID and LTBI are so similar, it will be less confusing
## to put them together.
group <- En$targets$group
group <- as.character(En$targets$group)
group[ group != "TB" ] <- "CTRL"
group <- factor(group)
pca3d(pca, group=group)

## if you cannot install pca3d, use the pcaplot function from the tmod
## package
# pca2d <- tmod::pcaplot


## If you inspect the above 3d graphics, you will notice that the main differences
## between the groups appears on component 3. We can show them on a 2D plot
## (in the regular plot window)
pca2d(pca, components=c(2,3), 
      group=group, legend="topleft", radius=2)

## we can add confidence ellipses. 95% will look very broad, we will only add 75%
pca2d(pca, components=c(2,3), group=group, legend="topleft", 
      ellipse.ci = 0.75, show.ellipses = TRUE, radius=2)


## There are several components. How do we know which of these components is related
## to the differences between the groups?
## The object  pca  is a list of matrices (and some other data types) that contain
## useful information. One of these matrices, pca$x, contains the coordinates of
## all samples in the new coordinate system -- that is, it contains the principal 
## components. We can use it to ask which of these is truly correlated with the
## group differences.

## For this, we fit a series of linear models (but not with limma!), in which the group is
## the predictor, and the components (each separately) are the response variables.

foo <- summary(lm(pca$x ~ group))
rs <- sapply(foo, function(x) x$r.squared)
plot(rs[1:15], type="b", bty="n", ylab="rÂ²")
## the plot shows a clear increase at the third component -- this is it!
## what about sex?

foo <- summary(lm(pca$x ~ En$targets$sex))
rs2 <- sapply(foo, function(x) x$r.squared)
lines(rs2[1:100], type="b", col="red")
legend("topright", c("group", "sex"), lty=1, pch=1, col=c("black", "red"), bty="n")
## clearly, the component 5 is interesting. We can update the pca plot accordingly

## since this is superuseful, we will define a function that quickly does exactly that
plotR2 <- function(x, group, n=min(25, ncol(x))) {
  foo <- summary(lm(x ~ group))
  rs2 <- sapply(foo, function(x) x$r.squared)
  plot(rs2[1:n], type="b", bty="n")
  abline(v=1:n, col="#cccccccc")
}

plotR2(pca$x, En$targets$sex)

pca2d(pca, components=c(3,5), group=En$targets$sex, legend="topleft", 
      ellipse.ci = 0.75, show.ellipses = TRUE, radius=2)

## OK, but what about the genes? Let's start with a biplot -- apart from samples,
## there will be arrows indicating the direction of the most important genes for each direction
pca2d(pca, components=c(3,5), group=group, legend="topleft", 
      ellipse.ci = 0.75, show.ellipses = TRUE, biplot=TRUE)
## ok, not very useful, since the probe identifiers are not self-explanatory.
## However, we can manually take a look at the genes with the greates weight in 
## each component. Specifically, we are interested in component 3

## we order the genes by the decreasing absolute weight in component 3
ord <- order(abs(pca$rotation[,3]), decreasing=TRUE)
## we will recognize some genes, e.g. CEACAM1
head(En$genes$GENE_SYMBOL[ord], 30)

## Clearly, this is now an ordered list of gene names.
## You know what we do to ordered lists of gene names.

res <- tmodCERNOtest(En$genes$GENE_SYMBOL[ord])
head(res)


## Actually, tmod has a facility to run PCA and gene set enrichment on the 
## components automagically. 
tmodPCA(pca, genes=En$genes$GENE_SYMBOL, components=c(3,5), 
        plot.params = list(radius=2, legend="topright", group=group)  )

## Unfortunately, we do not see any enrichment on the vertical axis which corresponds 
## to males / females. However, we may try another thing. Instead of ordering by the absolute 
## value of the loading, it is possible to sort the list simply by the weight and then 
## enrich twice: from top to bottom and from bottom to the top. This is achieved
## with the mode="c" parameter:
res <- tmodPCA(pca, genes=En$genes$GENE_SYMBOL, components=c(3,5), mode="c",
        plot.params = list(legend="topright", group=group, radius=2)  )


## What else can we do with a PCA? First of all, we can recalculate it for a group of
## points. PCA is a simple matrix calculation, and we have all the necessary elements.
## Here is how the principal components are calculated from the expression matrix.
## First, the data is converted to z-scores (scaled).

x <- En$E
x <- (x - pca$center)/pca$scale

## we transpose the matrix, such that the samples are in rows and genes in columns.
x <- t(x)

## Then, the scores are calculated by calculating the cross product between the 
## normalized expression values and the rotation matrix (matrix with loadings). The
## result are principal components. You will find that the result is identical to the
## pca$x
pc <- x %*% pca$rotation
identical(pca$x, pc)
par(mfrow=c(1,2))
pca2d(pca$x)
pca2d(pc)

## Now, imagine that we have two data sets, training and test. We calculate PCA for the
## first data set, but then we can apply the calculated scores to the second, training data 
## set.

n <- ncol(En)
sel <- sample(1:n, n*0.7)
length(sel)
train <- En[,sel]
test  <- En[, setdiff(1:n, sel)]

pca2 <- prcomp(t(train$E), scale.=TRUE)
plotR2(pca2$x, group[sel])
pca2d(pca2, components=c(1,3), group=group[sel], radius=2)

test.pc <- t( (test$E - pca2$center)/pca2$scale) %*% pca2$rotation
test.sel <- setdiff(1:n, sel)
points(test.pc[,1], test.pc[,3], pch=19, col=group[test.sel])

## we now show the test set points on the plot and demonstrate that the PCA, although
## calculated based on the training set only, also works for the test set! In fact, we could
## use this property to build a machine learning algorithm (but there are better options)
# text(test.pc[,4], test.pc[,3], test$targets$group, pos=2, cex=2)
text(test.pc[,4], test.pc[,3], test$targets$group, pos=2, cex=1)


## Alternatives to PCA.
## First alternative: ICA. ICA (independent component analysis) has a major advantage; instead
## of producing orthogonal (uncorrelated) components, it aims at producing components that are
## truly independent. This is actually quite a thing; however, ICA (despite the name of the
## package, "fastICA") is very slow and requires substantial amounts of memory. Therefore it is 
## necessary to severly limit the number of variables tested in ICA.

iqrs <- apply(En$E, 1, IQR)
cutoff <- quantile(iqrs, 0.98)
En.f2 <- En[ iqrs > cutoff, ]
library(fastICA)

## this is to ensure that repeating this command gives exactly the same results
set.seed(12345) 
fica <- fastICA(t(En.f2$E), n.comp=9)

## which of the 9 calculated components corresponds to TB and sex?
plotR2(fica$S, En$targets$sex)
plotR2(fica$S, group)

## please use the components that you have identified on the plot!!
grsex <- paste0(group, ".", En.f2$targets$sex)
## separating all four subgroups:
pca2d(fica$S, components=c(7, 8), group=grsex, legend="left", radius=2)

## which genes have the highest weights in components 1 and 9?
weights <- fica$K %*% fica$W
ord <- order(abs(weights[,7]), decreasing=T) # sort by weights in component 1
head(En.f2$genes[ord, ]) # yay!

ord <- order(abs(weights[,8]), decreasing=T) # sort by weights in component 1
head(En.f2$genes[ord, ]) # yay!


## t-SNE is a novel algorithm, an alternative to PCA. it seems to work nice 
## on single cell RNA samples. One of the many R implementations can be found in the package
## tsne.

library(tsne)
## calculate Euclidean distances between the samples
dd <- dist(t(En.f2$E))
ts <- tsne(dd, perplexity=15, max_iter=1000, k=9)
plotR2(ts, En$targets$sex)
pca2d(ts, components=c(2,4), group=grsex, legend="topleft", radius=2)

## the produced components, however, cannot be interpreted.
