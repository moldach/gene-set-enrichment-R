## load limma and create an EList
library(limma)
# remember: we created the object E using the following command line:
#E <- new("EListRaw", list(E=E.mat, genes=genes, targets=ph))

## are the arrays similar?
plotDensities(E, log=TRUE, legend=FALSE)
# Arrays with only background noise are worrysome

## the data we used here does not need to be background corrected first but usually you would do this normally
## background correct first 
# E.bg <- backgroundCorrect(E, method="normexp")
# plotDensities(E.bg, log=TRUE, legend=FALSE)

# But we do need to normalize between the arrays
En <- normalizeBetweenArrays(E, method="quantile")
plotDensities(En, log=TRUE, legend=FALSE)
# Shows all the denisities between arrays are on the same 
# Unlike RNA-seq microarray is an imaging (kind of like correcting for library size in RNA-seq)

## remove unnecessary objects
# rm(E, E.bg, geo)

## average per probe; some probes are spotted multiple times
## not required for RNASeq, of course
dim(En)
## how many unique probes?
length(unique(En$genes$NAME))
## You see there are a number of probes that are repeated

## in general, you need platform identifiers (platform specific probe ID), such as A_xx_xxxx for Agilent
## or ILMN_XXXX for Illumina
En <- avereps(En, ID=En$genes$NAME) # We are averaging for the duplicate of the probes
dim(En) # or: nrow(En)

## Remove the probes that are controls
En <- En[ En$genes$CONTROL_TYPE == "FALSE", ]
dim(En)
# Multiple different probes for each gene

# Affymatrix data requires a different pipeline so don't use this script for it

## do we see any apparent problems with PCA?
## January likes to do a PCA on every dataset he gets first
pca <- prcomp(t(En$E), scale.=TRUE)  # t is for transposing the data
plot(pca$x[,1], pca$x[,2], pch=19)
# use pca3d to get a 3D view of the PCA
library(pca3d)
pca3d(pca)

# alternatively: PCA for genes, not very interesting
pca <- prcomp(En$E, scale.=TRUE)
plot(pca$x[,1], pca$x[,2], pch=19)
pca3d(pca)

## lets go back to regular PCA
pca <- prcomp(t(En$E), scale.=TRUE)
plot(pca$x[,1], pca$x[,2], pch=19, col=factor(En$targets$group))
legend("topright", as.character(levels(factor(En$targets$group))), pch=19, col=1:3)
plot(pca$x[,3], pca$x[,4], pch=19, col=factor(En$targets$group))
legend("topright", as.character(levels(factor(En$targets$group))), pch=19, col=1:3)
# NID is non-infected disease group
# LTBI is 
# TB is tuberculoses
# You see the third component really distinguishes between TB and and the other two groups (non infected disease group)


## which genes matter for that component?
ord <- order(abs(pca$rotation[,3]), decreasing = T)
head(En$genes[ord,])

## if you wish to filter your data by variance, here is how
## we use the non-parametric interquartile range to do the filtering
iqrs <- apply(En$E, 1, IQR)
cutoff <- quantile(iqrs, 0.90) # number that is smaller than 25% largest IQRs
sum(iqrs < cutoff) / nrow(En) * 100 # 90%
En.f <- En[ iqrs > cutoff, ] # retain only 25% probes with the largest IQR

## create a basic design. We first test for differences between the main experimental groups.
## there are many ways to do this.
En$targets$group <- factor(En$targets$group, levels=c("NID", "LTBI", "TB"))
d <- model.matrix(~ group, data=En$targets)
# View(d) 
# When we see this we notice the first group NID is not here but the intercept (grand average) and the second and third group is difference of LTBI and first group and the next is TB to the first group.
fit1 <- lmFit(En, d)
fit1 <- eBayes(fit1)
tt <- topTable(fit1, coef="groupTB")  # After the group here you need to put in one of the levels
#View(tt)
# Minus logFC (TB-NBID) means its higher in NBID
# Plus logFC (TB-NBID) means its higher in Tuberculosis
# adjP.Val is multiple testing with the benjamini-hofberg adjustment

## you can play with topTable parameters. How to show more than 10 lines? 
## you use the parameter number=x to show how many lines you want to show
## How to show all genes with q < 0.0001? Read the manual (?topTable)
## You can set the p.value<0.0001 command


## what is the impact of filtering on the results?
#fit1b <- eBayes(lmFit(En.f, d))

#nrow(topTable(fit1, coef="groupTB", number=Inf, p.value=0.001, lfc=1))
#nrow(topTable(fit1b, coef="groupTB", number=Inf, p.value=0.001, lfc=1))

#tt1 <- topTable(fit1, coef="groupTB", number=Inf, sort.by="none")
#tt1b <- topTable(fit1b, coef="groupTB", number=Inf, sort.by="none")
#tt1 <- tt1[ rownames(tt1b), ] # only genes present in tt1b as well
#plot(tt1$logFC, tt1b$logFC)
#plot(tt1$P.Value, tt1b$P.Value, log="xy", pch=19, col="#33333311")
#abline(0,1, col="red")
## bottom line: the filtering did not have a positive impact!

## way #2. More complex, but also more explicit.
d <- model.matrix(~ 0 + group, data=En$targets)
# View(d)
# Now we can define the contrasts ourselves
colnames(d) <- levels(En$targets$group)
fit2 <- lmFit(En, d)

c <- makeContrasts(TBvsNID="TB-NID", LTBIvsNID="LTBI-NID", TBvsLTBI="TB-LTBI", 
                   TBvsAll="TB-(NID+LTBI)/2", # TB by the average of the two
                   levels=d)
fit2 <- contrasts.fit(fit2, c)
fit2 <- eBayes(fit2) # Calculate bayes factors and p-value
# head(fit1$coefficients) # we see that we have groupTB
# head(fit2$coefficients) # 

## the results are actually the same
cor(fit1$coefficients[,3], fit2$coefficients[,1])
smoothScatter(fit1$coefficients[,3], fit2$coefficients[,1]) # The p-values between fit1 and fit2 should be identical

## or, compare the results of the following topTable call with the previous one
topTable(fit2, coef="TBvsNID")


## one of the most common ways of visualizing the results is a "volcano plot"
volcanoplot(fit2, coef="TBvsNID")
# Things on left and right are the high and low expression wih highest p-values

# we can do the same using plot()
tt <- topTable(fit2, coef="TBvsNID", number=Inf)
with(tt, plot(logFC, -log10(P.Value), pch=19))
# show top 50 genes
with(tt[1:50,], points(logFC, -log10(P.Value), pch=19, col="red"))
with(tt[1:50,], text(logFC, -log10(P.Value), labels=GENE_SYMBOL, col="red"))


## the advantages of the second method are apparent if you want to know whether there is
## an interaction between sex and TB

En$targets$gr.sex <- paste0(En$targets$group, ".", En$targets$sex)
En$targets$gr.sex <- factor(En$targets$gr.sex)

d <- model.matrix(~ 0 + gr.sex, data=En$targets)
colnames(d) <- levels(En$targets$gr.sex)
fit3 <- lmFit(En, d)

## we want to know now whether the differences between TB and NID are different in males
## than in females. Thus, we are asking about the difference of differences == interaction

c <- makeContrasts(int="(TB.Female-NID.Female)-(TB.Male-NID.Male)", levels=d)
fit3 <- eBayes(contrasts.fit(fit3, c))
topTable(fit3, coef=1)

## heatmaps
## first, select top 40 genes to show
tt <- topTable(fit2, coef="TBvsNID", number=Inf)
sel <- rownames(tt)[1:40]

## select data to show on heatmap
x <- data.matrix(En$E[sel,])
library(gplots)
genelabs <- En$genes$GENE_SYMBOL[ match(sel, En$genes$NAME)]
heatmap.2(x, trace="n", labRow = genelabs, scale="r", 
             labCol=En$targets$group)
## if you get error message "invalid graphics sate", 
## click on the red dot with an "x" above the plot

colf <- colorRampPalette(c("purple", "black", "cyan"))
heatmap.2(x, trace="n", labRow = genelabs, scale="r", 
             labCol=En$targets$group, col=colf)

## make a bit more contrast. Don't worry about the warning
heatmap.2(x, trace="n", labRow = genelabs, scale="r", 
             labCol=En$targets$group, col=colf, 
             breaks=seq(-2, 2, length.out=21))

## get rid of what we dont want to have
ord <- order(En$targets$group)
heatmap.2(x[,ord], trace="n", labRow = genelabs, scale="r", 
             labCol=En$targets$group[ord], col=colf, 
             breaks=seq(-2, 2, length.out=21),
             key = F, Colv = NULL, colsep=62)
