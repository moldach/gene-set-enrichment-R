# this script is for generating modules for tmod using custom annotations

#load the table with gene names and modules they belong to
# Column 1 should be gene IDs
# Column 2 should be module list

#firedls in table are tab separated and modules in the last field are ", " (comma space) separated
foo = read.table("NCBI.w.annotation.txt", header = T, sep = '\t', stringsAsFactors = F, fill = T)
head(foo)

#extract the gene to module association
g2m <- strsplit(foo$system, ", ")
head(g2m)
names(g2m) <- foo$gene
head(g2m)

#create the module definition 
m <- unique(unlist(g2m))
m <- data.frame(ID=paste0("M.", 1:length(m)), Title=m, stringsAsFactors = F)
head(m)

#create the module to gene association
m2g <- sapply(m$Title, function(mm) names(g2m)[ sapply(g2m, function(g) mm %in% g)], simplify=FALSE)
head(m2g)

foo[ grep("Proteasome bacterial", foo$system), ]
head(m)

#create the custom module in tmod 
library(tmod)
all(names(m2g) == m$Title)
names(m2g) <- m$ID

mymodules <- makeTmod(modules=m, modules2genes = m2g)
mymodules

#save it to a file
save(mymodules, file="mymodules.rda")
