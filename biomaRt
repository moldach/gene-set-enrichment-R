library(biomaRt)
View(listMarts()) # to see what is there

## there are also other marts available at other servers. For example, to 
## access the metazoa ensembl biomart, one can use
# listMarts(host="metazoa.ensembl.org")
# useMart(host="metazoa.ensembl.org")

## connect to biomart database
mart <- useMart("ensembl")

## example for metazoa:
# mart <- useMart("metazoa_mart", host="metazoa.ensembl.org")

## next, what data sets are there? Is there a data set for my organism?
View(listDatasets(mart))

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

## the two following commands show us what we can query (filters) and what we 
## can get (attributes)
View(listAttributes(mart))
View(listFilters(mart))

## for microarrays and other specialized platforms:
## it is possible sometimes to query directly with the platform IDs
## as in the example below - the Illumina microarray probe IDs

## this takes a lot of time and returns a lot of data
bm <- getBM(filters="illumina_humanwg_6_v3", 
            values = kaforou$genes$ID, 
            attributes = c( "hgnc_symbol", "description", "reactome", "illumina_humanwg_6_v3", "go_id", "name_1006", "go_linkage_type"), 
            mart=mart)

## how many probe IDs were recognized?
sum(kaforou$genes$ID %in% bm$illumina_humanwg_6_v3)
all(bm$illumina_humanwg_6_v3 %in% kaforou$genes$ID)

## in the next steps, we will build a tmod module set based on the GO annotations
## from biomaRt. We will use the Illumina probe ID as the main gene identifier.
## This is not necessary; the point is to demonstrate that it does not matter what 
## IDs we use; it only matter that the IDs provided to tmod (genes= option)
## and the IDs that were used to create the modules must match.

## first, what are the unique go ids:
goids <- unique(bm$go_id)
length(goids)

## prepare a mapping between go ids and illumina ids
## list with go ids as names, and vectors containing Illumina IDs as values
go2gene <- split(bm$illumina_humanwg_6_v3, f = bm$go_id)

head(lengths(go2gene))
max(lengths(go2gene))

## remove all GO terms with less than 20 genes or more than 500
go2gene <- go2gene[ lengths(go2gene) >= 20 ]
go2gene <- go2gene[ lengths(go2gene) <= 250 ]
length(go2gene)

goterms <- data.frame(ID=names(go2gene), 
   Title=bm$name_1006[ match(names(go2gene), bm$go_id)])

## gene2go is the reverse mapping - from genes to GO terms
gene2go <- split(rep(names(go2gene), lengths(go2gene)), unlist(go2gene))
mygo <- makeTmod(modules=goterms, modules2genes = go2gene, genes2modules = gene2go)

res <- tmodLimmaTest(fit, genes=fit$genes$ID, mset=mygo)
