# Data prep.
# David L Gibbs
# June 1, 2016

# this data comes from Array Express.
# https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-1908/
# E-MTAB-1908 - Transcription profiling by array of Saccharomyces cerevisiae cells to estimate labeled and total mRNA levels every 5 minutes for three complete cell cycles


# from bioconductor
library(affy)

# the cel file information matrix.
inf <- read.table("E-MTAB-1908.sdrf.txt", sep="\t", header=T, stringsAsFactors=F)

# make sure the cel file information matches the file listing.
list.celfiles() # this is what needs to line up with the cel information matrix.
all(inf$Array.Data.File == list.celfiles())
#[1] TRUE

# give row names to the information matrix
rownames(inf) <- inf$Array.Data.File

# read the cel files
eset <- justRMA(filenames=list.celfiles(), phenoData=inf)

# extract the matrix
exprmat <- exprs(eset)

# to get gene symbols...
# this is the affymetrix annotation for the Yeast2 array, release 31
# filtered down by Species.Scientific.Name == "Saccharomyces cerevisiae"
annot <- read.table("yeast2_array_affymetrix_annotation.csv.gz", sep=",", header=T, stringsAsFactors=F)

# now we filter down the expression matrix
exprmat <- exprmat[rownames(exprmat) %in% annot$Probe.Set.ID,]

# now getting the annotation into the same order as the expression matrix.
annot_sorted <- annot[match(table=annot$Probe.Set.ID, x=rownames(exprmat)),]

# for comparison
eser <- read.table(".../miergolf/data/Eser_Averaged_Expression.txt.gz", sep="\t", header=F, stringsAsFactors=F)
eser_genes <- reeser <- read.table(".../miergolf/data/Eser_Averaged_Expression.txt.gz", sep="\t", header=F, stringsAsFactors=F)ad.delim("~/Dropbox/Research/Projects/Influence_Maximization_Problem/miergolf/data/yeast_array_genesymbols.csv.gz", sep="\t", header=F, stringsAsFactors=F)
prev_annot <- read.table(".../miergolf/data/yeast_array_annotation.csv", sep=",", header=T, stringsAsFactors=F)

all(annot_sorted$Gene.Symbol == eser_genes$V1)
#[1] TRUE

all(annot_sorted$Probe.Set.ID == rownames(exprmat))
# True

sum(annot_sc$ProbeIDs %in% prev_annot$Probe.Set.ID)
#[1] 5814

idx1 <- 85:125
idx2 <- 126:166
x <- exprmat[,idx1]
y <- exprmat[,idx2]
eser_data_used_in_this_study <- (x+y)/2

#done#
