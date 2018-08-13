#==== required packages check and installation ====
list.of.packages <- c("BiocInstaller", "easycsv", "stringr")
new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[,"Package"])]
if(length(new.packages)) {
  install.packages(new.packages)
}
library(BiocInstaller)
library(easycsv)
library(stringr)
list.of.packages <- c("AnnotationDbi", "RColorBrewer", "xps", "annotate", "limma", "GenomeInfoDb", "DESeq2", "DEXSeq", "GenomicFeatures")
new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[,"Package"])]
if(length(new.packages)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite(new.packages)
}
library(xps)
library(XML)
library(limma)
library(GenomeInfoDb)
library(annotate)
library(DEXSeq)
library(GenomicFeatures)
source("functions.R")
##### First run #####
# ++ Create a ROOT scheme file ####
MoGene10st <- import.exon.scheme("AffymetrixMouseGene10STArray",
                                 filedir = file.path(getwd(), "scheme"),
                                 layoutfile=file.path(getwd(), "lib", "MoGene-1_0-st-v1.r4.clf"),
                                 schemefile=file.path(getwd(), "lib", "MoGene-1_0-st-v1.r4.pgf"),
                                 probeset=file.path(getwd(), "lib", "MoGene-1_0-st-v1.na36.mm10.probeset.csv"),
                                 transcript=file.path(getwd(),"lib", "MoGene-1_0-st-v1.na36.mm10.transcript.csv"))

# ++ upload the raw files ####
celdir <- file.path(choose_dir())
celfiles <- select.list(list.files(path = celdir, pattern = ".[Cc][Ee][Ll]"), multiple = TRUE)
data.raw <- import.data(xps.scheme = MoGene10st, 
                        filename = "GSE62128", 
                        celdir = celdir, 
                        celfiles = celfiles, 
                        verbose = TRUE)

# Or upload previous results #
MoGene10st <- root.scheme(file.path(getwd(), "scheme", "AffymetrixMouseGene10STArray.root"))
data.raw <- root.data(xps.scheme = MoGene10st, 
                      rootfile = file.path(getwd(), "GSE62128", "GSE62128_cel.root"))

# Get the list for converting probsetID to unitID #
data.raw <- attachInten(data.raw)
data.raw <- attachUnitNames(data.raw, treetype="pbs")
data.raw <- attachMask(data.raw)
MoGene10st_pbs2unitID<-data.raw@scheme@unitname

# ++ Background correction and expression normalize #####
data.rma.tran <- xpsRMA(data.raw, 
                     "GSE62128_rma_tran", 
                     background = "antigenomic", 
                     exonlevel = "all", 
                     normalize = TRUE,
                     option = "transcript",
                     add.data = TRUE,
                     verbose=TRUE)
data.rma.exon <- xpsRMA(data.raw, 
                     "GSE62128_rma_exon.root", 
                     background = "antigenomic", 
                     exonlevel = "all",
                     normalize = TRUE,
                     option = "exon", 
                     add.data = TRUE,
                     verbose=TRUE)
##### Convert transcript results to Biobase data set
res_tran<-data.rma.tran@data
rownames(res_tran)<-res_tran[,1]
res_tran<-res_tran[, - c(1,2)]
res_tran<-log2(res_tran)
#
GeneID<-as.vector(unitID2symbol(data.raw, rownames(res_tran), as.list=FALSE))
res_tran<-cbind(symbol = GeneID, res_tran)
res_tran<-res_tran[complete.cases(res_tran), ]
res_tran<-res_tran[, -1]
res_tran<-new("ExpressionSet", exprs = as.matrix(res_tran))

##### Exprimental Design ####
pd <- data.frame(population = c( 2, 2, 2, 1, 1, 1), replicate = c(1, 2, 3, 1, 2, 3))
pd
rownames(pd) <- sampleNames(res_tran)
print(pd)
v1 <- list(population = "1 is control, 2 is treatment", replicate = "arbitrary numbering")
pData(res_tran)<-pd
pData(res_tran)
varMetadata(res_tran)$labelDescription[1] <- v1[1]
varMetadata(res_tran)$labelDescription[2] <- v1[2]
varMetadata(res_tran)
design<-model.matrix(~factor(res_tran$population))
design
##### Statistical processing ####
fit <- lmFit(res_tran, design, ndups = 1, method = "robust")
ebayes <- eBayes(fit)
names(ebayes)
tableTop <- topTable(ebayes, coef=NULL, adjust.method = "BH", number = length(fit$coefficients))
p1 <- rownames(tableTop[abs(tableTop$logFC)>0.59, ])
length(p1)
p2 <- rownames(tableTop[tableTop$adj.P.Val<0.1, ])
length(p2)
p <- union(p1,p2)
length(p)
tableTop[p1,]
tableTop[p2,]
plot(tableTop[setdiff(rownames(tableTop), p), ]$logFC, 
     -log10(tableTop[setdiff(rownames(tableTop), p), ]$P.Value), 
     cex =0.25, xlim = c(-1, 1), ylim = c(0, 7), 
     xlab ="", ylab="", main = "GSE62128")
title(ylab = expression( ~ "-log"[10]^{ italic(P) ~ "-value"}), line = 2)
title(xlab = expression(~ log[2]^{"Fold Change"}), line = 3)
# Then add the points with the same colors as before:
points(tableTop[p1, ]$logFC, -log10(tableTop[p1, ]$P.Value), pch = 18, col = "blue", cex=0.7)
points(tableTop[p2, ]$logFC, -log10(tableTop[p2, ]$P.Value), pch = 1, col = "red", cex=0.7)

##### Annotation ####
MoGene10st_pbs<-read.table(file = file.path(getwd(), "GSE62128", "MoGene-1_0-st-v1_pbs.txt"), sep = "\t", header = TRUE)
anno_tran<-read.csv(file.path(getwd(), "lib", "MoGene-1_0-st-v1.na36.mm10.transcript.csv"), skip = 23, header = TRUE)
GeneID<-as.vector(unitID2symbol(data.raw, rownames(tableTop), as.list=FALSE))
length(GeneID) == length(rownames(tableTop))

##### Output ####
res<-cbind(UnitID = rownames(tableTop), symbol = GeneID, tableTop)
write.csv(res, file = "Results.csv", row.names = FALSE)
res_DEGs <- res[intersect(p1,p2),]
write.csv(res_DEGs, file = "Results_DEGs.csv", row.names = FALSE)

##### Annotate the exon results #####
anno_prob<-read.csv(file.path(getwd(), "lib", "MoGene-1_0-st-v1.na36.mm10.probeset.csv"), skip = 22, header = TRUE)
anno_prob$gene_assignment<-word(anno_prob$gene_assignment, 3)
anno_prob<-anno_prob[complete.cases(anno_prob), ]
anno_prob<-data.frame(ID=paste(anno_prob$seqname, "_", anno_prob$gene_assignment, anno_prob$strand, sep = ""), anno_prob)
anno_wzy<-anno_prob[!duplicated(anno_prob$exon_id), ]
rownames(anno_wzy)<-anno_wzy$exon_id
# Numbering the exons of each gene
anno_wzy<-wzy_numberExon(x=anno_wzy)
anno_wzy<-anno_wzy[complete.cases(anno_wzy), ]
# Read exon results and Remove duplicated rows
res_exon<-data.rma.exon@data
res_exon<-res_exon[, -1]
res_exon<-res_exon[!duplicated(res_exon$UnitName), ]
rownames(res_exon)<-res_exon$UnitName
# Annotate the exon level data
anno_exon<-wzy_anno(anno_wzy = anno_wzy, res_exon = res_exon)
length(res_exon$UnitName) == length(anno_exon@groupID_featureID$groupID)
res_exon<-data.frame(groupID=anno_exon@groupID_featureID$groupID, 
                     featureID=anno_exon@groupID_featureID$featureID, 
                     res_exon)
res_exon<-res_exon[complete.cases(res_exon), ]
length(res_exon$UnitName) == length(anno_exon@featureRanges)
names(anno_exon@featureRanges)<-paste(res_exon$groupID, res_exon$featureID, sep = ":")
##### Statistics by DEXSeq #####
# Build table of experiment condition
sampleTable = data.frame(
  row.names = c("Treat1","Treat2","Treat3", "Control1","Control2","Control3"),
  condition = c("High_Gravity", "High_Gravity", "High_Gravity", "Control","Control","Control"),
  libType = c("single-end","single-end","single-end","single-end","single-end","single-end")
)

# Build dataset
suppressPackageStartupMessages( library( "DEXSeq" ) )
dxd <- DEXSeqDataSet(
  countData = as.matrix(round(res_exon[, -c(1:3)]), rownames=FALSE),
  sampleData=sampleTable,
  design= ~ sample + exon + condition:exon,
  featureID = res_exon$featureID,
  groupID = res_exon$groupID,
  featureRanges = anno_exon@featureRanges)

genesForSubset = read.table(
  file.path(choose_dir(), "Subset8.txt"),
  stringsAsFactors=FALSE)[[1]]
dxd = dxd[geneIDs( dxd ) %in% genesForSubset,]

# Check the dataset
colData(dxd)
head( counts(dxd), 5 )
split( seq_len(ncol(dxd)), colData(dxd)$exon )
head( featureCounts(dxd), 5 )
head( rowRanges(dxd), 3 )
sampleAnnotation( dxd )

# Statistical analysis
dxd <- estimateSizeFactors( dxd )
dxd <- estimateDispersions( dxd, maxit=2000)
# Check fit
plotDispEsts( dxd )
# Test differential expression exon
dxd <- testForDEU( dxd )
# Calculate the fold change
dxd <- estimateExonFoldChanges( dxd, fitExpToVar="condition", independentFiltering=FALSE)
# Store the results to "dxr"
dxr <- DEXSeqResults( dxd, independentFiltering=FALSE)
# Output the results
write.csv(dxr, file = "Exon_Results.csv")

# Check statistic results of dxr
mcols(dxr)$description
table ( dxr$padj < 0.1 )
table ( tapply( dxr$padj < 0.1, dxr$groupID, any ) )
plotMA( dxr, cex=0.8 )
dxr@listData
# Visualize the target gene 
for (i in 1:length(genesForSubset)) {
  svg(
    filename = paste0(genesForSubset[i], ".svg"),
    width = 10,
    height = 5,
    pointsize = 12
  )
  # If you only want to check one of your target genes, please change the i value and  run the below one line code only
  plotDEXSeq( dxr, genesForSubset[i], FDR=0.1, displayTranscripts=FALSE, legend=TRUE, expression=FALSE, splicing=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
  dev.off()
}


      