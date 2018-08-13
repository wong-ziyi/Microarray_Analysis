#==== required packages check and installation ====
list.of.packages <- c("BiocInstaller", "easycsv")
new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[,"Package"])]
if(length(new.packages)) {
  install.packages(new.packages)
}
library(BiocInstaller)
library(easycsv)
list.of.packages <- c("xps", "annotate", "limma", "GenomeInfoDb")
new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[,"Package"])]
if(length(new.packages)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite(new.packages)
}
library(xps)
library(limma)
library(GenomeInfoDb)
library(XML)
library(annotate)
##### First run #####

# ++ Create a ROOT scheme file ####
MoGenome430A20 <- import.expr.scheme("AffymetrixMouseGenome430A20Array",
                                 filedir = file.path(getwd(), "scheme"),
                                 schemefile =file.path(getwd(), "lib", "Mouse430A_2.cdf"),
                                 probefile =file.path(getwd(), "lib", "Mouse430A_2_probe.tab"),
                                 annotfile =file.path(getwd(), "lib", "Mouse430A_2.na36.annot.csv"))

# ++ upload the raw files ####
celdir <- file.path(choose_dir())
celfiles <- select.list(list.files(path = celdir, pattern = ".[Cc][Ee][Ll]"), multiple = TRUE)
data.raw <- import.data(xps.scheme = MoGenome430A20, 
                        filename = "GSE42874", 
                        celdir = celdir, 
                        celfiles = celfiles, 
                        verbose = TRUE)

# Or upload previous results #
MoGenome430A20 <- root.scheme(file.path(getwd(), "scheme", "AffymetrixMouseGenome430A20Array.root"))
data.raw <- root.data(xps.scheme = MoGenome430A20, 
                      rootfile = file.path(getwd(), "GSE42874", "GSE42874_cel.root"))

# Get the list for converting probsetID to unitID #
data.raw <- attachInten(data.raw)
data.raw <- attachUnitNames(data.raw, treetype="idx")
data.raw <- attachMask(data.raw)

# ++ Background correction and expression normalize #####
data.rma.tran <- xpsRMA(data.raw, 
                     "GSE42874_rma_tran", 
                     background = "pmonly", 
                     exonlevel = "", 
                     normalize = TRUE,
                     option = "transcript",
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
pd <- data.frame(population = c( 1, 1, 1, 2, 2, 2), replicate = c(1, 2, 3, 1, 2, 3))
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
     cex =0.25, xlim = c(-2, 2), ylim = c(0, 7), 
     xlab ="", ylab="", main = "GSE42874")
title(ylab = expression( ~ "-log"[10]^{ italic(P) ~ "-value"}), line = 2)
title(xlab = expression(~ log[2]^{"Fold Change"}), line = 3)
# Then add the points with the same colors as before:
points(tableTop[p1, ]$logFC, -log10(tableTop[p1, ]$P.Value), pch = 18, col = "blue", cex=0.7)
points(tableTop[p2, ]$logFC, -log10(tableTop[p2, ]$P.Value), pch = 1, col = "red", cex=0.7)

##### Annotation ####
MoGenome430A20_pbs<-read.table(file = file.path(getwd(), "GSE42874", "Mouse430A_2_idx.txt"), sep = "\t", header = TRUE)
anno<-read.csv(file.path(getwd(), "lib", "Mouse430A_2.na36.annot.csv"), skip = 22, header = TRUE)
GeneID<-as.vector(unitID2symbol(data.raw, rownames(tableTop), as.list=FALSE))
length(GeneID) == length(rownames(tableTop))

##### Output ####
res<-cbind(UnitID = rownames(tableTop), symbol = GeneID, tableTop)
write.csv(res, file = "Results.csv", row.names = FALSE)
res_DEGs <- res[intersect(p1,p2),]
write.csv(res_DEGs, file = "Results_DEGs.csv", row.names = FALSE)
