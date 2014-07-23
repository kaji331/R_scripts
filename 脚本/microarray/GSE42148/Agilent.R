# Version info: R 2.14.1, Biobase 2.15.3, GEOquery 2.23.2, limma 3.10.1
# R scripts generated  Thu Jul 10 06:21:38 EDT 2014
# modified by kaji331 from GEO2R

################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

setwd("~/Projects/microarray_data/GSE42148/")

# load series and platform data from GEO
gset <- getGEO("GSE42148", GSEMatrix =TRUE)
if (length(gset) > 1) idx <- grep("GPL13607", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
sml <- as.character(factor(c(rep("Control",11),rep("CAD",13))))

#   Boxplot for selected GEO samples
# order samples by group
ex <- exprs(gset)[ , order(sml)]
sml <- sml[order(sml)]
fl <- as.factor(sml)
labels <- c("control","coronary artery disease")

# set parameters and draw the plot
palette(c("#dfeaf4","#f4dfdf", "#AABBCC"))
dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE42148", '/', annotation(gset), " selected samples", sep ='')
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, fill=palette(), bty="n")

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
            exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(CAD-Control, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
tT <- topTable(fit2, number=494)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","GeneName"))
# write.table(tT, file=stdout(), row.names=F, sep="\t")

# 上下调差异基因识别
tT_up <- tT[tT$logFC > 0,]
tT_down <- tT[tT$logFC < 0,]

# 热图
library(pheatmap)
selected <- gset[rownames(tT[str_length(tT$GeneName) < 9,]),]
rownames(selected) <- tT[str_length(tT$GeneName) < 9,]$GeneName
pheatmap(selected,color=colorRampPalette(c("green","black","red"))(100),border_color=NA)

# 换成clusterProfiler包，出柱状图
library(clusterProfiler)
ego_up <- enrichGO(gene=tT_up$ID,ont="BP",pvalueCutoff=0.01,readable=T)
ego_down <- enrichGO(gene=tT_down$ID,ont="BP",pvalueCutoff=0.01,readable=T)
