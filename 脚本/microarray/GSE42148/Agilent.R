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
tT <- topTable(fit2, number=1615)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","GeneName"))
tT <- tT[!is.na(tT$GeneName),]
dup <- duplicated(tT$GeneName)
tT <- tT[!dup,]
tT <- tT[!(regexpr("^LOC[0-9]{4}",tT$GeneName) == 1),]
tT <- tT[!(regexpr("^lincRNA",tT$GeneName) == 1),]
tT <- tT[!(regexpr("^ERCC-[0-9]",tT$GeneName) == 1),]
tT <- tT[!(regexpr("^ENST[0-9]",tT$GeneName) == 1),]
tT <- tT[!(regexpr("^A_[0-9]",tT$GeneName) == 1),]
tT <- tT[!(regexpr("^RP[0-9]+-[0-9]+",tT$GeneName) == 1),]
tT <- tT[!(regexpr("^BE[0-9]{4}",tT$GeneName) == 1),]
tT <- tT[!(regexpr("^FLJ[0-9]{4}",tT$GeneName) == 1),]
tT <- tT[!(regexpr("DarkCorner",tT$GeneName) == 1),]
tT <- tT[!(regexpr("^hCG_[0-9]{4}",tT$GeneName) == 1),]
tT <- tT[!(regexpr("^ETG[0-9]+_[0-9]+",tT$GeneName) == 1),]
tT <- tT[!(regexpr("^CB[0-9]{4}",tT$GeneName) == 1),]
tT <- tT[!(regexpr("^AK[0-9]{4}",tT$GeneName) == 1),]
tT <- tT[!(regexpr("^CR[0-9]{4}",tT$GeneName) == 1),]
tT <- tT[!(regexpr("^BQ[0-9]{4}",tT$GeneName) == 1),]
tT <- tT[!(regexpr("^NP[0-9]{4}",tT$GeneName) == 1),]
# write.table(tT, file=stdout(), row.names=F, sep="\t")

# 上下调差异基因识别
tT_up <- tT[tT$logFC > 0,]
tT_down <- tT[tT$logFC < 0,]


# 火山图
t <- topTable(fit2, number=99999)
t <- subset(t, select=c("ID","adj.P.Val","P.Value","t","B","logFC","GeneName"))
t <- t[!is.na(t$GeneName),]
color <- rep("lightblue",62972)
t <- cbind(t,color)
t$color <- as.character(t$color)
# 抽样示意，以加快绘图速度
# set.seed(1984)
# a <- sample(t$GeneName,3000)
# t <- t[t$GeneName[a],]
t[t$logFC > 1 & -log10(t$P.Value) > 2,]$color <- "red"
t[t$logFC < -1 & -log10(t$P.Value) > 2,]$color <- "green"
plot(t$logFC,-log10(t$P.Value),col=t$color,pch=20)
abline(h=2) # p<0.01，若p<0.05则h=1.30103
abline(v=1)
abline(v=-1)

# 热图
library(pheatmap)
selected <- ex[rownames(tT),]
rownames(selected) <- tT$GeneName
pheatmap(selected,color=colorRampPalette(c("green","black","red"))(100),border_color=NA)

# 换成clusterProfiler包，出柱状图
library(clusterProfiler)
ego_up <- enrichGO(gene=tT_up$ID,ont="BP",pvalueCutoff=0.05,readable=T)
ego_down <- enrichGO(gene=tT_down$ID,ont="BP",pvalueCutoff=0.05,readable=T)
barplot(ego_up)
barplot(ego_down)

library(GeneAnswers)
humanGeneInput_up <- tT_up[,c("ID","B","adj.P.Val")]
humanExpr_up <- ex[match(rownames(tT_up),rownames(ex)),]
humanExpr_up <- cbind(humanGeneInput_up[,"ID"],humanExpr_up)
humanGeneInput_up <- humanGeneInput_up[!is.na(humanGeneInput_up[,1]),]
humanExpr_up <- humanExpr_up[!is.na(humanExpr_up[,1]),]
y_up <- geneAnswersBuilder(humanGeneInput_up,"org.Hs.eg.db",categoryType="KEGG",testType="hyperG",pvalueT=0.05,geneExpressionProfile=humanExpr_up,verbose=F)
yy_up <- geneAnswersReadable(y_up,verbose=F)
geneAnswersConceptNet(yy_up,colorValueColumn="B",centroidSize="pvalue",output="interactive")
yyy_up <- geneAnswersSort(yy_up,sortBy="pvalue")
geneAnswersHeatmap(yyy_up)

humanGeneInput_down <- tT_down[,c("ID","B","adj.P.Val")]
humanExpr_down <- ex[match(rownames(tT_down),rownames(ex)),]
humanExpr_down <- cbind(humanGeneInput_down[,"ID"],humanExpr_down)
humanGeneInput_down <- humanGeneInput_down[!is.na(humanGeneInput_down[,1]),]
humanExpr_down <- humanExpr_down[!is.na(humanExpr_down[,1]),]
y_down <- geneAnswersBuilder(humanGeneInput_down,"org.Hs.eg.db",categoryType="KEGG",testType="hyperG",pvalueT=0.05,geneExpressionProfile=humanExpr_down,verbose=F)
yy_down <- geneAnswersReadable(y_down,verbose=F)
geneAnswersConceptNet(yy_down,colorValueColumn="B",centroidSize="pvalue",output="interactive")
yyy_down <- geneAnswersSort(yy_down,sortBy="pvalue")
geneAnswersHeatmap(yyy_down)
