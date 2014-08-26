# GSE32078
# 载入数据

setwd("~/Projects/microarray_data/GSE32078/CELs")
system("ls *.CEL.gz | xargs -n1 7z e")
library(affy)
gse <- ReadAffy()
gse_sample <- read.csv("../Samples.csv")
sampleNames(gse) <- sub("\\.CEL$","",sampleNames(gse))
mt <- match(gse_sample$SampleID,sampleNames(gse))
vmd <- data.frame(labelDescription=c("Sample ID","Time","Status","Repeat times"))
phenoData(gse) <- new("AnnotatedDataFrame",data=gse_sample[mt,],varMetadata=vmd)
sampleNames(gse) <- sampleNames(gse@protocolData)

group_color <- (factor(c(rep("lightblue",6),rep("pink",6)))) %>>% as.character

library(annotate)
gse_annotation <- annPkgName(gse@annotation,type="db")
if(!require(gse_annotation,character.only=T)) {
  biocLite(gse_annotation)
  library(gse_annotation,character.only=T)
}

# 质量检验

library(simpleaffy)
qc(gse) %>>% plot

library(affyPLM)
par(mar=c(8,3,2,2))
fitPLM(gse) %>>% {
  Mbox(.,main="RLE",col=group_color,las=3,ylim=c(-0.4,0.4),whisklty=0)
  boxplot(.,main="NUSE",col=group_color,las=3,ylim=c(0.95,1.1),whisklty=0)
}

par(mar=c(5,3,3,0))
plotRNAdeg <- function(x){
  pp <- function(x,col){
    plotAffyRNAdeg(x,col=col)
    legend("topleft",rownames(pData(gse)),
           col=col,lwd=1,ncol=4,inset=0.05,cex=0.5)
  }
  pp(x,rainbow(12))
}
AffyRNAdeg(gse) %>>% plotRNAdeg

# 均一化

gse_rma <- rma(gse)
par(mar=c(8,3,2,2))
boxplot(gse_rma,col=group_color,las=3,main="RMA")

# 选取差异表达基因

library(limma)
gse_eset <- exprs(gse_rma)
levels <- factor(paste(gse_sample$Time,gse_sample$Status,sep="_"))
# 这里可以做完了看一下就知道contrasts的意义了
contrasts(levels) <- cbind(Time=c(0,0,1,1),day_13.5=c(0,1,0,0),day_15.5=c(0,0,0,1))
design <- model.matrix(~levels)
colnames(design) <- c("Intercept","Time","day_13.5","day_15.5")
fit <- lmFit(gse_eset,design)
cont.matrix <- cbind(day_13.5=c(0,0,1,0),day_15.5=c(0,0,0,1))
fit2 <- contrasts.fit(fit,cont.matrix)
fit2 <- eBayes(fit2)
# F p-value选取，可以使用global或者nestedF，这里使用另外更加保守的方法不依赖数据的分布情况
# results <- decideTests(fit2,method="global")
# results <- decideTests(fit2,method="nestedF")
grep("AFFX",featureNames(gse_rma)) %>>% summary(fit2$F.p.value[.]) # Min=0.08336，因此p.value取0.08
results <- classifyTestsF(fit2,p.value=0.08)
# 各分组上下调检测数
summary(results)
# 分组间上下调检测交叉数
table(day_13.5=results[,1],day_15.5=results[,2])
# 分组间上下调检测交叉图
vennDiagram(results,include="up")
vennDiagram(results,include="down")

# 区分上调基因和下调基因并加上基因名注释
results_up <- getSYMBOL(rownames(results[results[,1]==1 & results[,2]==1,]),gse_annotation)
results_down <- getSYMBOL(rownames(results[results[,1]==-1 & results[,2]==-1,]),gse_annotation)
# results13.5_up <- getSYMBOL(rownames(results[results[,1]==1,]),gse_annotation)
# results13.5_down <- getSYMBOL(rownames(results[results[,1]==-1,]),gse_annotation)
# results15.5_up <- getSYMBOL(rownames(results[results[,2]==1,]),gse_annotation)
# results15.5_down <- getSYMBOL(rownames(results[results[,2]==-1,]),gse_annotation)
# 合并基因id注释
results_up <- cbind(results_up,getEG(names(results_up),gse_annotation))
results_down <- cbind(results_down,getEG(names(results_down),gse_annotation))
# results13.5_up <- cbind(results13.5_up,getEG(names(results13.5_up),gse_annotation))
# results13.5_down <- cbind(results13.5_down,getEG(names(results13.5_down),gse_annotation))
# results15.5_up <- cbind(results15.5_up,getEG(names(results15.5_up),gse_annotation))
# results15.5_down <- cbind(results15.5_down,getEG(names(results15.5_down),gse_annotation))
# 合并芯片检测标签
results_up <- cbind(results_up,rownames(results_up))
results_down <- cbind(results_down,rownames(results_down))
# results13.5_up <- cbind(results13.5_up,rownames(results13.5_up))
# results13.5_down <- cbind(results13.5_down,rownames(results13.5_down))
# results15.5_up <- cbind(results15.5_up,rownames(results15.5_up))
# results15.5_down <- cbind(results15.5_down,rownames(results15.5_down))

# 剔除没有注释基因名称和基因id的行
results_up <- results_up[!is.na(results_up[,1]) & !is.na(results_up[,2]),]
results_down <- results_down[!is.na(results_down[,1]) & !is.na(results_down[,2]),]
# results13.5_up <- results13.5_up[!is.na(results13.5_up[,1]) & !is.na(results13.5_up[,2]),]
# results13.5_down <- results13.5_down[!is.na(results13.5_down[,1]) & !is.na(results13.5_down[,2]),]
# results15.5_up <- results15.5_up[!is.na(results15.5_up[,1]) & !is.na(results15.5_up[,2]),]
# results15.5_down <- results15.5_down[!is.na(results15.5_down[,1]) & !is.na(results15.5_down[,2]),]
# 从矩阵转换为数据框
d_up <- as.data.frame(results_up)
d_down <- as.data.frame(results_down)
# d_13.5_up <- as.data.frame(results13.5_up)
# d_13.5_down <- as.data.frame(results13.5_down)
# d_15.5_up <- as.data.frame(results15.5_up)
# d_15.5_down <- as.data.frame(results15.5_down)
# 添加列名
colnames(d_up) <- c("Symbol","ID","Label")
colnames(d_down) <- c("Symbol","ID","Label")
# colnames(d_13.5_up) <- c("Symbol","ID","Label")
# colnames(d_13.5_down) <- c("Symbol","ID","Label")
# colnames(d_15.5_up) <- c("Symbol","ID","Label")
# colnames(d_15.5_down) <- c("Symbol","ID","Label")
# 提取表达平均值，F值，和F检验的p值
d_temp_up <- fit2[rownames(d_up),] %>>% as.data.frame
d_temp_down <- fit2[rownames(d_down),] %>>% as.data.frame
d_temp_up <- d_temp_up[,c("Amean","F","F.p.value")]
d_temp_down <- d_temp_down[,c("Amean","F","F.p.value")]

d_up <- cbind(d_up,d_temp_up)
d_up <- d_up[order(d_up$F.p.value),]
d_down <- cbind(d_down,d_temp_down)
d_down <- d_down[order(d_down$F.p.value),]
# 去掉上下调不确定的基因（包含了重复的symbol和id）
d <- rbind(d_up[!(d_up$Symbol %in% d_down$Symbol),],
           d_down[!(d_down$Symbol %in% d_up$Symbol),])
d <- d[order(d$F.p.value),]
# d_13.5 <- rbind(d_13.5_up[!(d_13.5_up$Symbol %in% d_13.5_down$Symbol),],
#                 d_13.5_down[!(d_13.5_down$Symbol %in% d_13.5_up$Symbol),])
# d_15.5 <- rbind(d_15.5_up[!(d_15.5_up$Symbol %in% d_15.5_down$Symbol),],
#                 d_15.5_down[!(d_15.5_down$Symbol %in% d_15.5_up$Symbol),])

# 热图

library(pheatmap)
selected <- gse_eset[rownames(d),]
rownames(selected) <- d$Symbol
pheatmap(selected,color=colorRampPalette(c("green","black","red"))(100),border_color=NA)
# # 13.5天基因差异表达热图
# selected <- gse_eset[rownames(d_13.5),
#                      c("GSM795342","GSM795343","GSM795344","GSM795345","GSM795346","GSM795347")]
# rownames(selected) <- d_13.5$Symbol
# pheatmap(selected,color=colorRampPalette(c("green","black","red"))(100),border_color=NA)
# # 15.5天基因差异表达热图
# selected <- gse_eset[rownames(d_15.5),
#                      c("GSM795348","GSM795349","GSM795350","GSM795351","GSM795352","GSM795353")]
# rownames(selected) <- d_15.5$Symbol
# pheatmap(selected,color=colorRampPalette(c("green","black","red"))(100),border_color=NA)

# 统计分析及可视化

# 换成clusterProfiler包，出柱状图
library(clusterProfiler)
library(DOSE)
ego_up <- enrichGO(gene=unique(as.character(d_up$ID)),organism="mouse",ont="BP",pvalueCutoff=0.05,pAdjustMethod="none",readable=T)
ego_down <- enrichGO(gene=unique(as.character(d_down$ID)),organism="mouse",ont="BP",pvalueCutoff=0.05,pAdjustMethod="none",readable=T)
plot(ego_up)
plot(ego_down)

# library(GOstats)
# entrezUniverse <- unique(unlist(mget(rownames(gse_eset),mouse4302ENTREZID)))
# entrezSelected <- unique(as.character(d$ID))
# params <- new("GOHyperGParams",geneIds=entrezSelected,universeGeneIds=entrezUniverse,annotation=gse_annotation,ontology="BP",pvalueCutoff=0.001,conditional=F,testDirection="over")
# hgOver <- hyperGTest(params)
# bp <- summary(hgOver)
# htmlReport(hgOver,file="../ALL_go.html")
# library(Rgraphviz)
# ghandle <- goDag(hgOver)
# subGHandle <- subGraph(snodes=as.character(summary(hgOver)[,1]),graph=ghandle)
# plot(subGHandle)

# 分别对上调和下调的基因进行KEGG分析
# 由于KEGG缺乏维护，现在开始流行Reactome分析。不过ReactomePA和clusterProfiler包分析enrichPathway和enrichKEGG结果数量明显少于GeneAnswers包，还有赖于未来的项目验证。
ekg_up <- enrichKEGG(gene=unique(as.character(d_up$ID)),organism="mouse",pvalueCutoff=1,minGSSize=1,qvalueCutoff=1,readable=T)
ekg_down <- enrichKEGG(gene=unique(as.character(d_down$ID)),organism="mouse",pvalueCutoff=1,minGSSize=1,qvalueCutoff=1,readable=T)
cnetplot(ekg_up,fixed=F)
cnetplot(ekg_down,fixed=F)

library(ReactomePA)
epa_up <- enrichPathway(gene=unique(as.character(d_up$ID)),organism="mouse",pvalueCutoff=1,qvalueCutoff=1,minGSSize=1,readable=T)
epa_down <- enrichPathway(gene=unique(as.character(d_down$ID)),organism="mouse",pvalueCutoff=1,qvalueCutoff=1,minGSSize=1,readable=T)
cnetplot(epa_up,fixed=F)
cnetplot(epa_down,fixed=F)

# library(GeneAnswers)
# humanGeneInput_up <- d_up[,c("ID","F","F.p.value")]
# colnames(humanGeneInput_up) <- c("GeneID","F","pvalue")
# humanExpr_up <- gse_eset[match(rownames(d_up),rownames(gse_eset)),]
# humanExpr_up <- cbind(d_up$ID,humanExpr_up)
# y_up <- geneAnswersBuilder(humanGeneInput_up,"org.Mm.eg.db",categoryType="KEGG",testType="hyperG",pvalueT=0.1,geneExpressionProfile=humanExpr_up,verbose=F)
# yy_up <- geneAnswersReadable(y_up,verbose=F)
# geneAnswersConceptNet(yy_up,colorValueColumn="F",centroidSize="pvalue",output="interactive")
