# GSE44565
# 载入数据

setwd("~/Projects/microarray_data/GSE44565/CELs")
library(limma)
gse <- readTargets("../GSE44565 sample.txt") %>>% 
  read.maimages(source="agilent",green.only=T) %>>%
  backgroundCorrect(method="normexp",offset=16) %>>%
  normalizeBetweenArrays(method="quantile") %>>%
  (avereps(.,ID=.$genes$ProbeName))

group_color <- (factor(c(rep("lightblue",3),rep("pink",13)))) %>>% as.character
boxplot(gse$E,col=group_color,las=3,main="quantile")

# 选取差异表达基因

f <- factor(gse@.Data[[2]]$Condition,levels=unique(gse@.Data[[2]]$Condition))
design <- model.matrix(~0+f)
colnames(design) <- levels(f)
contrast.matrix <- makeContrasts("infected - ctrl",levels=design)
fit <- lmFit(gse,design) %>>% contrasts.fit(contrast.matrix) %>>% eBayes
results <- topTable(fit,coef="infected - ctrl",genelist=gse$genes,p.value=0.01)
results <- results[,c("GeneName","logFC","AveExpr","t","adj.P.Val","B")]

# 热图

library(pheatmap)
selected <- gse[rownames(results),]
rownames(selected) <- gse$GeneName
pheatmap(selected,color=colorRampPalette(c("green","black","red"))(100),border_color=NA)

# 统计分析及可视化

# 换成clusterProfiler包，出柱状图
library(clusterProfiler)
library(DOSE)
ego_up <- enrichGO(gene=unique(as.character(d_up$ID)),organism="mouse",ont="BP",pvalueCutoff=0.05,pAdjustMethod="none",readable=T)
ego_down <- enrichGO(gene=unique(as.character(d_down$ID)),organism="mouse",ont="BP",pvalueCutoff=0.05,pAdjustMethod="none",readable=T)
plot(ego_up)
plot(ego_down)
temp <- d_up$Amean
names(temp) <- d_up$ID
cnetplot(ego_up,fixed=F,categorySize="pvalue",foldChange=temp)
# 13.5 cnet图太大太复杂
ego_13.5_up <- enrichGO(gene=unique(as.character(d_13.5_up$ID)),organism="mouse",ont="BP",pvalueCutoff=0.05,pAdjustMethod="none",readable=T)
ego_13.5_down <- enrichGO(gene=unique(as.character(d_13.5_down$ID)),organism="mouse",ont="BP",pvalueCutoff=0.05,pAdjustMethod="none",readable=T)
plot(ego_13.5_up)
plot(ego_13.5_down)
# 15.5 cnet图太大太复杂
ego_15.5_up <- enrichGO(gene=unique(as.character(d_15.5_up$ID)),organism="mouse",ont="BP",pvalueCutoff=0.05,pAdjustMethod="none",readable=T)
ego_15.5_down <- enrichGO(gene=unique(as.character(d_15.5_down$ID)),organism="mouse",ont="BP",pvalueCutoff=0.05,pAdjustMethod="none",readable=T)
plot(ego_15.5_up)
plot(ego_15.5_down)

# 分别对上调和下调的基因进行KEGG分析
# 由于KEGG缺乏维护，现在开始流行Reactome分析。不过ReactomePA和clusterProfiler包分析enrichPathway和enrichKEGG结果数量明显少于GeneAnswers包，还有赖于未来的项目验证。
ekg_up <- enrichKEGG(gene=unique(as.character(d_up$ID)),organism="mouse",pvalueCutoff=1,minGSSize=1,qvalueCutoff=1,readable=T)
ekg_down <- enrichKEGG(gene=unique(as.character(d_down$ID)),organism="mouse",pvalueCutoff=1,minGSSize=1,qvalueCutoff=1,readable=T)
temp <- d_up$Amean
names(temp) <- d_up$ID
cnetplot(ekg_up,fixed=F,categorySize="pvalue",foldChange=temp)
temp <- -d_down$Amean
names(temp) <- d_down$ID
cnetplot(ekg_down,fixed=F,categorySize="pvalue",foldChange=temp)

ekg_13.5_up <- enrichKEGG(gene=unique(as.character(d_13.5_up$ID)),organism="mouse",pvalueCutoff=1,minGSSize=1,qvalueCutoff=1,readable=T)
ekg_13.5_down <- enrichKEGG(gene=unique(as.character(d_13.5_down$ID)),organism="mouse",pvalueCutoff=1,minGSSize=1,qvalueCutoff=1,readable=T)
temp <- d_13.5_up$Amean
names(temp) <- d_13.5_up$ID
cnetplot(ekg_13.5_up,fixed=F,categorySize="pvalue",foldChange=temp)
temp <- -d_13.5_down$Amean
names(temp) <- d_13.5_down$ID
cnetplot(ekg_13.5_down,fixed=F,categorySize="pvalue",foldChange=temp)

ekg_15.5_up <- enrichKEGG(gene=unique(as.character(d_15.5_up$ID)),organism="mouse",pvalueCutoff=1,minGSSize=1,qvalueCutoff=1,readable=T)
ekg_15.5_down <- enrichKEGG(gene=unique(as.character(d_15.5_down$ID)),organism="mouse",pvalueCutoff=1,minGSSize=1,qvalueCutoff=1,readable=T)
temp <- d_15.5_up$Amean
names(temp) <- d_15.5_up$ID
cnetplot(ekg_15.5_up,fixed=F,categorySize="pvalue",foldChange=temp)
temp <- -d_15.5_down$Amean
names(temp) <- d_15.5_down$ID
cnetplot(ekg_15.5_down,fixed=F,categorySize="pvalue",foldChange=temp)

library(ReactomePA)
epa_up <- enrichPathway(gene=unique(as.character(d_up$ID)),organism="mouse",pvalueCutoff=1,qvalueCutoff=1,minGSSize=1,readable=T)
epa_down <- enrichPathway(gene=unique(as.character(d_down$ID)),organism="mouse",pvalueCutoff=1,qvalueCutoff=1,minGSSize=1,readable=T)
temp <- d_up$Amean
names(temp) <- d_up$ID
cnetplot(epa_up,fixed=F,categorySize="pvalue",foldChange=temp)
temp <- -d_down$Amean
names(temp) <- d_down$ID
cnetplot(epa_down,fixed=F,categorySize="pvalue",foldChange=temp)

epa_13.5_up <- enrichPathway(gene=unique(as.character(d_13.5_up$ID)),organism="mouse",pvalueCutoff=1,qvalueCutoff=1,minGSSize=1,readable=T)
epa_13.5_down <- enrichPathway(gene=unique(as.character(d_13.5_down$ID)),organism="mouse",pvalueCutoff=1,qvalueCutoff=1,minGSSize=1,readable=T)
temp <- d_13.5_up$Amean
names(temp) <- d_13.5_up$ID
cnetplot(epa_13.5_up,fixed=F,categorySize="pvalue",foldChange=temp)
temp <- -d_13.5_down$Amean
names(temp) <- d_13.5_down$ID
cnetplot(epa_13.5_down,fixed=F,categorySize="pvalue",foldChange=temp)

epa_15.5_up <- enrichPathway(gene=unique(as.character(d_15.5_up$ID)),organism="mouse",readable=T)
epa_15.5_down <- enrichPathway(gene=unique(as.character(d_15.5_down$ID)),organism="mouse",readable=T)
temp <- d_15.5_up$Amean
names(temp) <- d_15.5_up$ID
cnetplot(epa_15.5_up,fixed=F,categorySize="pvalue",foldChange=temp)
temp <- -d_15.5_down$Amean
names(temp) <- d_15.5_down$ID
cnetplot(epa_15.5_down,fixed=F,categorySize="pvalue",foldChange=temp)
