# GSE44565
# 载入数据

setwd("~/Projects/microarray_data/GSE44565/CELs")
library(limma)
gse <- readTargets("../GSE44565 sample.txt") %>>% 
  read.maimages(source="agilent",green.only=T) %>>%
  backgroundCorrect(method="normexp",offset=16) %>>%
  normalizeBetweenArrays(method="quantile") %>>%
  (avereps(.,ID=.$genes$ProbeName))
gse@.Data[[2]]$FileName <- sub("_(.*)","",gse@.Data[[2]]$FileName)
attr(gse@.Data[[1]],"dimnames")[[2]] <- gse@.Data[[2]]$FileName

group_color <- (factor(c(rep("lightblue",3),rep("pink",13)))) %>>% as.character
par(mar=c(8,3,2,2))
boxplot(gse$E,col=group_color,las=3,main="quantile")

# 选取差异表达基因

f <- factor(gse@.Data[[2]]$Condition,levels=unique(gse@.Data[[2]]$Condition))
design <- model.matrix(~0+f)
colnames(design) <- levels(f)
contrast.matrix <- makeContrasts("infected - ctrl",levels=design)
fit <- lmFit(gse,design) %>>% contrasts.fit(contrast.matrix) %>>% eBayes
results <- topTable(fit,coef="infected - ctrl",genelist=gse$genes,p.value=0.01)
results <- results[,c("ProbeUID","GeneName","logFC","AveExpr","t","adj.P.Val","B")]

# 热图

library(pheatmap)
selected <- gse[rownames(results),]
rownames(selected) <- results$GeneName
pheatmap(selected,color=colorRampPalette(c("green","black","red"))(100),border_color=NA)

# 统计分析及可视化

library(GOstats)
entrezUniverse <- unique(unlist(mget(rownames(gse_eset),hgu133aENTREZID)))
entrezSelected <- unique(d[!is.na(d$EntrezID),"EntrezID"])
params <- new("GOHyperGParams",geneIds=entrezSelected,universeGeneIds=entrezUniverse,annotation=gse_annotation,ontology="BP",pvalueCutoff=0.001,conditional=F,testDirection="over")
hgOver <- hyperGTest(params)
bp <- summary(hgOver)
htmlReport(hgOver,file="../ALL_go.html")
library(Rgraphviz)
ghandle <- goDag(hgOver)
subGHandle <- subGraph(snodes=as.character(summary(hgOver)[,1]),graph=ghandle)
plot(subGHandle)

library(GeneAnswers)
humanGeneInput_up <- d_up[,c("EntrezID","F","adj.P.Val")]
humanExpr_up <- gse_eset[match(rownames(d_up),rownames(gse_eset)),]
humanExpr_up <- cbind(humanGeneInput_up[,"EntrezID"],humanExpr_up)
humanGeneInput_up <- humanGeneInput_up[!is.na(humanGeneInput_up[,1]),]
humanExpr_up <- humanExpr_up[!is.na(humanExpr_up[,1]),]
y_up <- geneAnswersBuilder(humanGeneInput_up,"org.Hs.eg.db",categoryType="KEGG",testType="hyperG",pvalueT=0.01,geneExpressionProfile=humanExpr_up,verbose=F)
yy_up <- geneAnswersReadable(y_up,verbose=F)
geneAnswersConceptNet(yy_up,colorValueColumn="F",centroidSize="pvalue",output="interactive")
yyy_up <- geneAnswersSort(yy_up,sortBy="pvalue")
geneAnswersHeatmap(yyy_up)

# # 换成clusterProfiler包，出柱状图
# library(clusterProfiler)
# library(DOSE)
# ego <- enrichGO(gene=unique(as.character(results$ProbeUID)),organism="mouse",ont="BP",pvalueCutoff=0.05,pAdjustMethod="none",readable=T)
# ego_down <- enrichGO(gene=unique(as.character(d_down$ID)),organism="mouse",ont="BP",pvalueCutoff=0.05,pAdjustMethod="none",readable=T)
# plot(ego_up)
# plot(ego_down)
# temp <- d_up$Amean
# names(temp) <- d_up$ID
# cnetplot(ego_up,fixed=F,categorySize="pvalue",foldChange=temp)
# 
# # 分别对上调和下调的基因进行KEGG分析
# # 由于KEGG缺乏维护，现在开始流行Reactome分析。不过ReactomePA和clusterProfiler包分析enrichPathway和enrichKEGG结果数量明显少于GeneAnswers包，还有赖于未来的项目验证。
# ekg_up <- enrichKEGG(gene=unique(as.character(d_up$ID)),organism="mouse",pvalueCutoff=1,minGSSize=1,qvalueCutoff=1,readable=T)
# ekg_down <- enrichKEGG(gene=unique(as.character(d_down$ID)),organism="mouse",pvalueCutoff=1,minGSSize=1,qvalueCutoff=1,readable=T)
# temp <- d_up$Amean
# names(temp) <- d_up$ID
# cnetplot(ekg_up,fixed=F,categorySize="pvalue",foldChange=temp)
# temp <- -d_down$Amean
# names(temp) <- d_down$ID
# cnetplot(ekg_down,fixed=F,categorySize="pvalue",foldChange=temp)
# 
# library(ReactomePA)
# epa_up <- enrichPathway(gene=unique(as.character(d_up$ID)),organism="mouse",pvalueCutoff=1,qvalueCutoff=1,minGSSize=1,readable=T)
# epa_down <- enrichPathway(gene=unique(as.character(d_down$ID)),organism="mouse",pvalueCutoff=1,qvalueCutoff=1,minGSSize=1,readable=T)
# temp <- d_up$Amean
# names(temp) <- d_up$ID
# cnetplot(epa_up,fixed=F,categorySize="pvalue",foldChange=temp)
# temp <- -d_down$Amean
# names(temp) <- d_down$ID
# cnetplot(epa_down,fixed=F,categorySize="pvalue",foldChange=temp)
