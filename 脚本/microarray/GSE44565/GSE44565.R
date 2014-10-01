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
results13.5_up <- getSYMBOL(rownames(results[results[,1]==1,]),gse_annotation)
results13.5_down <- getSYMBOL(rownames(results[results[,1]==-1,]),gse_annotation)
results15.5_up <- getSYMBOL(rownames(results[results[,2]==1,]),gse_annotation)
results15.5_down <- getSYMBOL(rownames(results[results[,2]==-1,]),gse_annotation)
# 合并基因id注释
results_up <- cbind(results_up,getEG(names(results_up),gse_annotation))
results_down <- cbind(results_down,getEG(names(results_down),gse_annotation))
results13.5_up <- cbind(results13.5_up,getEG(names(results13.5_up),gse_annotation))
results13.5_down <- cbind(results13.5_down,getEG(names(results13.5_down),gse_annotation))
results15.5_up <- cbind(results15.5_up,getEG(names(results15.5_up),gse_annotation))
results15.5_down <- cbind(results15.5_down,getEG(names(results15.5_down),gse_annotation))
# 合并芯片检测标签
results_up <- cbind(results_up,rownames(results_up))
results_down <- cbind(results_down,rownames(results_down))
results13.5_up <- cbind(results13.5_up,rownames(results13.5_up))
results13.5_down <- cbind(results13.5_down,rownames(results13.5_down))
results15.5_up <- cbind(results15.5_up,rownames(results15.5_up))
results15.5_down <- cbind(results15.5_down,rownames(results15.5_down))

# 剔除没有注释基因名称和基因id的行
results_up <- results_up[!is.na(results_up[,1]) & !is.na(results_up[,2]),]
results_down <- results_down[!is.na(results_down[,1]) & !is.na(results_down[,2]),]
results13.5_up <- results13.5_up[!is.na(results13.5_up[,1]) & !is.na(results13.5_up[,2]),]
results13.5_down <- results13.5_down[!is.na(results13.5_down[,1]) & !is.na(results13.5_down[,2]),]
results15.5_up <- results15.5_up[!is.na(results15.5_up[,1]) & !is.na(results15.5_up[,2]),]
results15.5_down <- results15.5_down[!is.na(results15.5_down[,1]) & !is.na(results15.5_down[,2]),]
# 从矩阵转换为数据框
d_up <- as.data.frame(results_up)
d_down <- as.data.frame(results_down)
d_13.5_up <- as.data.frame(results13.5_up)
d_13.5_down <- as.data.frame(results13.5_down)
d_15.5_up <- as.data.frame(results15.5_up)
d_15.5_down <- as.data.frame(results15.5_down)
# 添加列名
colnames(d_up) <- c("Symbol","ID","Label")
colnames(d_down) <- c("Symbol","ID","Label")
colnames(d_13.5_up) <- c("Symbol","ID","Label")
colnames(d_13.5_down) <- c("Symbol","ID","Label")
colnames(d_15.5_up) <- c("Symbol","ID","Label")
colnames(d_15.5_down) <- c("Symbol","ID","Label")
# 提取表达平均值，F值，和F检验的p值
d_temp_up <- fit2[rownames(d_up),] %>>% as.data.frame
d_temp_down <- fit2[rownames(d_down),] %>>% as.data.frame
d_temp_up <- d_temp_up[,c("Amean","F","F.p.value")]
d_temp_down <- d_temp_down[,c("Amean","F","F.p.value")]

d_13.5_temp_up <- fit2[rownames(d_13.5_up),] %>>% as.data.frame
d_13.5_temp_down <- fit2[rownames(d_13.5_down),] %>>% as.data.frame
d_13.5_temp_up <- d_13.5_temp_up[,c("Amean","F","F.p.value")]
d_13.5_temp_down <- d_13.5_temp_down[,c("Amean","F","F.p.value")]

d_15.5_temp_up <- fit2[rownames(d_15.5_up),] %>>% as.data.frame
d_15.5_temp_down <- fit2[rownames(d_15.5_down),] %>>% as.data.frame
d_15.5_temp_up <- d_15.5_temp_up[,c("Amean","F","F.p.value")]
d_15.5_temp_down <- d_15.5_temp_down[,c("Amean","F","F.p.value")]
# ===
d_up <- cbind(d_up,d_temp_up)
d_up <- d_up[order(d_up$F.p.value),]
d_down <- cbind(d_down,d_temp_down)
d_down <- d_down[order(d_down$F.p.value),]

d_13.5_up <- cbind(d_13.5_up,d_13.5_temp_up)
d_13.5_up <- d_13.5_up[order(d_13.5_up$F.p.value),]
d_13.5_down <- cbind(d_13.5_down,d_13.5_temp_down)
d_13.5_down <- d_13.5_down[order(d_13.5_down$F.p.value),]

d_15.5_up <- cbind(d_15.5_up,d_15.5_temp_up)
d_15.5_up <- d_15.5_up[order(d_15.5_up$F.p.value),]
d_15.5_down <- cbind(d_15.5_down,d_15.5_temp_down)
d_15.5_down <- d_15.5_down[order(d_15.5_down$F.p.value),]
# 去掉上下调不确定的基因（包含了重复的symbol和id）
d <- rbind(d_up[!(d_up$Symbol %in% d_down$Symbol),],
           d_down[!(d_down$Symbol %in% d_up$Symbol),])
d <- d[order(d$F.p.value),]
d_13.5 <- rbind(d_13.5_up[!(d_13.5_up$Symbol %in% d_13.5_down$Symbol),],
                d_13.5_down[!(d_13.5_down$Symbol %in% d_13.5_up$Symbol),])
d_13.5 <- d_13.5[order(d_13.5$F.p.value),]
d_15.5 <- rbind(d_15.5_up[!(d_15.5_up$Symbol %in% d_15.5_down$Symbol),],
                d_15.5_down[!(d_15.5_down$Symbol %in% d_15.5_up$Symbol),])
d_15.5 <- d_15.5[order(d_15.5$F.p.value),]

# 热图

library(pheatmap)
# 13.5天和15.5天共有差异表达基因热图
selected <- gse_eset[rownames(d),]
rownames(selected) <- d$Symbol
pheatmap(selected,color=colorRampPalette(c("green","black","red"))(100),border_color=NA)
# 13.5天基因差异表达热图
selected <- gse_eset[rownames(d_13.5),
                     c("GSM795342","GSM795343","GSM795344","GSM795345","GSM795346","GSM795347")]
rownames(selected) <- d_13.5$Symbol
pheatmap(selected,color=colorRampPalette(c("green","black","red"))(100),border_color=NA)
# 15.5天基因差异表达热图
selected <- gse_eset[rownames(d_15.5),
                     c("GSM795348","GSM795349","GSM795350","GSM795351","GSM795352","GSM795353")]
rownames(selected) <- d_15.5$Symbol
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
