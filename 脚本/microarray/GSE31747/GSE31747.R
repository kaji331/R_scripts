# 数据读入

setwd("~/Projects/microarray_data/GSE31747/CELs")
system("ls *.CEL.gz | xargs -n1 7z e")
library(affy)
gse <- ReadAffy()
gse_sample <- read.csv("../Samples.csv")
sampleNames(gse) <- sub("\\.CEL$","",sampleNames(gse))
mt <- match(gse_sample$SampleID,sampleNames(gse))
vmd <- data.frame(labelDescription=c("Sample ID","Individuals","Condition","Time"))
phenoData(gse) <- new("AnnotatedDataFrame",data=gse_sample[mt,],varMetadata=vmd)
sampleNames(gse) <- sampleNames(gse@protocolData)

group_color <- factor(c(rep(c("lightblue","pink"),6))) %>>% as.character

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

# 怀疑降解严重，RNA降解图就不放了
# par(mar=c(5,3,3,0))
# plotRNAdeg <- function(x){
#   pp <- function(x,col){
#     plotAffyRNAdeg(x,col=col)
#     legend("topleft",rownames(pData(gse)),
#            col=col,lwd=1,ncol=4,inset=0.05,cex=0.5)
#   }
#   pp(x,rainbow(12))
# }
# AffyRNAdeg(gse) %>>% plotRNAdeg

# 均一化
gse_rma <- rma(gse)
par(mar=c(8,3,2,2))
boxplot(gse_rma,col=group_color,las=3,main="RMA")

# 注释
fData(gse_rma) <- featureNames(gse_rma) %>>% getSYMBOL(gse_annotation) %>>%
  (data.frame(Symbol=.))

# 选取差异表达基因

library(limma)
gse_eset <- exprs(gse_rma)
levels <- factor(paste(gse_sample$Time,gse_sample$Condition,sep="_"))
# 执行后看看contrasts(levels)后理解意义并修正
contrasts(levels) <- cbind(Time=c(0,0,1,1),one_hr=c(1,0,0,0),six_hr=c(0,0,1,0))
design <- model.matrix(~levels)
colnames(design) <- c("Intercept","Time","Ebola_infected_1hr",
                      "Ebola_infected_6hr")
fit <- lmFit(gse_eset,design)
cont.matrix <- cbind(Ebola_infected_1hr=c(0,0,1,0),
                     Ebola_infected_6hr=c(0,0,0,1))
fit2 <- contrasts.fit(fit,cont.matrix)
fit3 <- eBayes(fit2)

grep("AFFX",featureNames(gse_rma)) %>>% (summary(fit3$F.p.value[.]))
results <- classifyTestsF(fit3,p.value=0.05)

# 区分上调基因和下调基因并加上基因名注释
results_up <- getSYMBOL(rownames(results[results[,1]==1 & results[,2]==1,]),gse_annotation) %>>%
  (.[!is.na(.)])
results_down <- getSYMBOL(rownames(results[results[,1]==-1 & results[,2]==-1,]),gse_annotation) %>>%
  (.[!is.na(.)])
results_1hr_up <- getSYMBOL(rownames(results[results[,1]==1,]),gse_annotation) %>>%
  (.[!is.na(.)])
results_1hr_down <- getSYMBOL(rownames(results[results[,1]==-1,]),gse_annotation) %>>%
  (.[!is.na(.)])
results_6hr_up <- getSYMBOL(rownames(results[results[,2]==1,]),gse_annotation) %>>%
  (.[!is.na(.)])
results_6hr_down <- getSYMBOL(rownames(results[results[,2]==-1,]),gse_annotation) %>>%
  (.[!is.na(.)])
# 合并基因id注释
results_up <- cbind(results_up,getEG(names(results_up),gse_annotation))
results_down <- cbind(results_down,getEG(names(results_down),gse_annotation))
results_1hr_up <- cbind(results_1hr_up,getEG(names(results_1hr_up),gse_annotation))
results_1hr_down <- cbind(results_1hr_down,getEG(names(results_1hr_down),gse_annotation))
results_6hr_up <- cbind(results_6hr_up,getEG(names(results_6hr_up),gse_annotation))
results_6hr_down <- cbind(results_6hr_down,getEG(names(results_6hr_down),gse_annotation))
# 合并芯片检测标签
results_up <- cbind(results_up,rownames(results_up))
results_down <- cbind(results_down,rownames(results_down))
results_1hr_up <- cbind(results_1hr_up,rownames(results_1hr_up))
results_1hr_down <- cbind(results_1hr_down,rownames(results_1hr_down))
results_6hr_up <- cbind(results_6hr_up,rownames(results_6hr_up))
results_6hr_down <- cbind(results_6hr_down,rownames(results_6hr_down))
# 从矩阵转换为数据框
d_up <- as.data.frame(results_up)
d_down <- as.data.frame(results_down)
d_1hr_up <- as.data.frame(results_1hr_up)
d_1hr_down <- as.data.frame(results_1hr_down)
d_6hr_up <- as.data.frame(results_6hr_up)
d_6hr_down <- as.data.frame(results_6hr_down)
# 添加列名
colnames(d_up) <- c("Symbol","ID","Label")
colnames(d_down) <- c("Symbol","ID","Label")
colnames(d_1hr_up) <- c("Symbol","ID","Label")
colnames(d_1hr_down) <- c("Symbol","ID","Label")
colnames(d_6hr_up) <- c("Symbol","ID","Label")
colnames(d_6hr_down) <- c("Symbol","ID","Label")
# 提取表达平均值，F值，和F检验的p值
d_temp_up <- fit3[rownames(d_up),] %>>% as.data.frame
d_temp_down <- fit3[rownames(d_down),] %>>% as.data.frame
d_temp_up <- d_temp_up[,c("Amean","F","F.p.value")]
d_temp_down <- d_temp_down[,c("Amean","F","F.p.value")]

d_1hr_temp_up <- fit3[rownames(d_1hr_up),] %>>% as.data.frame
d_1hr_temp_down <- fit3[rownames(d_1hr_down),] %>>% as.data.frame
d_1hr_temp_up <- d_1hr_temp_up[,c("Amean","F","F.p.value")]
d_1hr_temp_down <- d_1hr_temp_down[,c("Amean","F","F.p.value")]

d_6hr_temp_up <- fit3[rownames(d_6hr_up),] %>>% as.data.frame
d_6hr_temp_down <- fit3[rownames(d_6hr_down),] %>>% as.data.frame
d_6hr_temp_up <- d_6hr_temp_up[,c("Amean","F","F.p.value")]
d_6hr_temp_down <- d_6hr_temp_down[,c("Amean","F","F.p.value")]
# ===
d_up <- cbind(d_up,d_temp_up)
d_up <- d_up[order(d_up$F.p.value),]
d_down <- cbind(d_down,d_temp_down)
d_down <- d_down[order(d_down$F.p.value),]

d_1hr_up <- cbind(d_1hr_up,d_1hr_temp_up)
d_1hr_up <- d_1hr_up[order(d_1hr_up$F.p.value),]
d_1hr_down <- cbind(d_1hr_down,d_1hr_temp_down)
d_1hr_down <- d_1hr_down[order(d_1hr_down$F.p.value),]

d_6hr_up <- cbind(d_6hr_up,d_6hr_temp_up)
d_6hr_up <- d_6hr_up[order(d_6hr_up$F.p.value),]
d_6hr_down <- cbind(d_6hr_down,d_6hr_temp_down)
d_6hr_down <- d_6hr_down[order(d_6hr_down$F.p.value),]
# 去掉上下调不确定的基因（包含了重复的symbol和id）
d <- rbind(d_up[!(d_up$Symbol %in% d_down$Symbol),],
           d_down[!(d_down$Symbol %in% d_up$Symbol),])
d <- d[order(d$F.p.value),]
d_1hr <- rbind(d_1hr_up[!(d_1hr_up$Symbol %in% d_1hr_down$Symbol),],
              d_1hr_down[!(d_1hr_down$Symbol %in% d_1hr_up$Symbol),])
d_1hr <- d_1hr[order(d_1hr$F.p.value),]
d_6hr <- rbind(d_6hr_up[!(d_6hr_up$Symbol %in% d_6hr_down$Symbol),],
               d_6hr_down[!(d_6hr_down$Symbol %in% d_6hr_up$Symbol),])
d_6hr <- d_6hr[order(d_6hr$F.p.value),]

# 整理去掉重复值
duplicated(d_1hr$Symbol) %>>% (d_1hr[.,])
duplicated(d_6hr$Symbol) %>>% (d_6hr[.,])
# 1hr中TNF重复，6hr中IL1B重复
d_1hr <- d_1hr[d_1hr$Label != "259_s_at",]
d_1hr_up <- d_1hr_up[d_1hr_up$Label != "259_s_at",]
d_6hr <- d_6hr[d_6hr$Label != "39402_at",]
d_6hr_up <- d_6hr_up[d_6hr_up$Label != "39402_at",]

# 输出为EXCEL
library(XLConnect)
writeWorksheetToFile("../dif.xlsx",data=d,sheet="Sheet1",startRow=1,startCol=1)
writeWorksheetToFile("../dif_1hr.xlsx",data=d_1hr,sheet="Sheet1",startRow=1,startCol=1)
writeWorksheetToFile("../dif_6hr.xlsx",data=d_6hr,sheet="Sheet1",startRow=1,startCol=1)

# 热图

library(pheatmap)
# 感染1hr和6hr共有差异表达基因热图
selected <- gse_eset[rownames(d),]
rownames(selected) <- d$Symbol
pheatmap(selected,color=colorRampPalette(c("green","black","red"))(100),border_color=NA)
# 感染1hr差异表达基因热图
selected <- gse_eset[rownames(d_1hr),]
rownames(selected) <- d_1hr$Symbol
pheatmap(selected,color=colorRampPalette(c("green","black","red"))(100),border_color=NA)
# 感染6hr差异表达基因热图
selected <- gse_eset[rownames(d_6hr),]
rownames(selected) <- d_6hr$Symbol
pheatmap(selected,color=colorRampPalette(c("green","black","red"))(100),border_color=NA)

# 统计分析及可视化

# 添加logFC值
t <- topTable(fit3,coef=c("Ebola_infected_1hr","Ebola_infected_6hr"),number=100)
fold_d <- t[rownames(t) %in% rownames(d),c("Ebola_infected_1hr","Ebola_infected_6hr")]
rownames(fold_d) <- getEG(rownames(fold_d),gse_annotation)
fold_d <- apply(fold_d,1,mean)

fold_d_1hr <- t[rownames(t) %in% rownames(d_1hr),c("Ebola_infected_1hr","Ebola_infected_6hr")]
rownames(fold_d_1hr) <- getEG(rownames(fold_d_1hr),gse_annotation)
fold_d_1hr <- apply(fold_d_1hr,1,mean)

fold_d_6hr <- t[rownames(t) %in% rownames(d_6hr),c("Ebola_infected_1hr","Ebola_infected_6hr")]
rownames(fold_d_6hr) <- getEG(rownames(fold_d_6hr),gse_annotation)
fold_d_6hr <- apply(fold_d_6hr,1,mean)


library(clusterProfiler)
library(DOSE)
# 感染1hr和6hr共有差异表达基因GO分析
ego_up <- enrichGO(gene=unique(as.character(d_up$ID)),organism="human",ont="BP",readable=T)
ego_down <- enrichGO(gene=unique(as.character(d_down$ID)),organism="human",ont="BP",readable=T,
                     minGSSize=1)
cnetplot(ego_up,fixed=F,categorySize="pvalue",foldChange=fold_d)
cnetplot(ego_down,fixed=F,categorySize="pvalue",foldChange=fold_d)
# 感染1hr差异表达基因GO分析
ego_1hr_up <- enrichGO(gene=unique(as.character(d_1hr_up$ID)),organism="human",ont="BP",readable=T)
ego_1hr_down <- enrichGO(gene=unique(as.character(d_1hr_down$ID)),organism="human",ont="BP",readable=T)
cnetplot(ego_1hr_up,fixed=F,categorySize="pvalue",foldChange=fold_d_1hr)
cnetplot(ego_1hr_down,fixed=F,categorySize="pvalue",foldChange=fold_d_1hr)
# 感染6hr差异表达基因GO分析
ego_6hr_up <- enrichGO(gene=unique(as.character(d_6hr_up$ID)),organism="human",ont="BP",readable=T)
ego_6hr_down <- enrichGO(gene=unique(as.character(d_6hr_down$ID)),organism="human",ont="BP",readable=T,
                         pvalueCutoff=0.05,pAdjustMethod="none",qvalueCutoff=1,minGSSize=3)
cnetplot(ego_6hr_up,fixed=F,categorySize="pvalue",foldChange=fold_d_6hr)
cnetplot(ego_6hr_down,fixed=F,categorySize="pvalue",foldChange=fold_d_6hr)

summary(ego_up) %>>% 
  (writeWorksheetToFile("../ego_up.xlsx",data=.,sheet="Sheet1",startRow=1,startCol=1))
summary(ego_down) %>>% 
  (writeWorksheetToFile("../ego_down.xlsx",data=.,sheet="Sheet1",startRow=1,startCol=1))
summary(ego_1hr_up) %>>% 
  (writeWorksheetToFile("../ego_1hr_up.xlsx",data=.,sheet="Sheet1",startRow=1,startCol=1))
summary(ego_1hr_down) %>>% 
  (writeWorksheetToFile("../ego_1hr_down.xlsx",data=.,sheet="Sheet1",startRow=1,startCol=1))
summary(ego_6hr_up) %>>% 
  (writeWorksheetToFile("../ego_6hr_up.xlsx",data=.,sheet="Sheet1",startRow=1,startCol=1))
summary(ego_6hr_down) %>>% 
  (writeWorksheetToFile("../ego_6hr_down.xlsx",data=.,sheet="Sheet1",startRow=1,startCol=1))

# 分别对上调和下调的基因进行KEGG分析
# 由于KEGG缺乏维护，现在开始流行Reactome分析。
# 感染1hr和6hr共有差异表达基因KEGG分析
ekg_up <- enrichKEGG(gene=unique(as.character(d_up$ID)),organism="human",readable=T,
                     pAdjustMethod="none",minGSSize=1)
# ekg_down <- enrichKEGG(gene=unique(as.character(d_down$ID)),organism="human",readable=T,
#                        pAdjustMethod="none",minGSSize=1,qvalueCutoff=1,pvalueCutoff=1)
cnetplot(ekg_up,fixed=F,categorySize="pvalue",foldChange=fold_d)
# 感染1hr差异表达基因KEGG分析
ekg_1hr_up <- enrichKEGG(gene=unique(as.character(d_1hr_up$ID)),organism="human",readable=T,
                         pAdjustMethod="none",minGSSize=1)
ekg_1hr_down <- enrichKEGG(gene=unique(as.character(d_1hr_down$ID)),organism="human",readable=T,
                           pAdjustMethod="none",minGSSize=1)
cnetplot(ekg_1hr_up,fixed=F,categorySize="pvalue",foldChange=fold_d_1hr)
cnetplot(ekg_1hr_down,fixed=F,categorySize="pvalue",foldChange=fold_d_1hr)
# 感染6hr差异表达基因KEGG分析
ekg_6hr_up <- enrichKEGG(gene=unique(as.character(d_6hr_up$ID)),organism="human",readable=T,
                         pAdjustMethod="none",minGSSize=1)
# ekg_6hr_down <- enrichKEGG(gene=unique(as.character(d_6hr_down$ID)),organism="human",readable=T,
#                            pAdjustMethod="none",minGSSize=1,qvalueCutoff=1,pvalueCutoff=1)
cnetplot(ekg_6hr_up,fixed=F,categorySize="pvalue",foldChange=fold_d_6hr)

summary(ekg_up) %>>% 
  (writeWorksheetToFile("../ekg_up.xlsx",data=.,sheet="Sheet1",startRow=1,startCol=1))
summary(ekg_1hr_up) %>>% 
  (writeWorksheetToFile("../ekg_1hr_up.xlsx",data=.,sheet="Sheet1",startRow=1,startCol=1))
summary(ekg_1hr_down) %>>% 
  (writeWorksheetToFile("../ekg_1hr_down.xlsx",data=.,sheet="Sheet1",startRow=1,startCol=1))
summary(ekg_6hr_up) %>>% 
  (writeWorksheetToFile("../ekg_6hr_up.xlsx",data=.,sheet="Sheet1",startRow=1,startCol=1))

library(ReactomePA)
# 感染1hr和6hr共有差异表达基因Reactome分析
epa_up <- enrichPathway(gene=unique(as.character(d_up$ID)),organism="human",readable=T,
                        pAdjustMethod="none",minGSSize=1)
# epa_down <- enrichPathway(gene=unique(as.character(d_down$ID)),organism="human",readable=T,
#                           pAdjustMethod="none",minGSSize=1,qvalueCutoff=1,pvalueCutoff=1)
cnetplot(epa_up,fixed=F,categorySize="pvalue",foldChange=fold_d)
# 感染1hr差异表达基因Reactome分析
epa_1hr_up <- enrichPathway(gene=unique(as.character(d_1hr_up$ID)),organism="human",readable=T,
                            pAdjustMethod="none",minGSSize=1)
epa_1hr_down <- enrichPathway(gene=unique(as.character(d_1hr_down$ID)),organism="human",readable=T,
                              pAdjustMethod="none",minGSSize=1)
cnetplot(epa_1hr_up,fixed=F,categorySize="pvalue",foldChange=fold_d_1hr)
cnetplot(epa_1hr_down,fixed=F,categorySize="pvalue",foldChange=fold_d_1hr)
# 感染6hr差异表达基因Reactome分析
epa_6hr_up <- enrichPathway(gene=unique(as.character(d_6hr_up$ID)),organism="human",readable=T,
                            pAdjustMethod="none",minGSSize=1)
# epa_6hr_down <- enrichPathway(gene=unique(as.character(d_6hr_down$ID)),organism="human",readable=T,
#                               pAdjustMethod="none",minGSSize=1,qvalueCutoff=1,pvalueCutoff=1)
cnetplot(epa_6hr_up,fixed=F,categorySize="pvalue",foldChange=fold_d_6hr)

summary(epa_up) %>>% 
  (writeWorksheetToFile("../epa_up.xlsx",data=.,sheet="Sheet1",startRow=1,startCol=1))
summary(epa_1hr_up) %>>% 
  (writeWorksheetToFile("../epa_1hr_up.xlsx",data=.,sheet="Sheet1",startRow=1,startCol=1))
summary(epa_1hr_down) %>>% 
  (writeWorksheetToFile("../epa_1hr_down.xlsx",data=.,sheet="Sheet1",startRow=1,startCol=1))
summary(epa_6hr_up) %>>% 
  (writeWorksheetToFile("../epa_6hr_up.xlsx",data=.,sheet="Sheet1",startRow=1,startCol=1))
