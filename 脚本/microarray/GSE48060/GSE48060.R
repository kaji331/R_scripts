# GSE48060
# 载入数据

setwd("~/Projects/microarray_data/GSE48060/CELs")
system("ls *.cel.gz | xargs -n1 7z e")
library(affy)
gse <- ReadAffy()
gse_sample <- read.csv("../GSE48060 Samples.csv")
gse_sample <- cbind(gse_sample,Group=c("yes",rep("no",5),"yes",rep("no",2),
                                       "yes",rep("no",6),"yes",rep("no",8),
                                       "yes",rep("no",4),rep("yes",10),
                                       rep("no",10),"no","yes"))
gse_sample$Condition <- c(rep("Disease",30),rep("Control",20),"Disease","Control")
sampleNames(gse) <- sub("\\.cel$","",sampleNames(gse))
sampleNames(gse)[21] <- "GSM1167092_Nel010142-HG-U133Plus2"
mt <- match(gse_sample$SampleID,sampleNames(gse))
vmd <- data.frame(labelDescription=c("Sample ID","Tissue","Condition","Label","Group"))
phenoData(gse) <- new("AnnotatedDataFrame",data=gse_sample[mt,],varMetadata=vmd)
sampleNames(gse) <- sampleNames(gse@protocolData)
sampleNames(gse) <- sub("_Nel0101(.*)","",sampleNames(gse))

group_color <- (factor(c(rep("lightblue",30),rep("pink",20),"lightblue","pink"))) %>>%
  as.character

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

# 注释
fData(gse_rma) <- featureNames(gse_rma) %>>% getSYMBOL(gse_annotation) %>>%
  (data.frame(Symbol=.))

# 选取差异表达基因

library(limma)
gse_eset <- exprs(gse_rma)
levels <- factor(paste(gse_sample$Group,gse_sample$Condition,sep="_"))
contrasts(levels) <- cbind(Group=c(0,0,1,1),no=c(0,1,0,0),yes=c(0,0,0,1))
design <- model.matrix(~levels)
colnames(design) <- c("Intercept","Group","patient_without_recurrent_events",
                      "patient_with_recurrent_events")
fit <- lmFit(gse_eset,design)
cont.matrix <- cbind(patient_without_recurrent_events=c(0,0,1,0),
                     patient_with_recurrent_events=c(0,0,0,1))
fit2 <- contrasts.fit(fit,cont.matrix)
fit3 <- eBayes(fit2)
# F p-value选取，可以使用global或者nestedF，或者使用另外更加保守的方法不依赖数据的分布情况
# results <- decideTests(fit3,method="global")
# results <- decideTests(fit3,method="nestedF")
grep("AFFX",featureNames(gse_rma)) %>>% (summary(fit3$F.p.value[.]))
# Min=0.0000010，因此p.value取0.01
results <- classifyTestsF(fit3,p.value=0.01)
# 各分组上下调检测数
summary(results)
# 分组间上下调检测交叉数
table(patient_without_recurrent_events=results[,1],
      patient_with_recurrent_events=results[,2])
# 分组间上下调检测交叉图
vennDiagram(results)
vennDiagram(results,include="up")
vennDiagram(results,include="down")

# 区分上调基因和下调基因并加上基因名注释
results_up <- getSYMBOL(rownames(results[results[,1]==1 & results[,2]==1,]),gse_annotation) %>>%
  (.[!is.na(.)])
results_down <- getSYMBOL(rownames(results[results[,1]==-1 & results[,2]==-1,]),gse_annotation) %>>%
  (.[!is.na(.)])
results_no_up <- getSYMBOL(rownames(results[results[,1]==1,]),gse_annotation) %>>%
  (.[!is.na(.)])
results_no_down <- getSYMBOL(rownames(results[results[,1]==-1,]),gse_annotation) %>>%
  (.[!is.na(.)])
results_yes_up <- getSYMBOL(rownames(results[results[,2]==1,]),gse_annotation) %>>%
  (.[!is.na(.)])
results_yes_down <- getSYMBOL(rownames(results[results[,2]==-1,]),gse_annotation) %>>%
  (.[!is.na(.)])
# 合并基因id注释
results_up <- cbind(results_up,getEG(names(results_up),gse_annotation))
results_down <- cbind(results_down,getEG(names(results_down),gse_annotation))
results_no_up <- cbind(results_no_up,getEG(names(results_no_up),gse_annotation))
results_no_down <- cbind(results_no_down,getEG(names(results_no_down),gse_annotation))
results_yes_up <- cbind(results_yes_up,getEG(names(results_yes_up),gse_annotation))
results_yes_down <- cbind(results_yes_down,getEG(names(results_yes_down),gse_annotation))
# 合并芯片检测标签
results_up <- cbind(results_up,rownames(results_up))
results_down <- cbind(results_down,rownames(results_down))
results_no_up <- cbind(results_no_up,rownames(results_no_up))
results_no_down <- cbind(results_no_down,rownames(results_no_down))
results_yes_up <- cbind(results_yes_up,rownames(results_yes_up))
results_yes_down <- cbind(results_yes_down,rownames(results_yes_down))
# 从矩阵转换为数据框
d_up <- as.data.frame(results_up)
d_down <- as.data.frame(results_down)
d_no_up <- as.data.frame(results_no_up)
d_no_down <- as.data.frame(results_no_down)
d_yes_up <- as.data.frame(results_yes_up)
d_yes_down <- as.data.frame(results_yes_down)
# 添加列名
colnames(d_up) <- c("Symbol","ID","Label")
colnames(d_down) <- c("Symbol","ID","Label")
colnames(d_no_up) <- c("Symbol","ID","Label")
colnames(d_no_down) <- c("Symbol","ID","Label")
colnames(d_yes_up) <- c("Symbol","ID","Label")
colnames(d_yes_down) <- c("Symbol","ID","Label")
# 提取表达平均值，F值，和F检验的p值
d_temp_up <- fit3[rownames(d_up),] %>>% as.data.frame
d_temp_down <- fit3[rownames(d_down),] %>>% as.data.frame
d_temp_up <- d_temp_up[,c("Amean","F","F.p.value")]
d_temp_down <- d_temp_down[,c("Amean","F","F.p.value")]

d_no_temp_up <- fit3[rownames(d_no_up),] %>>% as.data.frame
d_no_temp_down <- fit3[rownames(d_no_down),] %>>% as.data.frame
d_no_temp_up <- d_no_temp_up[,c("Amean","F","F.p.value")]
d_no_temp_down <- d_no_temp_down[,c("Amean","F","F.p.value")]

d_yes_temp_up <- fit3[rownames(d_yes_up),] %>>% as.data.frame
d_yes_temp_down <- fit3[rownames(d_yes_down),] %>>% as.data.frame
d_yes_temp_up <- d_yes_temp_up[,c("Amean","F","F.p.value")]
d_yes_temp_down <- d_yes_temp_down[,c("Amean","F","F.p.value")]
# ===
d_up <- cbind(d_up,d_temp_up)
d_up <- d_up[order(d_up$F.p.value),]
d_down <- cbind(d_down,d_temp_down)
d_down <- d_down[order(d_down$F.p.value),]

d_no_up <- cbind(d_no_up,d_no_temp_up)
d_no_up <- d_no_up[order(d_no_up$F.p.value),]
d_no_down <- cbind(d_no_down,d_no_temp_down)
d_no_down <- d_no_down[order(d_no_down$F.p.value),]

d_yes_up <- cbind(d_yes_up,d_yes_temp_up)
d_yes_up <- d_yes_up[order(d_yes_up$F.p.value),]
d_yes_down <- cbind(d_yes_down,d_yes_temp_down)
d_yes_down <- d_yes_down[order(d_yes_down$F.p.value),]
# 去掉上下调不确定的基因（包含了重复的symbol和id）
d <- rbind(d_up[!(d_up$Symbol %in% d_down$Symbol),],
           d_down[!(d_down$Symbol %in% d_up$Symbol),])
d <- d[order(d$F.p.value),]
d_no <- rbind(d_no_up[!(d_no_up$Symbol %in% d_no_down$Symbol),],
                d_no_down[!(d_no_down$Symbol %in% d_no_up$Symbol),])
d_no <- d_no[order(d_no$F.p.value),]
d_yes <- rbind(d_yes_up[!(d_yes_up$Symbol %in% d_yes_down$Symbol),],
                d_yes_down[!(d_yes_down$Symbol %in% d_yes_up$Symbol),])
d_yes <- d_yes[order(d_yes$F.p.value),]

# 输出为EXCEL
library(XLConnect)

# 热图

library(pheatmap)
# patient_with_recurrent_events相对于与patient_without_recurrent_events差异表达基因热图
dif <- d_yes[!(d_yes$Symbol %in% d$Symbol),]
dif_up <- d_yes_up[!(d_yes_up$Symbol %in% d$Symbol),]
dif_up <- dif_up[dif_up$Symbol != "OSGIN2",] #排除上下调重复的
dif_down <- d_yes_down[!(d_yes_down$Symbol %in% d$Symbol),]
dif_down <- dif_down[dif_down$Symbol != "OSGIN2",]

writeWorksheetToFile("../dif.xlsx",data=dif,sheet="Sheet1",startRow=1,startCol=1)
writeWorksheetToFile("../dif_up.xlsx",data=dif_up,sheet="Sheet1",startRow=1,startCol=1)
writeWorksheetToFile("../dif_down.xlsx",data=dif_down,sheet="Sheet1",startRow=1,startCol=1)

selected <- gse_eset[rownames(dif),]
rownames(selected) <- NULL
pheatmap(selected,color=colorRampPalette(c("green","black","red"))(100),border_color=NA)

# 统计分析及可视化

# 换成clusterProfiler包，出柱状图
library(clusterProfiler)
library(DOSE)

# 复发较未复发相对对照表达差异显著基因GO分析
ego_dif <- enrichGO(gene=unique(as.character(dif$ID)),organism="human",ont="BP",readable=T)
ego_dif_up <- enrichGO(gene=unique(as.character(dif_up$ID)),organism="human",ont="BP",readable=T)
ego_dif_down <- enrichGO(gene=unique(as.character(dif_down$ID)),organism="human",ont="BP",readable=T)
plot(ego_dif)
plot(ego_dif_up)
plot(ego_dif_down)
summary(ego_dif) %>>% 
  (writeWorksheetToFile("../ego_dif.xlsx",data=.,sheet="Sheet1",startRow=1,startCol=1))
summary(ego_dif_up) %>>%
  (writeWorksheetToFile("../ego_dif_up.xlsx",data=.,sheet="Sheet1",startRow=1,startCol=1))
summary(ego_dif_down) %>>%
  (writeWorksheetToFile("../ego_dif_down.xlsx",data=.,sheet="Sheet1",startRow=1,startCol=1))

cnetplot(ego_dif,fixed=F,categorySize="pvalue")
cnetplot(ego_dif_up,fixed=F,categorySize="pvalue")
cnetplot(ego_dif_down,fixed=F,categorySize="pvalue")

# KEGG和Reactome分析
# 有问题
# t <- topTable(fit3,coef="patient_with_recurrent_events",number=4100)
# fold <- t[rownames(t) %in% rownames(dif),"logFC"]
# names(fold) <- rownames(t[rownames(t) %in% rownames(dif),])
# names(fold) <- getEG(names(fold),gse_annotation)
fold_up <- fold[fold > 0]
fold_down <- fold[fold < 0]

ekg_dif <- enrichKEGG(gene=unique(as.character(dif$ID)),organism="human",readable=T)
ekg_dif_up <- enrichKEGG(gene=unique(as.character(dif_up$ID)),organism="human",readable=T,minGSSize=2)
ekg_dif_down <- enrichKEGG(gene=unique(as.character(dif_down$ID)),organism="human",readable=T)
plot(ekg_dif)
plot(ekg_dif_up)
plot(ekg_dif_down)
summary(ekg_dif) %>>%
  (writeWorksheetToFile("../ekg_dif.xlsx",data=.,sheet="Sheet1",startRow=1,startCol=1))
summary(ekg_dif_up) %>>%
  (writeWorksheetToFile("../ekg_dif_up.xlsx",data=.,sheet="Sheet1",startRow=1,startCol=1))
summary(ekg_dif_down) %>>%
  (writeWorksheetToFile("../ekg_dif_down.xlsx",data=.,sheet="Sheet1",startRow=1,startCol=1))

cnetplot(ekg_dif,fixed=F,categorySize="pvalue",foldChange=fold)
cnetplot(ekg_dif_up,fixed=F,categorySize="pvalue",foldChange=fold_up)
cnetplot(ekg_dif_down,fixed=F,categorySize="pvalue",foldChange=fold_down)

library(ReactomePA)
epa_dif <- enrichPathway(gene=unique(as.character(dif$ID)),organism="human",readable=T,pAdjustMethod="none")
epa_dif_up <- enrichPathway(gene=unique(as.character(dif_up$ID)),organism="human",readable=T,pAdjustMethod="none")
epa_dif_down <- enrichPathway(gene=unique(as.character(dif_down$ID)),organism="human",readable=T,pAdjustMethod="none")
plot(epa_dif)
plot(epa_dif_up)
plot(epa_dif_down)
summary(epa_dif) %>>%
  (writeWorksheetToFile("../epa_dif.xlsx",data=.,sheet="Sheet1",startRow=1,startCol=1))
summary(epa_dif_up) %>>%
  (writeWorksheetToFile("../epa_dif_up.xlsx",data=.,sheet="Sheet1",startRow=1,startCol=1))
summary(epa_dif_down) %>>%
  (writeWorksheetToFile("../epa_dif_down.xlsx",data=.,sheet="Sheet1",startRow=1,startCol=1))

cnetplot(epa_dif,fixed=F,categorySize="pvalue",foldChange=fold)
cnetplot(epa_dif_up,fixed=F,categorySize="pvalue",foldChange=fold_up)
cnetplot(epa_dif_down,fixed=F,categorySize="pvalue",foldChange=fold_down)
