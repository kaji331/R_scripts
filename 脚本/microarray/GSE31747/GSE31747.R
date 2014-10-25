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

# 输出为EXCEL
library(XLConnect)

