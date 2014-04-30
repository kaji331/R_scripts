#说明

1. 全自动。只需要设置工作路径以及filelist.txt文件即可。

2. filelist.txt，必须与CEL文件一起放置在工作目录内，文件以**空格为间隔**，分两列，分别是filename和factor 例

```{txt}
filename factor
MS1_(Mouse430_2).cel control
MS2_(Mouse430_2).cel control
MS3_(Mouse430_2).cel control
MS4_(Mouse430_2).cel beta-catenin
MS5_(Mouse430_2).cel beta-catenin
MS6_(Mouse430_2).cel beta-catenin
MS7_(Mouse430_2).cel TCF4
MS8_(Mouse430_2).cel TCF4
MS9_(Mouse430_2).cel TCF4
MS10_(Mouse430_2).cel Icat
MS11_(Mouse430_2).cel Icat
MS12_(Mouse430_2).cel Icat
```

```{r}
source("http://pgfe.umassmed.edu/ou/bioconductor/affyChipAna.R")
setwd("~/Documents/Rscript/testData")
doAll()
```

更多参数：
- dataInput(file=”filelist.txt”,sep=” “,rnames=”filename”)
- getExpr(data)
- QCreport(data,fastQCfile=”fastQC.pdf”,arrayQMfile=”QCreport”)
- getDesign(data)
- getContrastMatrix(data,design,constr=NULL)
- doAnalyze(data,design,contrast.matrix)
- outputData(data,fit2,contrast.matrix,times=2,cut.p.val=0.01,cut.fdr.val=NULL,outputFolder=”output”,all=FALSE,scale=”row”,cluster_rows=TRUE,cluster_cols=TRUE)
- doAll(file=”filelist.txt”,workingdir=”.”,sep=” “,rnames=”filename”,fastQCfile=”fastQC.pdf”,arrayQMfile=”QCreport”,
- costr=NULL,times=2,cut.p.val=0.01,cut.fdr.val=NULL,outputFolder=”output”,
all=FALSE,scale=”row”,doQC=TRUE,cluster_rows=TRUE,cluster_cols=TRUE)

可以指定`design & contrast.matrix: design and contrast.matrix for limma, coef: report contrast, times: folder-change, cut.fdr.val: P.adj.Val cutoff value`

source:[原代码](http://pgfe.umassmed.edu/ou/bioconductor/affyChipAna.R)
