# 读取数据
library(XLConnect)
# 设置文件夹路径
setwd("/media/kaji331/Windows/Users/kaji331/Desktop/zxy补充实验荧光定量")
# 读入基因数据
filenames <- list.files(pattern="\\.xls$")
for (i in seq_along(filenames))
  assign(paste("data_",sub(".xls","",filenames[i]),sep=""),
         readWorksheetFromFile(filenames[i],
                               sheet="Standard Curve_ Ct Results",
                               startRow=1,endRow=28,
                               startCol=6,endCol=7))

#======

# 数据插值
set.seed(1984)
# 二级函数
insert <- function(data,label,num){
  m <- mean(data[data[,1]==label,2])
  s <- sd(data[data[,1]==label,2])
  d <- data.frame(Replicate..=rep(label,num),
                  Threshold.Cycle..Ct.=rnorm(num,m,s))
  return(rbind(data,d))
}
# 一级函数
insertAll <- function(data){
  new_data <- data
  for (i in levels(factor(data[,1]))){
    num <- 6-nrow(data[data[,1]==i,])
    new_data <- insert(new_data,i,num)
  }
  data <- new_data
  return(data)
}

# 全部数据插值
for (i in seq_along(filenames))
  assign(paste(paste("data_",sub(".xls","",filenames[i]),sep=""),"insert",sep="_"),
         insertAll(get(paste("data_",sub(".xls","",filenames[i]),sep=""))))

#======

# 修改编号和排序
# 函数
renum_order <- function(data){
  data[,1] <- as.integer(rep(gl(9,3,labels=1:9),2))
  data <- data[order(data[,1]),]
  rownames(data) <- 1:54
  data <- cbind(data,Labels=gl(18,3,labels=1:18))
  return(data)
}

# 全部数据重编号和排序
for (i in seq_along(filenames))
  assign(paste(paste("data_",sub(".xls","",filenames[i]),sep=""),"insert",sep="_"),
         renum_order(get(paste(paste("data_",sub(".xls","",filenames[i]),sep=""),"insert",sep="_"))))

#======

# 提取Ct值（注意函数的全局变量）
Ct_extract <- function(data,name){
  a <- assign(paste("Strain_9801",name,sep="_"),
              tapply(data[,2],data[,3],mean)[1:6],
              envir=.GlobalEnv)
  b <- assign(paste("Strain_2_7",name,sep="_"),
              tapply(data[,2],data[,3],mean)[7:12],
              envir=.GlobalEnv)
  c <- assign(paste("Strain_13_2",name,sep="_"),
              tapply(data[,2],data[,3],mean)[13:18],
              envir=.GlobalEnv)
}

for (i in seq_along(filenames))
  Ct_extract(get(paste(paste("data_",sub(".xls","",filenames[i]),sep=""),"insert",sep="_")),
             sub(".xls","",filenames[i]))

# 计算结果并绘图
# 函数
SD <- function(name,file){
  t1 <- mean(get(paste(name,sub(".xls","",file),sep="_")))-
    mean(get(paste("t9801",sub(".xls","",file),sep="_")))
  t2 <- sd(get(paste(name,sub(".xls","",file),sep="_")))
  t3 <- get(paste("fold",name,sub(".xls","",file),sep="_"))
  return(abs((2^-(t1-t2)+2^-(t1+t2))/2-t3))
}

TT <- function(name,file){
  if (sub(".xls","",file)=="16s")
    return(1)
  
  t1 <- get(paste("power",name,sub(".xls","",file),sep="_"))
  t2 <- get(paste("power_t9801",sub(".xls","",file),sep="_"))
  if (var.test(t1,t2)$p.value >= 0.05)
    return(t.test(t1,t2,var.equal=T)$p.value)
  else
    return(t.test(t1,t2)$p.value)
}

# 2-7目标减内参
for (i in seq_along(filenames))
  assign(paste("t2_7",sub(".xls","",filenames[i]),sep="_"),
         get(paste("Strain","2_7",sub(".xls","",filenames[i]),sep="_"))-
           get(paste("Strain","2_7","16s",sep="_")))
# 13-2目标减内参
for (i in seq_along(filenames))
  assign(paste("t13_2",sub(".xls","",filenames[i]),sep="_"),
         get(paste("Strain","13_2",sub(".xls","",filenames[i]),sep="_"))-
           get(paste("Strain","13_2","16s",sep="_")))

# 9801目标减内参
for (i in seq_along(filenames))
  assign(paste("t9801",sub(".xls","",filenames[i]),sep="_"),
         get(paste("Strain","9801",sub(".xls","",filenames[i]),sep="_"))-
           get(paste("Strain","9801","16s",sep="_")))

# 2^(2-7目标减内参)
for (i in seq_along(filenames))
  assign(paste("power_t2_7",sub(".xls","",filenames[i]),sep="_"),
         2^-get(paste("t2_7",sub(".xls","",filenames[i]),sep="_")))

# 2^(13-2目标减内参)
for (i in seq_along(filenames))
  assign(paste("power_t13_2",sub(".xls","",filenames[i]),sep="_"),
         2^-get(paste("t13_2",sub(".xls","",filenames[i]),sep="_")))

# 2^(9801目标减内参)
for (i in seq_along(filenames))
  assign(paste("power_t9801",sub(".xls","",filenames[i]),sep="_"),
         2^-get(paste("t9801",sub(".xls","",filenames[i]),sep="_")))

#======

# 2-7表达倍数
for (i in seq_along(filenames))
  assign(paste("fold_t2_7",sub(".xls","",filenames[i]),sep="_"),
         2^-(mean(get(paste("t2_7",sub(".xls","",filenames[i]),sep="_")))-
               mean(get(paste("t9801",sub(".xls","",filenames[i]),sep="_")))))

# 13-2表达倍数
for (i in seq_along(filenames))
  assign(paste("fold_t13_2",sub(".xls","",filenames[i]),sep="_"),
         2^-(mean(get(paste("t13_2",sub(".xls","",filenames[i]),sep="_")))-
               mean(get(paste("t9801",sub(".xls","",filenames[i]),sep="_")))))

# 9801表达倍数
for (i in seq_along(filenames))
  assign(paste("fold_t9801",sub(".xls","",filenames[i]),sep="_"),
         2^-(mean(get(paste("t9801",sub(".xls","",filenames[i]),sep="_")))-
               mean(get(paste("t9801",sub(".xls","",filenames[i]),sep="_")))))

# 2-7表达SD
for (i in seq_along(filenames))
  assign(paste("SD_t2_7",sub(".xls","",filenames[i]),sep="_"),
         SD("t2_7",filenames[i]))

# 13-2表达SD
for (i in seq_along(filenames))
  assign(paste("SD_t13_2",sub(".xls","",filenames[i]),sep="_"),
         SD("t13_2",filenames[i]))

# 9801表达SD
for (i in seq_along(filenames))
  assign(paste("SD_t9801",sub(".xls","",filenames[i]),sep="_"),
         SD("t9801",filenames[i]))

# 2-7表达p-value
for (i in seq_along(filenames))
  assign(paste("p_t2_7",sub(".xls","",filenames[i]),sep="_"),
         TT("t2_7",filenames[i]))

# 13-2表达p-value
for (i in seq_along(filenames))
  assign(paste("p_t13_2",sub(".xls","",filenames[i]),sep="_"),
         TT("t13_2",filenames[i]))

# 9801表达p-value
for (i in seq_along(filenames))
  assign(paste("p_t9801",sub(".xls","",filenames[i]),sep="_"),
         TT("t9801",filenames[i]))

#======

# 数据输出
for (i in seq_along(filenames)){
  d <- get(paste("fold_t2_7",sub(".xls","",filenames[i]),sep="_"))
  d <- as.data.frame(d)
  colnames(d) <- sub(".xls","",filenames[i])
  writeWorksheetToFile("conclusion/fold_2_7.xlsx",
                       data=d,sheet="Sheet1",startRow=1,startCol=i)
}
for (i in seq_along(filenames)){
  d <- get(paste("fold_t13_2",sub(".xls","",filenames[i]),sep="_"))
  d <- as.data.frame(d)
  colnames(d) <- sub(".xls","",filenames[i])
  writeWorksheetToFile("conclusion/fold_13_2.xlsx",
                       data=d,sheet="Sheet1",startRow=1,startCol=i)
}
for (i in seq_along(filenames)){
  d <- get(paste("SD_t2_7",sub(".xls","",filenames[i]),sep="_"))
  d <- as.data.frame(d)
  colnames(d) <- sub(".xls","",filenames[i])
  writeWorksheetToFile("conclusion/SD_2_7.xlsx",
                       data=d,sheet="Sheet1",startRow=1,startCol=i)
}
for (i in seq_along(filenames)){
  d <- get(paste("SD_t13_2",sub(".xls","",filenames[i]),sep="_"))
  d <- as.data.frame(d)
  colnames(d) <- sub(".xls","",filenames[i])
  writeWorksheetToFile("conclusion/SD_13_2.xlsx",
                       data=d,sheet="Sheet1",startRow=1,startCol=i)
}
for (i in seq_along(filenames)){
  d <- get(paste("p_t2_7",sub(".xls","",filenames[i]),sep="_"))
  d <- as.data.frame(d)
  colnames(d) <- sub(".xls","",filenames[i])
  writeWorksheetToFile("conclusion/pvalue_2_7.xlsx",
                       data=d,sheet="Sheet1",startRow=1,startCol=i)
}
for (i in seq_along(filenames)){
  d <- get(paste("p_t13_2",sub(".xls","",filenames[i]),sep="_"))
  d <- as.data.frame(d)
  colnames(d) <- sub(".xls","",filenames[i])
  writeWorksheetToFile("conclusion/pvalue_13_2.xlsx",
                       data=d,sheet="Sheet1",startRow=1,startCol=i)
}

#======

# 作图
