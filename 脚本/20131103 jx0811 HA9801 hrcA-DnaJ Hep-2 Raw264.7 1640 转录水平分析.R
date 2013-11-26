data_rnorm <- function(x) { ##x的第一列为均值，第二列为sd，每行数据求6个随机数
  if (!is.data.frame(x))
    stop("invalid data.frame")
  
  d <- data.frame()
  d <- rbind(rnorm(6,x[1,1],x[1,2]))
  for (i in 1:(nrow(x)-1)) {
    d <- rbind(d,rnorm(6,x[i+1,1],x[i+1,2]))
  }
  
  return(d)
}

data <- read.csv("~/Projects/R/数据/20131103 jx0811 HA9801 hrcA-DnaJ Hep-2 Raw264.7 1640 转录水平分析.csv")
data1 <- subset(data,select=c(Mean,SD))
data2 <- data_rnorm(data1)
data2 <- cbind(subset(data,select=c(Bac,Gene,Cell)),data2)
data_Hep_2 <- subset(data2,Cell=="Hep-2")
data_Raw264.7 <- subset(data2,Cell=="Raw264.7")
data_1640 <- subset(data2,Cell=="1640")

d_Hep_2 <- data.frame()
d_Raw264.7 <- data.frame()
d_1640 <- data.frame()

for (i in 1:nrow(data_Hep_2))
  for (j in 4:9)
    d_Hep_2 <- rbind(d_Hep_2,cbind(data_Hep_2[i,1:3],data_Hep_2[i,j]))

for (i in 1:nrow(data_Raw264.7))
  for (j in 4:9)
    d_Raw264.7 <- rbind(d_Raw264.7,cbind(data_Raw264.7[i,1:3],data_Raw264.7[i,j]))

for (i in 1:nrow(data_1640))
  for (j in 4:9)
    d_1640 <- rbind(d_1640,cbind(data_1640[i,1:3],data_1640[i,j]))

names(d_Hep_2) <- c("Bac","Gene","Cell","Value")
names(d_Raw264.7) <- c("Bac","Gene","Cell","Value")
names(d_1640) <- c("Bac","Gene","Cell","Value")