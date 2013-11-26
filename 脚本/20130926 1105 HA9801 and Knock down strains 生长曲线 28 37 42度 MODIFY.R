library(drc)
Data <-read.csv("~/Projects/R/数据/20130926 1105 HA9801 and Knock down strains 生长曲线 28 37 42度 MODIFY.csv")
data1 <- subset(Data,Temperature==28)
data1 <- data1[,1:4]
data2 <- subset(Data,Temperature==37)
data2 <- subset(data2,Temperature==37,select=c(Time:Temperature))
data3 <- subset(Data,Temperature==42)
data3 <- subset(data3,Temperature==42,select=c(Time,Names,Data,Temperature))
#data1 <- read.csv("~/Projects/R/数据/20130926 1105 HA9801 and Knock down strains 生长曲线 28 37 42度 MODIFY--28.csv")
#data2 <- read.csv("~/Projects/R/数据/20130926 1105 HA9801 and Knock down strains 生长曲线 28 37 42度 MODIFY--37.csv")
#data3 <- read.csv("~/Projects/R/数据/20130926 1105 HA9801 and Knock down strains 生长曲线 28 37 42度 MODIFY--42.csv")
#data1 <- data.frame(data1$Time,data1$Names,data1$Data)
#data2 <- data.frame(data2$Time,data2$Names,data2$Data)
#data3 <- data.frame(data3$Time,data3$Names,data3$Data)
#names(data1) <- c("Time","Names","Data")
#names(data2) <- c("Time","Names","Data")
#names(data3) <- c("Time","Names","Data")

used <- NA
for (x in levels(data1$Names))
{
  used <- cbind(used,x)
  for (y in levels(data1$Names))
    if (!is.element(y,used))
    {
      d1 <- subset(data1,Names==x)
      d2 <- subset(data1,Names==y)
      data <- rbind(d1,d2)
      p1 <- drm(Data~Time,Names,fct=l4(names=c("slope","lower","upper","ed50")),data=data)
      p2 <- drm(Data~Time,fct=l4(names=c("slope","lower","upper","ed50")),data=data)
      print(c(x,y,"28 degree"))
      a <- anova(p2,p1)$"p value"[2]
      print(a)
      cat('\n\n\n')
    }
}
rm(used,a)

used <- NA
for (x in levels(data2$Names))
{
  used <- cbind(used,x)
  for (y in levels(data2$Names))
    if (!is.element(y,used))
    {
      d1 <- subset(data2,Names==x)
      d2 <- subset(data2,Names==y)
      data <- rbind(d1,d2)
      p1 <- drm(Data~Time,Names,fct=l4(names=c("slope","lower","upper","ed50")),data=data)
      p2 <- drm(Data~Time,fct=l4(names=c("slope","lower","upper","ed50")),data=data)
      print(c(x,y,"37 degree"))
      a <- anova(p2,p1)$"p value"[2]
      print(a)
      cat('\n\n\n')
    }
}
rm(used,a)

used <- NA
for (x in levels(data3$Names))
{
  used <- c(used,x)
  for (y in levels(data3$Names))
    if (!is.element(y,used))
    {
      d1 <- subset(data3,Names==x)
      d2 <- subset(data3,Names==y)
      data <- rbind(d1,d2)
      p1 <- drm(Data~Time,Names,fct=l4(names=c("slope","lower","upper","ed50")),data=data)
      p2 <- drm(Data~Time,fct=l4(names=c("slope","lower","upper","ed50")),data=data)
      print(c(x,y,"42 degree"))
      a <- anova(p2,p1)$"p value"[2]
      print(a)
      cat('\n\n\n')
    }
}
rm(used,a)