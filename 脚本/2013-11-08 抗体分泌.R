data <- read.csv("~/Projects/R/数据/2013-11-08 抗体分泌.csv")

## 数据整理
Value <- data$C
for (i in 2:10)
  Value <- c(Value,data[,i])
Group <- gl(10,8,labels=names(data))
data <- data.frame(Value,Group)

write.csv(data,"~/Projects/R/数据/2013-11-08 抗体分泌_processed.csv",
          row.names=FALSE)

## 数据作图
boxplot(Value~Group,data=data)
library(ggplot2)
m <- ggplot(data,aes(x=Group,y=Value))
m + geom_bar(stat="identity")
m + geom_boxplot()
library(sciplot)
bargraph.CI(Group,Value,data=data,xlab="Proteins",ylab="Mean of counts",
            space=0.1,ci.fun=function(x) c(mean(x,na.rm=TRUE)-sd(x,na.rm=TRUE),
                                           mean(x,na.rm=TRUE)+sd(x,na.rm=TRUE)),
            col=gray(seq(0.1,0.9,length=10)))
bargraph.CI(Group,Value,data=data,xlab="Proteins",ylab="Mean of counts",
            space=0.1,ci.fun=function(x) c(mean(x,na.rm=TRUE)-sd(x,na.rm=TRUE),
                                           mean(x,na.rm=TRUE)+sd(x,na.rm=TRUE)),
            col="black",angle = 45,density=c(0,20))
bargraph.CI(Group,Value,data=data,xlab="Proteins",ylab="Mean of counts",
            space=0.1,ci.fun=function(x) c(mean(x,na.rm=TRUE)-sd(x,na.rm=TRUE),
                                           mean(x,na.rm=TRUE)+sd(x,na.rm=TRUE)),
            col=rainbow(10))
lineplot.CI(Group,Value,data=data,xlab="Proteins",ylab="Mean of counts",
            ci.fun=function(x) c(mean(x,na.rm=TRUE)-sd(x,na.rm=TRUE),
                                 mean(x,na.rm=TRUE)+sd(x,na.rm=TRUE)),
            col=rainbow(10))