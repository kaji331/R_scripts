#example from stackoverflow.com 
#"Barplot with significant differences and interactions"
library(gplots)
hh <- t(VADeaths)[1:2, 5:1]
mybarcol <- "gray20"
ci.l <- hh * 0.85
ci.u <- hh * 1.15
mp <- barplot2(hh, beside = TRUE,
               col = c("grey12", "grey82"),
               legend = colnames(VADeaths)[1:2], ylim = c(0, 100),
               cex.names = 1.5, plot.ci = TRUE, ci.l = ci.l, ci.u = ci.u)
y.cord<-rbind(c(ci.u[1,]+1),c(apply(ci.u,2,max)+5),
              c(apply(ci.u,2,max)+5),c(ci.u[2,]+1))
x.cord<-apply(mp,2,function(x) rep(x,each=2))
sapply(1:5,function(x) lines(x.cord[,x],y.cord[,x]))
x.text<-colMeans(mp)
y.text<-apply(ci.u,2,max)+7
text(c("*","**","***","NS","***"),x=x.text,y=y.text)

#example by kaji331
#===========未分组数据============
set.seed(1984)
category <- gl(3,3,labels=c("a","b","c"))
value <- c(runif(3,min=0.8,max=1.2),rnorm(3,mean=0.2,sd=0.1),
           rnorm(3,mean=1.5,sd=0.1))
data <- data.frame(category,value)

library(sciplot)
p <- bargraph.CI(category,value,data=data,col=gray(seq(0.1,0.9,length=3)),
                 ylim=c(0,2),lc=F,family="Times")
box()

#根据每个图来微调参数，然后话连接线
y.cord <- c(p$CI[2,1,1]+0.1,max(p$CI[2,1,1:2])+0.2,max(p$CI[2,1,1:2])+0.2,
            p$CI[2,1,2]+0.1)
x.cord <- c(rep(p$xvals[1,1],each=2),rep(p$xvals[2,1],each=2))
lines(x.cord,y.cord)

#画星星
x.text <- mean(p$xvals[1:2,1])
y.text <- max(p$CI[2,1,1:2])+0.3
text("**",x=x.text,y=y.text)

#===========分组数据============
library(sciplot)
p <- bargraph.CI(dose,len,group=supp,data=ToothGrowth,xlab="Dose",ylab="Growth",
                 cex.lab=1.5,x.leg=1,col="black",angle=45,cex.names=1.25,
                 density=c(0,20),legend=T,lc=F,ylim=c(0,30),family="Times")
box()

#根据每个图来微调参数，然后画连接线
y.cord <- matrix(nrow=0,ncol=4)
for (i in 1:2)
  y.cord <-rbind(y.cord,c(p$CI[2,1,i]+1,max(p$CI[2,1:2,i])+2,max(p$CI[2,1:2,i])+2,
             p$CI[2,2,i]+1))

x.cord <- matrix(nrow=0,ncol=4)
for (i in 1:2)
  x.cord <- rbind(x.cord,c(rep(p$xvals[1,i],each=2),rep(p$xvals[2,i],each=2)))
for (i in 1:2)
  lines(x.cord[i,],y.cord[i,])

#画星星
x.text <- colMeans(p$xvals[,1:2])
y.text <- apply(p$CI[2,1:2,],2,max)+3
text("**",x=x.text,y=y.text)