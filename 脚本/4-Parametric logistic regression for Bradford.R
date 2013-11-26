library(drc)
##OD已经平均
data <- read.csv("/home/kaji331/Projects/R/数据/Bradford-standard.csv")
p <- drm(Mean.OD.Value~Concentration, data = data, fct = 
           l4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
summary(p)
drc:::Rsq(p)
#R^2不适合于非线性回归，应该用相对标准离差（相对标准差）RSD=SD/Mean
abs(sqrt(summary(p)$"resVar") / mean(fitted(p)))
ED(p, c(5,50,95), interval="delta")
plot(p)

##OD没有平均
data2 <- read.csv("/home/kaji331/Projects/R/数据/Bradford-standard-2.csv")
p2 <- drm(OD.Values~Concentration, data = data2, fct = 
           l4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
summary(p2)
ED(p2, c(5,50,95), interval="delta")
attributes(p2)
#可以做方差检验
modelFit(p2)
abs(sqrt(summary(p2)$"resVar") / mean(fitted(p2)))
plot(p2)
#通过浓度预测OD
Conc <- c(seq(0,10,0.1))
newdata <- data.frame(Conc)
p3 <- predict(p2, newdata, type="response")
plot(p3~Conc)
#通过OD预测浓度
p4 <- drm(Concentration~OD.Values, data = data2, fct = 
            l4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
OD <- c(seq(0.4,1.6,0.02))
newdata <- data.frame(OD)
p5 <- predict(p4, newdata, type="response")
plot(p5~OD)

##The four-parameter Weibull functions
data3 <- read.csv("/home/kaji331/Projects/R/数据/Bradford-standard-2.csv")
p6 <- drm(OD.Values~Concentration, data = data2, fct = 
            w4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
summary(p6)
ED(p6, c(5,50,95), interval="delta")
attributes(p6)
#可以做方差检验
modelFit(p6)
abs(sqrt(summary(p6)$"resVar") / mean(fitted(p6)))
plot(p6,family="Times New Roman",font.lab=2)
axis(4,at=p6$dataList$resp,col.axis="blue",las=2,family="sans")