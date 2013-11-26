library(ggplot2)

data <- read.csv("/home/kaji331/Projects/R/数据/2013-06-10 ELISA His.csv")

ggplot(data, aes(x=Group, y=OD, color=Group)) + geom_boxplot() + geom_jitter()

data.aov <- aov(OD~Group, data=data)
summary(data.aov)

pairwise.t.test(data$OD, data$Group, p.adj="none")
pairwise.t.test(data$OD, data$Group, p.adj="bonferroni") #more conservative
pairwise.t.test(data$OD, data$Group, p.adj="holm") #generally considered superior to the Bonferroni
TukeyHSD(data.aov)