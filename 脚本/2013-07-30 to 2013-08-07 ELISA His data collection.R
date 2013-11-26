library(lattice)
data <- read.csv("/home/kaji331/Projects/R/数据/2013-07-30 to 2013-08-07 ELISA His data collection.csv")
xyplot(value[1:162]~volume[1:162]|group[1:162], data=data)
xyplot(value[1:162]~volume[1:162]|as.factor(group[1:162]), data=data)