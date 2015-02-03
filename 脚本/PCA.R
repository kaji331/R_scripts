# 根据系统生物学的网络分析课程中Matlab程序，总结的R进行PCA分析的方法，但最后出图不大对

data <- read.table("~/Projects/R/R_scripts/数据/organs_repeats.csv",sep="\t",header=T)
expressions <- data[,2:10] # 这里忽略了基因名，因为有重复不能作为行名，需要进行去重

pc <- prcomp(t(expressions)) # 这个函数结果与Matlab的较接近

# 不知道下面应该是选择x还是应该rotation（似乎应该是rotation）

# 这个是x的版本
x <- scale(pc$x[,1]) # R中scale就是Matlab的zscore
y <- scale(pc$x[,2])
z <- scale(pc$x[,3])

library(rgl)
plot3d(x,y,z) # 点太少

# 这个是rotation的版本
x <- scale(pc$rotation[,1]) # R中scale就是Matlab的zscore
y <- scale(pc$rotation[,2])
z <- scale(pc$rotation[,3])

library(rgl)
plot3d(x,y,z) # 点太多，没有清晰分成三块