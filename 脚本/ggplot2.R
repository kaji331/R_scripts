library(ggplot2)
##
data <- read.csv("/home/kaji331/Projects/R/R_scripts/数据/data.csv")
P <- ggplot(data, aes(x = name, y = value, color = name))
P + geom_boxplot() + geom_density(stat = "identity")
##
d <- data.frame(tapply(data$value, data$name, mean)[], names(tapply(data$value, data$name, mean)))
names(d) <- c("concentration", "protein")
p <- ggplot(d, aes(x = protein, y = concentration, color = protein, fill = protein))
p + geom_bar(alpha = 0.5, color = "black", stat = "identity") + ggtitle("Bradford Assay of E2")
plot(density(d$concentration))
##
data2 <- read.csv("/home/kaji331/Projects/R/R_scripts/数据/2013-06-03 ELISA.csv")
#
a <- data.frame(data2$X467[1:8], data2$X.3[1:8])
names(a) <- c("value", "name")
A <- ggplot(a, aes(x = name, y = value, group = 1))
A + geom_point(color = "red", size = 2, shape = 10, fill = "white") + geom_line(color = "blue", linetype = "dotted", size = 1.5)
#
b <- data.frame(data2$X467[19:26], data2$X.3[19:26])
names(b) <- c("value", "name")
B <- ggplot(b, aes(x = name, y = value, color = name, group = 1))
B + geom_point() + geom_line()
#
c <- data.frame(data2$X467[37:44], data2$X.3[37:44])
names(c) <- c("value", "name")
C <- ggplot(c, aes(x = name, y = value, color = name, group = 1))
C + geom_point() + geom_line()