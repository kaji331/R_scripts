#读取抗体分泌原始数据
data_antibody <- read.csv(
  "/home/kaji331/Projects/R/数据/ELISPOT多次实验数据预处理-antibody.csv")

#计算综合得分
z <- scale(data_antibody[,2:3])
score <- apply(z,1,mean)
data_antibody <- cbind(data_antibody,score)

data_antibody$cs <- data_antibody$score+(1-data_antibody$score[1])

#输出到新文件
names(data_antibody) <- c("Protein","1","2","Score","Adjusted_score")
write.csv(data_antibody,
      "/home/kaji331/Projects/R/数据/ELISPOT多次实验数据预处理-antibody score.csv")

#================卖萌的分割线==================

#读取IFN-γ分泌原始数据
data_IFN <- read.csv(
  "/home/kaji331/Projects/R/数据/ELISPOT多次实验数据预处理-IFN-γ.csv")

#计算综合得分
z <- scale(data_IFN[,2:5])
score <- apply(z,1,mean)
data_IFN <- cbind(data_IFN,score)

data_IFN$cs <- data_IFN$score+(3-data_IFN$score[1])

#输出到新文件
names(data_IFN) <- c("Protein","1","2","3","4","Score","Adjusted_score")
write.csv(data_IFN,
          "/home/kaji331/Projects/R/数据/ELISPOT多次实验数据预处理-IFN-γ score.csv")