#读取抗体分泌原始数据
data_antibody <- read.csv(
  "/home/kaji331/Projects/R/数据/ELISPOT多次实验数据mean预处理-antibody.csv")

#计算综合得分
z <- scale(data_antibody[,2:3])
score <- apply(z,1,mean)
data_antibody <- cbind(data_antibody,score)
data_antibody$cs <- data_antibody$score+(1-data_antibody$score[1])

#sd_z <- scale(data_antibody[,4:5])
#sd_score <- apply(sd_z,1,mean)
#data_antibody <- cbind(data_antibody,sd_score)
#data_antibody$sd_cs <- data_antibody$sd_score+(1-data_antibody$sd_score[1])

#se_z <- scale(data_antibody[,6:7])
#se_score <- apply(se_z,1,mean)
#data_antibody <- cbind(data_antibody,se_score)
#data_antibody$se_cs <- data_antibody$se_score+(1-data_antibody$se_score[1])

#输出到新文件
names(data_antibody) <- c("Protein","1","2","SD1","SD2","SE2","SE2","Score","Adjusted_score")
write.csv(data_antibody,
          "/home/kaji331/Projects/R/数据/ELISPOT多次实验数据mean预处理-antibody score.csv")

#================卖萌的分割线==================

#读取IFN-γ分泌原始数据
data_IFN <- read.csv(
  "/home/kaji331/Projects/R/数据/ELISPOT多次实验数据mean预处理-IFN-γ.csv")

#计算综合得分
z <- scale(data_IFN[,2:5])
score <- apply(z,1,mean)
data_IFN <- cbind(data_IFN,score)
data_IFN$cs <- data_IFN$score+(3-data_IFN$score[1])

#输出到新文件
names(data_IFN) <- c("Protein","1","2","3","4","SD1","SD2","SD3","SD4","SE1",
                     "SE2","SE3","SE4","Score","Adjusted_score")
write.csv(data_IFN,
          "/home/kaji331/Projects/R/数据/ELISPOT多次实验数据mean预处理-IFN-γ score.csv")