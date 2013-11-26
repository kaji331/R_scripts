library(drc)
head(S.alba)
S.alba.m1 <- drm(DryMatter~Dose,Herbicide,
                 fct=l4(names=c("Slope","Lower Limits","Upper Limits","ED50")),data=S.alba)
summary(S.alba.m1) #可以看到上下限比较接近
modelFit(S.alba.m1) #评估曲线拟合质量
S.alba.m2 <- drm(DryMatter~Dose,Herbicide,
                 fct=l4(names=c("Slope","Lower Limits","Upper Limits","ED50")),
                 pmodels=data.frame(Herbicide,1,1,Herbicide),data=S.alba) #归一化上下限
anova(S.alba.m1,S.alba.m2) #比较归一化后与原始回归是否无显著差异
SI(S.alba.m2,c(50,50)) #比较归一化后两个曲线ED50差距
SI(S.alba.m2,c(50,50),reverse=TRUE) #倒数结果
plot(S.alba.m2,broken=TRUE,legend="") #legend有问题，会中断脚本执行，但是仍然可以去除自动图例
plot(S.alba.m1,broken=TRUE,col="red",lty=c(3,4),add=TRUE,legend="")
relpot(S.alba.m2, interval = "delta") #相对功效