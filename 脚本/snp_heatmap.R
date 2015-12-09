library(VariantAnnotation)
library(pheatmap)

fl <- system.file("extdata","chr22.vcf.gz",package="VariantAnnotation")
vcf <- readVcf(fl,"hg19")

temp <- geno(vcf)[[1]]
temp[temp == "0|0"] <- "0"
temp[temp == "1|1"] <- "1"
temp[temp == "0|1" | temp == "1|0"] <- "0.5"

temp <- as.data.frame(temp)
for (i in 1:ncol(temp)) temp[,i] <- as.numeric(temp[,i])
temp <- as.matrix(temp)
temp[temp == 1] <- 0
temp[temp == 2] <- 0.5
temp[temp == 3] <- 1

pheatmap(temp,color=colorRampPalette(c("gray","blue","red"))(3),
		 show_rownames=F,legend_breaks=c(0,0.5,1))
