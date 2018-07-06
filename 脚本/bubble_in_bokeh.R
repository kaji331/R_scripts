library(ggplot2)
library(cowplot)

set.seed(1984)
d <- data.frame(x=runif(4000,0,1) * 100,y=runif(4000,0,1) * 100,
				radii=runif(4000,0,1) * 5)
colors <- apply(d,1,function(x) 
				sprintf("#%02x%02x%02x",as.integer(50 + 2 * x[1]),
						as.integer(30 + 2 * x[2]),150))
d <- cbind(d,colors)
ggplot(d,aes(x,y,size=as.factor(radii),color=colors)) + geom_point(alpha=0.6) + 
	scale_color_manual(values=levels(as.factor(colors))) + 
	scale_size_manual(values=as.numeric(levels(as.factor(d$radii)))) + 
	theme_bw() + theme(legend.position="none") -> g
show(g)
