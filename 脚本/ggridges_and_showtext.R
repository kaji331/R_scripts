library(cowplot)
library(ggsci)
library(ggridges)
library(showtext)

dev.new(width=6,height=3)
font_add("csm",
		 regular="/usr/share/fonts/truetype/msttcorefonts/Comic_Sans_MS.ttf")

showtext_auto()
ggplot(iris,aes(Sepal.Length,Species,color=Species,fill=Species)) + 
	geom_density_ridges(alpha=0.75) + scale_color_nejm() + scale_fill_nejm() + 
	theme_ridges(font_family="csm",line_size=0.1) + 
	theme(plot.background=element_rect(fill="#F4D7B2",color="#F4D7B2"))
showtext_auto(F)
