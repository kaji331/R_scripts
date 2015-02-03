#dynamic 3D surface

x <- seq(-4,4,0.1)
y <- seq(-10,10,0.1)
sc <- function(x,y) {
  sin(x)*cos(sqrt(y))
}
z <- outer(x,y,sc)

library(rgl)
persp3d(x,y,z,col=rainbow(max(nrow(z),ncol(z))))