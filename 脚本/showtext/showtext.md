showtext package
========================================================

example of showtext package


```r
library(showtext)
```

```
## Loading required package: sysfonts
## Loading fonts...
## Loading fonts finished
```

```r
font.add("hkshaonv", "华康少女体 .ttf")

dev.set()
```

```
## pdf 
##   2
```

```r
showtext.begin()
a <- seq(0, 10, length.out = 1000)
plot(sin(a), type = "l", main = "正弦函数", family = "hkshaonv")
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1.png) 

```r
showtext.end()
dev.off()
```

```
## null device 
##           1
```

