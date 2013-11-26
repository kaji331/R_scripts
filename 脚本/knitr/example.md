Title
========================================================

This is an R Markdown document. Markdown is a simple formatting syntax for authoring web pages (click the **MD** toolbar button for help on Markdown).

When you click the **Knit HTML** button a web page will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:


```r
library(xtable)
c <- xtable(summary(cars))
print(c, type = "html")
```

<!-- html table generated in R 3.0.2 by xtable 1.7-1 package -->
<!-- Tue Oct 29 18:31:55 2013 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH>     speed </TH> <TH>      dist </TH>  </TR>
  <TR> <TD align="right"> 1 </TD> <TD> Min.   : 4.0   </TD> <TD> Min.   :  2   </TD> </TR>
  <TR> <TD align="right"> 2 </TD> <TD> 1st Qu.:12.0   </TD> <TD> 1st Qu.: 26   </TD> </TR>
  <TR> <TD align="right"> 3 </TD> <TD> Median :15.0   </TD> <TD> Median : 36   </TD> </TR>
  <TR> <TD align="right"> 4 </TD> <TD> Mean   :15.4   </TD> <TD> Mean   : 43   </TD> </TR>
  <TR> <TD align="right"> 5 </TD> <TD> 3rd Qu.:19.0   </TD> <TD> 3rd Qu.: 56   </TD> </TR>
  <TR> <TD align="right"> 6 </TD> <TD> Max.   :25.0   </TD> <TD> Max.   :120   </TD> </TR>
   </TABLE>


You can also embed plots, for example:


```r
plot(cars)
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2.png) 


$$
\mu = \beta_0 + {\beta_1}^2
$$
