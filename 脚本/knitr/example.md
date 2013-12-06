Title
========================================================

This is an R Markdown document. Markdown is a simple formatting syntax for authoring web pages (click the **MD** toolbar button for help on Markdown).

When you click the **Knit HTML** button a web page will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:


```r
c <- summary(cars)
kable(c)  #使用knitr自带的kable在Rmd用RStudio输出html时与xtable没有区别；在ReText中kable
```

|id  |    speed       |     dist      |
|:---|:---------------|:--------------|
|    |Min.   : 4.0    |Min.   :  2    |
|    |1st Qu.:12.0    |1st Qu.: 26    |
|    |Median :15.0    |Median : 36    |
|    |Mean   :15.4    |Mean   : 43    |
|    |3rd Qu.:19.0    |3rd Qu.: 56    |
|    |Max.   :25.0    |Max.   :120    |

```r
# 表保持与html一致但xtable为更漂亮的html代码；kable和xtable用ReText都能
# 正常导出且与ReText预览结果一致；kable表格能被pandoc正常转变为pdf的三线
# 表，xtable转变后格式不对；xtable在lyx中为三线表，kable则只是类似Excel
# 的网格线表而已。
library(xtable)
print(xtable(c), type = "html")
```

<!-- html table generated in R 3.0.2 by xtable 1.7-1 package -->
<!-- Fri Dec  6 17:12:07 2013 -->
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


```r
plot(sin(seq(0, 6.2831, length = 60)), type = "b")
```

![plot of chunk 360度sin图](figure/360度sin图.png) 

