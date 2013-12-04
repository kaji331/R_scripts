R中快速方便的数据结构调整
========================================================

R中需要将数据整合成长列表的形式才能进行很多统计分析，但是数据来自各种检测仪器或
Excel以及为了方便书写记录，往往都以直观简单的横列表或纵列表形式存在。为了除去在
Excel等其他工具中手动调整数据结构的麻烦，R其实提供了非常强大的工具和方法，本文以
简单例子来说明。

## 生成例子横列表

```r
set.seed(1984)
a <- runif(10, min = 1, max = 1.5)
b <- rbinom(10, size = 1, prob = 0.5)
c <- rnorm(10, mean = 1, sd = 0.5)
data <- rbind(a, b, c)
```


<!-- html table generated in R 3.0.2 by xtable 1.7-1 package -->
<!-- Wed Dec  4 17:05:28 2013 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> 1 </TH> <TH> 2 </TH> <TH> 3 </TH> <TH> 4 </TH> <TH> 5 </TH> <TH> 6 </TH> <TH> 7 </TH> <TH> 8 </TH> <TH> 9 </TH> <TH> 10 </TH>  </TR>
  <TR> <TD align="right"> a </TD> <TD align="right"> 1.33 </TD> <TD align="right"> 1.22 </TD> <TD align="right"> 1.19 </TD> <TD align="right"> 1.17 </TD> <TD align="right"> 1.37 </TD> <TD align="right"> 1.43 </TD> <TD align="right"> 1.02 </TD> <TD align="right"> 1.22 </TD> <TD align="right"> 1.41 </TD> <TD align="right"> 1.11 </TD> </TR>
  <TR> <TD align="right"> b </TD> <TD align="right"> 1.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 1.00 </TD> <TD align="right"> 1.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 1.00 </TD> <TD align="right"> 1.00 </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> c </TD> <TD align="right"> 1.14 </TD> <TD align="right"> 1.14 </TD> <TD align="right"> 1.42 </TD> <TD align="right"> 1.51 </TD> <TD align="right"> 1.01 </TD> <TD align="right"> 1.86 </TD> <TD align="right"> 1.06 </TD> <TD align="right"> 0.11 </TD> <TD align="right"> 1.06 </TD> <TD align="right"> 0.14 </TD> </TR>
   </TABLE>


## 转置成纵列表

```r
data_t <- t(data)
```


<!-- html table generated in R 3.0.2 by xtable 1.7-1 package -->
<!-- Wed Dec  4 17:05:28 2013 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> a </TH> <TH> b </TH> <TH> c </TH>  </TR>
  <TR> <TD align="right"> 1 </TD> <TD align="right"> 1.33 </TD> <TD align="right"> 1.00 </TD> <TD align="right"> 1.14 </TD> </TR>
  <TR> <TD align="right"> 2 </TD> <TD align="right"> 1.22 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 1.14 </TD> </TR>
  <TR> <TD align="right"> 3 </TD> <TD align="right"> 1.19 </TD> <TD align="right"> 1.00 </TD> <TD align="right"> 1.42 </TD> </TR>
  <TR> <TD align="right"> 4 </TD> <TD align="right"> 1.17 </TD> <TD align="right"> 1.00 </TD> <TD align="right"> 1.51 </TD> </TR>
  <TR> <TD align="right"> 5 </TD> <TD align="right"> 1.37 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 1.01 </TD> </TR>
  <TR> <TD align="right"> 6 </TD> <TD align="right"> 1.43 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 1.86 </TD> </TR>
  <TR> <TD align="right"> 7 </TD> <TD align="right"> 1.02 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 1.06 </TD> </TR>
  <TR> <TD align="right"> 8 </TD> <TD align="right"> 1.22 </TD> <TD align="right"> 1.00 </TD> <TD align="right"> 0.11 </TD> </TR>
  <TR> <TD align="right"> 9 </TD> <TD align="right"> 1.41 </TD> <TD align="right"> 1.00 </TD> <TD align="right"> 1.06 </TD> </TR>
  <TR> <TD align="right"> 10 </TD> <TD align="right"> 1.11 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.14 </TD> </TR>
   </TABLE>


## 如果首列为变量名可以用类似下列命令来转置后再重命名

```r
data_t <- t(data[, 2:10])
names(data_t) <- c("a", "b", "c")
```


<!-- html table generated in R 3.0.2 by xtable 1.7-1 package -->
<!-- Wed Dec  4 17:05:28 2013 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> a </TH> <TH> b </TH> <TH> c </TH>  </TR>
  <TR> <TD align="right"> 1 </TD> <TD align="right"> 1.22 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 1.14 </TD> </TR>
  <TR> <TD align="right"> 2 </TD> <TD align="right"> 1.19 </TD> <TD align="right"> 1.00 </TD> <TD align="right"> 1.42 </TD> </TR>
  <TR> <TD align="right"> 3 </TD> <TD align="right"> 1.17 </TD> <TD align="right"> 1.00 </TD> <TD align="right"> 1.51 </TD> </TR>
  <TR> <TD align="right"> 4 </TD> <TD align="right"> 1.37 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 1.01 </TD> </TR>
  <TR> <TD align="right"> 5 </TD> <TD align="right"> 1.43 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 1.86 </TD> </TR>
  <TR> <TD align="right"> 6 </TD> <TD align="right"> 1.02 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 1.06 </TD> </TR>
  <TR> <TD align="right"> 7 </TD> <TD align="right"> 1.22 </TD> <TD align="right"> 1.00 </TD> <TD align="right"> 0.11 </TD> </TR>
  <TR> <TD align="right"> 8 </TD> <TD align="right"> 1.41 </TD> <TD align="right"> 1.00 </TD> <TD align="right"> 1.06 </TD> </TR>
  <TR> <TD align="right"> 9 </TD> <TD align="right"> 1.11 </TD> <TD align="right"> 0.00 </TD> <TD align="right"> 0.14 </TD> </TR>
   </TABLE>


## *融合成真正需要的长列表形式*

```r
library(reshape2)
data_t <- melt(data_t, value.name = "OD等")
names(data_t)[2] <- "Protein等"
```


<!-- html table generated in R 3.0.2 by xtable 1.7-1 package -->
<!-- Wed Dec  4 17:05:28 2013 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> Var1 </TH> <TH> Protein等 </TH> <TH> OD等 </TH>  </TR>
  <TR> <TD align="right"> 1 </TD> <TD align="right">   1 </TD> <TD> a </TD> <TD align="right"> 1.22 </TD> </TR>
  <TR> <TD align="right"> 2 </TD> <TD align="right">   2 </TD> <TD> a </TD> <TD align="right"> 1.19 </TD> </TR>
  <TR> <TD align="right"> 3 </TD> <TD align="right">   3 </TD> <TD> a </TD> <TD align="right"> 1.17 </TD> </TR>
  <TR> <TD align="right"> 4 </TD> <TD align="right">   4 </TD> <TD> a </TD> <TD align="right"> 1.37 </TD> </TR>
  <TR> <TD align="right"> 5 </TD> <TD align="right">   5 </TD> <TD> a </TD> <TD align="right"> 1.43 </TD> </TR>
  <TR> <TD align="right"> 6 </TD> <TD align="right">   6 </TD> <TD> a </TD> <TD align="right"> 1.02 </TD> </TR>
  <TR> <TD align="right"> 7 </TD> <TD align="right">   7 </TD> <TD> a </TD> <TD align="right"> 1.22 </TD> </TR>
  <TR> <TD align="right"> 8 </TD> <TD align="right">   8 </TD> <TD> a </TD> <TD align="right"> 1.41 </TD> </TR>
  <TR> <TD align="right"> 9 </TD> <TD align="right">   9 </TD> <TD> a </TD> <TD align="right"> 1.11 </TD> </TR>
  <TR> <TD align="right"> 10 </TD> <TD align="right">   1 </TD> <TD> b </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> 11 </TD> <TD align="right">   2 </TD> <TD> b </TD> <TD align="right"> 1.00 </TD> </TR>
  <TR> <TD align="right"> 12 </TD> <TD align="right">   3 </TD> <TD> b </TD> <TD align="right"> 1.00 </TD> </TR>
  <TR> <TD align="right"> 13 </TD> <TD align="right">   4 </TD> <TD> b </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> 14 </TD> <TD align="right">   5 </TD> <TD> b </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> 15 </TD> <TD align="right">   6 </TD> <TD> b </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> 16 </TD> <TD align="right">   7 </TD> <TD> b </TD> <TD align="right"> 1.00 </TD> </TR>
  <TR> <TD align="right"> 17 </TD> <TD align="right">   8 </TD> <TD> b </TD> <TD align="right"> 1.00 </TD> </TR>
  <TR> <TD align="right"> 18 </TD> <TD align="right">   9 </TD> <TD> b </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> 19 </TD> <TD align="right">   1 </TD> <TD> c </TD> <TD align="right"> 1.14 </TD> </TR>
  <TR> <TD align="right"> 20 </TD> <TD align="right">   2 </TD> <TD> c </TD> <TD align="right"> 1.42 </TD> </TR>
  <TR> <TD align="right"> 21 </TD> <TD align="right">   3 </TD> <TD> c </TD> <TD align="right"> 1.51 </TD> </TR>
  <TR> <TD align="right"> 22 </TD> <TD align="right">   4 </TD> <TD> c </TD> <TD align="right"> 1.01 </TD> </TR>
  <TR> <TD align="right"> 23 </TD> <TD align="right">   5 </TD> <TD> c </TD> <TD align="right"> 1.86 </TD> </TR>
  <TR> <TD align="right"> 24 </TD> <TD align="right">   6 </TD> <TD> c </TD> <TD align="right"> 1.06 </TD> </TR>
  <TR> <TD align="right"> 25 </TD> <TD align="right">   7 </TD> <TD> c </TD> <TD align="right"> 0.11 </TD> </TR>
  <TR> <TD align="right"> 26 </TD> <TD align="right">   8 </TD> <TD> c </TD> <TD align="right"> 1.06 </TD> </TR>
  <TR> <TD align="right"> 27 </TD> <TD align="right">   9 </TD> <TD> c </TD> <TD align="right"> 0.14 </TD> </TR>
   </TABLE>


## 后续操作等，如……

```r
data <- data_t[, 2:3]
```


<!-- html table generated in R 3.0.2 by xtable 1.7-1 package -->
<!-- Wed Dec  4 17:05:28 2013 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> Protein等 </TH> <TH> OD等 </TH>  </TR>
  <TR> <TD align="right"> 1 </TD> <TD> a </TD> <TD align="right"> 1.22 </TD> </TR>
  <TR> <TD align="right"> 2 </TD> <TD> a </TD> <TD align="right"> 1.19 </TD> </TR>
  <TR> <TD align="right"> 3 </TD> <TD> a </TD> <TD align="right"> 1.17 </TD> </TR>
  <TR> <TD align="right"> 4 </TD> <TD> a </TD> <TD align="right"> 1.37 </TD> </TR>
  <TR> <TD align="right"> 5 </TD> <TD> a </TD> <TD align="right"> 1.43 </TD> </TR>
  <TR> <TD align="right"> 6 </TD> <TD> a </TD> <TD align="right"> 1.02 </TD> </TR>
  <TR> <TD align="right"> 7 </TD> <TD> a </TD> <TD align="right"> 1.22 </TD> </TR>
  <TR> <TD align="right"> 8 </TD> <TD> a </TD> <TD align="right"> 1.41 </TD> </TR>
  <TR> <TD align="right"> 9 </TD> <TD> a </TD> <TD align="right"> 1.11 </TD> </TR>
  <TR> <TD align="right"> 10 </TD> <TD> b </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> 11 </TD> <TD> b </TD> <TD align="right"> 1.00 </TD> </TR>
  <TR> <TD align="right"> 12 </TD> <TD> b </TD> <TD align="right"> 1.00 </TD> </TR>
  <TR> <TD align="right"> 13 </TD> <TD> b </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> 14 </TD> <TD> b </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> 15 </TD> <TD> b </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> 16 </TD> <TD> b </TD> <TD align="right"> 1.00 </TD> </TR>
  <TR> <TD align="right"> 17 </TD> <TD> b </TD> <TD align="right"> 1.00 </TD> </TR>
  <TR> <TD align="right"> 18 </TD> <TD> b </TD> <TD align="right"> 0.00 </TD> </TR>
  <TR> <TD align="right"> 19 </TD> <TD> c </TD> <TD align="right"> 1.14 </TD> </TR>
  <TR> <TD align="right"> 20 </TD> <TD> c </TD> <TD align="right"> 1.42 </TD> </TR>
  <TR> <TD align="right"> 21 </TD> <TD> c </TD> <TD align="right"> 1.51 </TD> </TR>
  <TR> <TD align="right"> 22 </TD> <TD> c </TD> <TD align="right"> 1.01 </TD> </TR>
  <TR> <TD align="right"> 23 </TD> <TD> c </TD> <TD align="right"> 1.86 </TD> </TR>
  <TR> <TD align="right"> 24 </TD> <TD> c </TD> <TD align="right"> 1.06 </TD> </TR>
  <TR> <TD align="right"> 25 </TD> <TD> c </TD> <TD align="right"> 0.11 </TD> </TR>
  <TR> <TD align="right"> 26 </TD> <TD> c </TD> <TD align="right"> 1.06 </TD> </TR>
  <TR> <TD align="right"> 27 </TD> <TD> c </TD> <TD align="right"> 0.14 </TD> </TR>
   </TABLE>

