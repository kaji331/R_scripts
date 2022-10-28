#!Rscript --vanilla

library(lubridate)

x <- Sys.Date()
y <- date("2001-10-29")
print(paste0("晓燕的年龄是17岁零",x-y,"天"))
y <- date("2001-12-29")
print(paste0("我的年龄是17岁零",x-y,"天"))
