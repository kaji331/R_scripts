#!Rscript --vanilla

library(lubridate)

x <- Sys.Date()
y <- date("2002-12-29")

print(paste0("你的年龄是18岁零",x-y,"天"))
