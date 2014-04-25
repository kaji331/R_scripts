#R脚本

平时写的各种各样乱七八糟的R脚本

##git的使用
1. 基本上是git init，然后git config --global user.name "\*\*\*"，在然后git config 
--global user.email "\*\*\*"就可以git add .和git commit -m "\*\*\*"了。如果要送到
github.com，就设置git remote add origin https://github.com/kaji331/$$$.git，其中
$$$是要在github.com先创立的，然后git push origin master就可以了。

2. 要克隆就使用git clone https://github.com/kaji331/\*\*\*.git，必须使用https，然后
同样设置git config --global user.name和user.email，接着就可以用了。

3. 用github.com更新本地仓库使用git pull origin master就可以了。

4. 显示所有远程仓库git remote -v，重命名git remote rename 原名 新名，删除远程仓
库git remote rm \*\*\*。

> http://bqnw.me/post/dplyr-note#11-filter
