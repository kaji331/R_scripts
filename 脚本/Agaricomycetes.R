# 将一个含10万多条序列的fasta文件中的2万多条伞菌纲序列提取出来成一个新文件用于比对
Agaricomycetes <- readLines("release11_3_Fungi_unaligned.fa")
Agaricomycetes[1]
a <- grep("^>(.*?)Agaricomycetes;class",Agaricomycetes)
b <- grep("^>(.*?)",Agaricomycetes)
c <- match(a,b)
for (i in 1:length(c)) {
  write(Agaricomycetes[b[c[i]]:(b[c[i]+1]-1)],file="rel.fa",append=T,sep="\n")
}
