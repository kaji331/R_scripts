library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)


genome <- BSgenome.Hsapiens.UCSC.hg19
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# After digested by restricted enzyme (pattern), bases of exons in fragments longer than 
# limited length divided by bases of whole exons in one chromsome (chr).
rate_exons <- function(pattern,chr,limit) {
  enzyme <- matchPDict(PDict(pattern),genome[[chr]])
  exons <- exons(txdb) %>>% (.[seqnames(.) == chr])
  strand(exons[strand(exons)=="-"]) <- "+" # Change "-" strand to "+"
  exons <- reduce(exons) # reduce function combined different exons in overlapped regions
  # from same strand
  
  s <- c(enzyme[[1]],enzyme[[2]],enzyme[[3]],enzyme[[4]]) %>>% start %>>% sort
  e <- c(enzyme[[1]],enzyme[[2]],enzyme[[3]],enzyme[[4]]) %>>% end %>>% sort
  
  t <- GRanges(Rle(chr),IRanges(s[1:(length(s)-1)],e[2:length(e)]),
               Rle(rep("+"),length(s)-1))
  t <- t[width(t) >= limit]
  return(intersect(exons,t,ignore.strand=T) %>>% width %>>% (sum(.)/sum(width(exons))))
}

# Whole average rate digested by Ptel, Eco52I, Pdil, and Cfr42I from chr1 to chr22, chrX, 
# chrY.
mean(sapply(names(genome)[1:24],
            function(x) rate_exons(c("GCGCGC","CGGCCG","GCCGGC","CCGCGG"),x,3500)))

# ======
library(rtracklayer)

# Download bed from 
# http://sourceforge.net/projects/cohcap/files/COHCAP_BSSEQ_anno.zip/download
# or from
# http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cpgIslandExt.txt.gz
# to make bed file (second column to fifth column-chrom,start,end,name).
CGI <- import("~/Projects/R/R_scripts/数据/UCSC_CpG_Islands.bed")

# non-coveraged eCGI rate after digested
nc_rate_ecgi <- function(pattern,chr,limit) {
  enzyme <- matchPDict(PDict(pattern),genome[[chr]])
  cgi <- CGI[seqnames(CGI) == chr]
  start(cgi) <- start(cgi)-1000
  end(cgi) <- end(cgi)+1000
  ecgi <- reduce(cgi)
  
  s <- c(enzyme[[1]],enzyme[[2]],enzyme[[3]],enzyme[[4]]) %>>% start %>>% sort
  e <- c(enzyme[[1]],enzyme[[2]],enzyme[[3]],enzyme[[4]]) %>>% end %>>% sort
  
  t <- GRanges(Rle(chr),IRanges(s[1:(length(s)-1)],e[2:length(e)]),
               Rle(rep("+"),length(s)-1))
  t <- t[width(t) >= limit]
  return(intersect(ecgi,t,ignore.strand=T) %>>% width %>>% (1-sum(.)/sum(width(ecgi))))
}

mean(sapply(names(genome)[1:24],
            function(x) nc_rate_ecgi(c("GCGCGC","CGGCCG","GCCGGC","CCGCGG"),x,3500)))

detect_ecgi <- function(pattern,chr) {
  enzyme <- matchPDict(PDict(pattern),genome[[chr]])
  cgi <- CGI[seqnames(CGI) == chr]
  start(cgi) <- start(cgi)-1000
  end(cgi) <- end(cgi)+1000
  ecgi <- reduce(cgi)
  
  s <- c(enzyme[[1]],enzyme[[2]],enzyme[[3]],enzyme[[4]]) %>>% start %>>% sort
  e <- c(enzyme[[1]],enzyme[[2]],enzyme[[3]],enzyme[[4]]) %>>% end %>>% sort

  t <- GRanges(Rle(chr),IRanges(s[1:(length(s)-1)],e[2:length(e)]),
               Rle(rep("+"),length(s)-1))
  result <- countOverlaps(ecgi,t,ignore.strand=T) %>>% table %>>% data.frame
  colnames(result) <- c("Detected","Freq.")
  result$Detected <- as.integer(result$Detected)
  
  return(result)
}

count_ecgi <- function(chr) {
  cgi <- CGI[seqnames(CGI) == chr]
  start(cgi) <- start(cgi)-1000
  end(cgi) <- end(cgi)+1000
  ecgi <- reduce(cgi)
  return(length(ecgi))
}

detected <- rbind(sapply(names(genome)[1:24],
                         function(x) detect_ecgi(c("GCGCGC","CGGCCG","GCCGGC","CCGCGG"),x))) 
count <- sapply(names(genome)[1:24], function(x) count_ecgi(x)) %>>% sum
d_1 <- apply(detected,2,function(x) sum(x$Freq.[which(x$Detected == 1)])) %>>% sum
d_2to5 <- apply(detected,2,function(x) 
  sum(x$Freq.[which(x$Detected >= 2 & x$Detected <= 5)])) %>>% sum
d_6to30 <- apply(detected,2,function(x) 
  sum(x$Freq.[which(x$Detected >= 6 & x$Detected <= 30)])) %>>% sum
d_30p <- apply(detected,2,function(x) sum(x$Freq.[which(x$Detected > 30)])) %>>% sum
d_0 <- (count-d_1-d_2to5-d_6to30-d_30p)

# ======
# You can get TSS from genome.ucsc.edu/cgi-bin/hgTables using parameters:Mammal, Human, 
# Feb. 2009 (GRCh37/hg19), Genes and Gene Predictions, UCSC Genes, knownGene, 
# (output format)selected fields from primary and related tables, (output file)<blank>.
# Or transcripts(txdb) %>>% start.
transcripts <- transcripts(txdb)

detect_promoters <- function(pattern,chr) {
  enzyme <- matchPDict(PDict(pattern),genome[[chr]])
  promoters <- promoters(transcripts[seqnames(transcripts) == chr],
                         upstream=3000,downstream=3000)
  promoters <- reduce(promoters)

  s <- c(enzyme[[1]],enzyme[[2]],enzyme[[3]],enzyme[[4]]) %>>% start %>>% sort
  e <- c(enzyme[[1]],enzyme[[2]],enzyme[[3]],enzyme[[4]]) %>>% end %>>% sort
  
  t <- GRanges(Rle(chr),IRanges(s[1:(length(s)-1)],e[2:length(e)]),
               Rle(rep("+"),length(s)-1))
  result <- countOverlaps(promoters,t,ignore.strand=T) %>>% table %>>% data.frame
  colnames(result) <- c("Detected","Freq.")
  result$Detected <- as.integer(result$Detected)

  return(result)
}

count_promoters <- function(chr) {
  promoters <- promoters(transcripts[seqnames(transcripts) == chr],
                         upstream=3000,downstream=3000)
  promoters <- reduce(promoters)
  return(length(promoters))
}

detected <- rbind(sapply(names(genome)[1:24],
                         function(x) detect_promoters(c("GCGCGC","CGGCCG","GCCGGC","CCGCGG"),x))) 
count <- sapply(names(genome)[1:24], function(x) count_promoters(x)) %>>% sum
d_1 <- apply(detected,2,function(x) sum(x$Freq.[which(x$Detected == 1)])) %>>% sum
d_2to5 <- apply(detected,2,function(x) 
  sum(x$Freq.[which(x$Detected >= 2 & x$Detected <= 5)])) %>>% sum
d_6to30 <- apply(detected,2,function(x) 
  sum(x$Freq.[which(x$Detected >= 6 & x$Detected <= 30)])) %>>% sum
d_30p <- apply(detected,2,function(x) sum(x$Freq.[which(x$Detected > 30)])) %>>% sum
d_0 <- (count-d_1-d_2to5-d_6to30-d_30p)

# ======
# How many cutting sites of one enzyme not in "eCGIs" and "Promoters".
count_no_inside <- function(pattern,chr) {
  enzyme <- matchPattern(pattern,genome[[chr]])
  cgi <- CGI[seqnames(CGI) == chr]
  start(cgi) <- start(cgi)-1000
  end(cgi) <- end(cgi)+1000
  ecgi <- reduce(cgi)
  promoters <- promoters(transcripts[seqnames(transcripts) == chr],
                         upstream=3000,downstream=3000)
  strand(promoters[strand(promoters)=="-"]) <- "+"
  promoters <- reduce(promoters)
  combined_sites <- c(ecgi,promoters)
  strand(combined_sites) <- "+"
  combined_sites <- reduce(combined_sites)
  strand(combined_sites) <- "*"

  s <- start(enzyme)
  e <- end(enzyme)
  
  t <- GRanges(Rle(chr),IRanges(start(enzyme),end(enzyme)))
  # %within% can get ranges which were completely included in larger ranges without specific 
  # strands.
  count_inside <- sapply(combined_sites,function(x) t %within% x) %>>%  sum
  count_no_inside <- length(enzyme)-count_inside
  return(count_no_inside)
}

mclapply(names(genome)[1:24],function(x) length(matchPattern("CCGCGG",genome[[x]])),
         mc.cores=4) %>>%
  simplify2array %>>% sum
mclapply(names(genome)[1:24],function(x) count_no_inside("CCGCGG",x),mc.cores=4) %>>% 
  simplify2array %>>% sum