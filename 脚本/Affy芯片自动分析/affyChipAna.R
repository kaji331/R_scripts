# source("http://bioconductor.org/biocLite.R")
# biocLite(c("affy","annotate","annaffy","affyQCReport","arrayQualityMetrics","limma","pheatmap"))

dataInput<-function(file="filelist.txt",sep=" ",rnames="filename"){
	require(limma,quietly=TRUE)
	require(affy,quietly=TRUE)
	require(annotate,quietly=TRUE)
	targets<-readTargets(file=file,sep=sep,row.names=rnames) #limma
	data <- ReadAffy(filenames=targets$filename) #affy
	affydb<-annPkgName(data@annotation,type="db")
	list("targets"=targets,"data"=data,"affydb"=affydb)
}
getExpr<-function(data){
	require(affy,quietly=TRUE)
	require(annotate,quietly=TRUE)
	require(annaffy,quietly=TRUE)
	require(data$affydb,character.only = TRUE,quietly=TRUE)
	data$eset<-rma(data$data,verbose=FALSE)
	try({
		APInfo<-mas5calls(data$data,verbose=FALSE)
		present.call<-exprs(APInfo)
		})
	eset.e<-exprs(data$eset)
	if(exists("present.call")){
		data$RPMA<-merge(eset.e,present.call,by="row.names")
		colnames(data$RPMA)<-make.names(gsub("_\\(.*?\\)\\.CEL","",colnames(data$RPMA),ignore.case=TRUE))
		colnames(data$RPMA)<-gsub("\\.y$",".present",colnames(data$RPMA))
		colnames(data$RPMA)<-gsub("\\.x$",".rma",colnames(data$RPMA))
	}else{
		colnames(eset.e)<-paste(colnames(eset.e),"rma",sep=".")
		Row.names=rownames(eset.e)
		data$RPMA<-cbind(Row.names,eset.e)
		data$RPMA<-data.frame(data$RPMA)
	}
	data$symbols<-as.character(aafSymbol(as.character(data$RPMA$Row.names),data$affydb))
	names(data$symbols)<-as.character(data$RPMA$Row.names)
	data$genes<-as.character(aafUniGene(as.character(data$RPMA$Row.names),data$affydb))
	names(data$genes)<-as.character(data$RPMA$Row.names)
	data
}
QCreport<-function(data,fastQCfile="fastQC.pdf",arrayQMfile="QCreport"){
	require(affyQCReport,quietly=TRUE)
	require(arrayQualityMetrics,quietly=TRUE)
	QCReport(data$data,file=fastQCfile)
	arrayQualityMetrics(expressionset = data$eset, outdir = paste(arrayQMfile,"RMA",sep="_"),force=T)
	arrayQualityMetrics(expressionset = data$data, outdir = paste(arrayQMfile,"Raw",sep="_"),force=T,do.logtransform=T)
	eset.e<-exprs(data$eset)
	d<-dist(t(eset.e))
	hc<-hclust(d,method="complete")
	pdf("hclust.complete.pdf")
	plot(hc)
	dev.off()
	cv.y<-numeric()
	y.sd<-apply(eset.e,1,sd)
	y.mean<-apply(eset.e,1,mean)
	cv.y<-y.sd/y.mean
	x<-eset.e
	y2<-cbind(x,cv.y)
	y2<-y2[y2[,ncol(y2)]>0.10,1:(ncol(y2)-1)]
	d<-dist(t(y2))
	hc<-hclust(d,"ave")
	pdf("hclust.ave.pdf")
	plot(hc)
	dev.off()
}
getDesign<-function(data){
	targets<-data$targets
	g<-factor(targets[,2])
	design<-model.matrix(~-1+g)
	if((ncol(targets)>2)){
		fo<-"~-1+g"
		for(i in 3:ncol(targets)){
			assign(paste("v",colnames(targets)[i],sep=""),factor(targets[,i]))
			fo<-paste(fo,paste("v",colnames(targets)[i],sep=""),sep="+")
		}
		design<-model.matrix(as.formula(fo))
		colnames(design)<-gsub("^v","",colnames(design))
	}
	colnames(design)<-gsub("^g","",colnames(design))
	colnames(design)<-gsub("-","_",colnames(design))
	colnames(design)<-make.names(colnames(design))
	design
}
getContrastMatrix<-function(data,design,constr=NULL){
	require(limma,quietly=TRUE)
	if(is.null(constr)){
		g<-factor(data$targets[,2])
		design.col.names<-levels(g)
		if(length(design.col.names)>2){
			for(i in 1:(length(design.col.names)-1)){
				for(j in (i+1):length(design.col.names)){
					constr<-c(constr,paste(design.col.names[i],design.col.names[j],sep="-"))
				}
			}
		}else{
			constr<-paste(design.col.names[2],design.col.names[1],sep="-")
		}
	}
	contrast.matrix<-makeContrasts(contrasts=constr,levels=design)
	contrast.matrix
}
doAnalyze<-function(data,design,contrast.matrix){
	require(limma,quietly=TRUE)
	fit<-lmFit(data$eset,design)
	fit1<-contrasts.fit(fit,contrast.matrix)
	fit2<-eBayes(fit1)
	fit2
}
outputData<-function(data,fit2,contrast.matrix,times=2,cut.p.val=0.01,cut.fdr.val=NULL,outputFolder="output",all=FALSE,scale="row",cluster_rows=TRUE,cluster_cols=TRUE){
	require(limma,quietly=TRUE)
	require(pheatmap,quietly=TRUE)
	dir.create(outputFolder)
	if(file.exists(outputFolder)) setwd(outputFolder)
	results<-NULL
	results.s<-NULL
	constr<-colnames(contrast.matrix)
	len<-length(constr)
	times<-log(times,2)
	symbols<-as.matrix(data$symbols)
	colnames(symbols)<-"symbol"
	genes<-as.matrix(data$genes)
	colnames(genes)<-"gene"
	for(i in 1:len){
		results[[i]]<-topTable(fit2,coef=i,n=nrow(fit2))
		results[[i]]$AveExpr<-NULL
		try({
		if(colnames(results[[i]])[1]!="ID") results[[i]] <- cbind(ID=rownames(results[[i]]), results[[i]])
		index<-grepl("AFFX",results[[i]]$ID)
		results[[i]]<-results[[i]][!index,]
		})
		if(colnames(results[[i]])[1]!="ID") stop("expect first column to be ID")
		results.s[[i]]<-results[[i]]
		colnames(results.s[[i]])[2:6]<-paste(constr[i],colnames(results.s[[i]])[2:6],sep=".")
	}
	merge.results.s<-results.s[[1]]
	if(len>1){
		for(i in 2:len){
			merge.results.s<-merge(results.s[[i]],merge.results.s,by="ID")
		}
	}
	all.results.s<-merge(merge.results.s,data$RPMA,by.x="ID",by.y="Row.names")
	all.results.s<-merge(genes,all.results.s,by.x="row.names",by.y="ID")
	all.results.s<-merge(symbols,all.results.s,by.x="row.names",by.y="Row.names")
	colnames(all.results.s)[1]<-"ID"
	write.csv(all.results.s,"all.results.csv",row.names=FALSE)
	for(i in 1:len){
		tmp<-results[[i]]
		tmp<-merge(tmp,data$RPMA,by.x="ID",by.y="Row.names")
		tmp<-merge(genes,tmp,by.x="row.names",by.y="ID")
		tmp<-merge(symbols,tmp,by.x="row.names",by.y="Row.names")
		colnames(tmp)[1]<-"ID"
		tmp<-tmp[tmp$P.Value<cut.p.val,]
		if(!is.null(cut.fdr.val)) tmp<-tmp[tmp$adj.P.Val<cut.fdr.val,]
		up<-tmp[tmp$logFC>times,]
		lo<-tmp[tmp$logFC<(-times),]
		up<-up[order(-up$logFC),]
		lo<-lo[order(lo$logFC),]
		write.csv(up,paste("up",constr[i],"csv",sep="."),row.names=FALSE)
		write.csv(lo,paste("lo",constr[i],"csv",sep="."),row.names=FALSE)
		if(all) selected<-data$RPMA[data$RPMA$Row.names %in% c(up$ID,lo$ID),]
		else selected<-rbind(up,lo)
		if(nrow(selected)>2){
			try({
			colnames(selected)[1]<-"ID"
			sname<-symbols[selected$ID,"symbol"]
			sname[sname=="character(0)"]<-"---"
			fmt<-paste("%-",max(nchar(sname)),"s",sep="")
			bname<-sprintf(fmt,sname)
			sname<-paste(bname,names(sname),sep=" ")
			rownames(selected)<-sname
			selected<-selected[,grepl(".rma$",colnames(selected))]
			colnames(selected)<-gsub(".rma$","",colnames(selected))
			if(all){
				annotation<-data$targets[,2:ncol(data$targets)]
				rownames(annotation)<-make.names(gsub("_\\(.*?\\)\\.CEL","",rownames(annotation),ignore.case=TRUE))
			}else{
				gp<-unlist(strsplit(constr[i],"-"))
				gp<-gsub("[\\(\\)]","",gp)
				targets<-data$targets[data$targets$factor %in% gp,]
				selcol<-make.names(gsub("_\\(.*?\\)\\.CEL","",rownames(targets),ignore.case=TRUE))
				selected<-selected[,selcol]
				annotation<-targets[,2:ncol(targets),drop=FALSE]
				rownames(annotation)<-selcol
			}
			pdf(paste("heatmap",constr[i],"pdf",sep="."),width=8+ncol(selected)/10,height=2+nrow(selected)/10)
			pheatmap(selected, annotation=annotation, scale=scale, fontsize=9, fontsize_row=6, fontfamily="mono", cluster_rows=cluster_rows, cluster_cols=cluster_cols)
			dev.off()
				})
		}
	}
}

doAll<-function(file="filelist.txt",workingdir=".",sep=" ",rnames="filename",fastQCfile="fastQC.pdf",arrayQMfile="QCreport",
constr=NULL,times=2,cut.p.val=0.01,cut.fdr.val=NULL,outputFolder="output",
all=FALSE,scale="row",doQC=TRUE,cluster_rows=TRUE,cluster_cols=TRUE){
	setwd(workingdir)
	cat("Data input\n")
	data<-dataInput(file,sep,rnames)
	cat("Get Expression\n")
	data<-getExpr(data)
	if(doQC){
		cat("Do QC\n")
		QCreport(data,fastQCfile,arrayQMfile)
	}
	cat("Do limma analysis\n")
	design<-getDesign(data)
	contrast.matrix<-getContrastMatrix(data,design,constr)
	fit2<-doAnalyze(data,design,contrast.matrix)
	cat("Data output\n")
	outputData(data,fit2,contrast.matrix,times,cut.p.val,cut.fdr.val,outputFolder,all,scale,cluster_rows,cluster_cols)
	cat("..Done..\n")
}
