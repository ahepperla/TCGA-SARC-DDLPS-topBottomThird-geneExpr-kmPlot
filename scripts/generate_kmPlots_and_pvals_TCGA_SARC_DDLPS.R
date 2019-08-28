###################################################################################################
#Code to generate KM plots and generate log rank p-values given gene expression and survival data
#
#
#Usage: Rscript generate_kmPlots_and_pvalsC_TCGA_SARC_DDLPS.R /list/of/genes.txt /matrix/of/expression/values.txt /path/to/survival/data.txt /path/to/outdir
#***Must have pamr and qvalue packages installed, and /path/to/outdir must already exist
###################################################################################################

#necessary libraries
library(pamr)
library(qvalue)

#define function to plot KM and return p-value
kmPlot<-function(x,event,stime,varName="",ymin=0,lineColors=NA,nclasses=NA,overall=F){ #accept the labels (x), event (Disease specific death or not), number of days overall survival (stime), and gene name (varName)
	
	mainLabel <- varName
	event<-as.numeric(as.vector(t(event)))
	stime<-as.numeric(as.vector(t(stime)))
	if(is.numeric(x)){
		tclass <- factor(x[!is.na(x)])
		event <- event[!is.na(x)]
		stime <- stime[!is.na(x)]
	}else{
		tclass <- factor(x[x!=""])
		event <- event[x!=""]
		stime <- stime[x!=""]
	}

	if(is.na(nclasses)){
		nclasses<-nlevels(tclass)
	}
	
	if(length(lineColors)<=1){
		lineColors<-seq(1,nclasses)
	}
	
	#fit KM Plot
	y<-survfit(Surv(time=stime, event=event, type="right")~tclass)
	#plot KM
	plot(y,col=lineColors,main=mainLabel,ylim=c(0,1),ylab="DSS Probability",xlab="Days",xlim=c(0, max(stime, na.rm=T)), mark.time=T)
	
	if(overall==T){
		yo<-survfit(Surv(stime, event)~rep(1,length(tclass)))
		lines(yo,col=1,lwd=2)
		legend("bottomleft",legend=c(levels(tclass),"All"),col=c(lineColors,1),bty="n",lty=rep(1,(nclasses+1)),lwd=c(rep(1,nclasses),2.5))
	}else{
		legend("bottomleft",legend=levels(tclass),col=lineColors,bty="n",lty=rep(1,nclasses))
	}
	
	#calculate logodds p-val
	pvalue<-1-pchisq(survdiff(Surv(stime, event)~tclass)$chisq,nclasses-1)
	val<-signif(pvalue,3)
	if(val==0){
		val<-"1e-22"
	}
	legend("bottomright",legend=paste("Log Rank p=",val,sep=""),cex=1.2,,bty="n")
	
	return(val)
}

#args[1] = gene list; each row is the name of one gene located in args[2]
#args[2] = tumor expression matrix; tab-delimited text file where first row must be TCGA barcode ID and first column must be gene names corresponding to args[1]
#args[3] = survival matrix; tab-delimited text file where first column is TCGA barcode ID, second is histology (unused), third is disease specific death (0=alive or unrelated death, 1=disease specific death), and fourth is overall days of survival
#args[4] = directory to put all pdfs in (must already be made)
args = commandArgs(TRUE)

genes = scan(args[1], character()) #read in genes to make KM on
expr = read.table(args[2], sep="\t", row.names=1, header=T) #read in expression matrix
surv.df = read.table(args[3], header=T, sep="\t", row.names=1) #read in survival data matrix

one.third = floor(ncol(expr)/3) #calculate 1/3 of data and floor it
pvals = c()

#for each gene given:
for(i in genes){
	#set top and bottom third vectors w/ minimal values to be replaced
	top.third.vals = rep(-Inf, one.third)
	top.third.names = rep("-", one.third)
	bottom.third.vals = rep(Inf, one.third)
	bottom.third.names = rep("-", one.third)
	
	gene.expr.per.tumor = expr[i,] #get expression data for specific gene
	
	#for each tumor expression value for the specific gene
	for(j in 1:length(gene.expr.per.tumor)){
		#compare to the top third values and see if current value is greater than any of the current top third (meaning it belongs in the top 1/3 over another value already there)
		if(as.numeric(as.character(unlist(gene.expr.per.tumor[j]))) > min(top.third.vals)){ #if so, replace the old, lower value w/ the current, higher value
			min.index = match(min(top.third.vals), top.third.vals)
			
			top.third.vals[min.index] = as.numeric(as.character(unlist(gene.expr.per.tumor[j])))
			top.third.names[min.index] = substr(gsub(".", "-", names(gene.expr.per.tumor[j]), fixed=TRUE), 1, 12)
		}
		#compare to the bottom third values and see if current value is less than any of the current bottom third (meaning it belongs in the bottom 1/3 over another value already there)
		if(as.numeric(as.character(unlist(gene.expr.per.tumor[j]))) < max(bottom.third.vals)){ #if so, replace the old, higher value w/ the current, lower value
			max.index = match(max(bottom.third.vals), bottom.third.vals)
			
			bottom.third.vals[max.index] = as.numeric(as.character(unlist(gene.expr.per.tumor[j])))
			bottom.third.names[max.index] = substr(gsub(".", "-", names(gene.expr.per.tumor[j]), fixed=TRUE), 1, 12)
		}
	}
	
	#set up vector w/ lower or upper 33% labels
	low.or.high.expr = c(rep("Bottom 33%", one.third), rep("Top 33%", one.third))
	#set up vector containing DSS values (1 = disease specific death; 0 = survival or non-DS death)
	DSS.status = surv.df[c(bottom.third.names, top.third.names), 2]
	#set up vector containing number of survival days
	surv.time.days = surv.df[c(bottom.third.names, top.third.names), 3]

	#make pdf file
	pdf(paste0(args[4],"/",i,"_kmPlot_topBottomThird.pdf"))
	#print KM plot and save pvalue
	val = kmPlot(low.or.high.expr, DSS.status, surv.time.days, varName = i)
	dev.off()
	
	#capture pvalue and add to list
	pvals = c(pvals, val)
}

#write out table of matched pvalues to their genes
write.table(cbind(genes, pvals), file=paste0(args[4],"/kmPlot_pvals.txt"), sep="\t", quote=F, row.names=F, col.names=F)

