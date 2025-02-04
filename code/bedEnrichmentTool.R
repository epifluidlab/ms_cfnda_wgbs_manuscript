# TODO: Add comment
# 
# Author: yaping
###############################################################################

for (e in commandArgs(TRUE)) {
	ta = strsplit(e,"=",fixed=TRUE)
	if(! is.na(ta[[1]][2])) {
		if(ta[[1]][1] == "wd"){
			wd<-ta[[1]][2]  ## directory to make plot file
			
		}
		if(ta[[1]][1] == "input"){
			input<-ta[[1]][2]
		}
		if(ta[[1]][1] == "output"){
			output<-ta[[1]][2]
		}
		if(ta[[1]][1] == "outpic"){
			outpic<-ta[[1]][2]
		}
	}
}

setwd(wd)

#print OUT "$file\t$overlap\t$not_overlap\t$overlap_r_mean\t$overlap_r_sd\t$not_overlap_r_mean\t$not_overlap_r_sd\n";



##input the output file into R, calculate ratio, Odds ratio, p value (chi-square, fisher test), output barplot figure..

inputData<-read.table(input,sep="\t",header=F)
ratio=inputData[,2]/inputData[,4]
odd_ratio=(inputData[,2]*inputData[,6])/(inputData[,3]*inputData[,4])
chisquare<-NULL
fisher<-NULL
for(i in 1:length(inputData[,1])){
	cp<-chisq.test(rbind(c(inputData[i,2],inputData[i,3]),c(inputData[i,4],inputData[i,6])))$p.value
	if((inputData[i,2]+inputData[i,3])<30000){
		fp<-fisher.test(rbind(c(inputData[i,2],inputData[i,3]),c(inputData[i,4],inputData[i,6])))$p.value
	}else{
		fp<-NA
	}
	
	chisquare<-c(chisquare,cp)
	fisher<-c(fisher,fp)
}
d<-cbind(inputData,ratio,odd_ratio,chisquare,fisher)
colnames(d)<-c("fileName","overlap","not_overlap","random.overlap.mean","random.overlap.sd","random.not_overlap.mean","random.not_overlap.sd","ratio","odd_ratio","pvalue.chisqaure","pvalue.fisher")
write.table(d,output,quote=F,row.names=F,col.names=T,sep="\t")

name.file<-gsub("\\w+$","",inputData[,1])


#error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
#	if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
#		stop("vectors must be same length")
#	arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
#}

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
	if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
		stop("vectors must be same length")
	arrows(x-lower,y, x+upper, y, angle=180, code=3, length=length, ...)
}

pdf(outpic, paper="special", height=1.5*dim(inputData)[1], width=4+0.2*max(nchar(as.character(name.file))))
par(oma=c(0, 0, 0, 0),mar=c(4, 7+0.4*max(nchar(as.character(name.file))), 3, 1))
barx <- barplot(rbind(inputData[,2],inputData[,4]), beside = TRUE,horiz=TRUE,
         col = c("red", "grey"),
         names.arg = name.file,ylab="",xlab="Count",main="",xlim=c(0,max(rbind(inputData[,2],inputData[,4]+5*1.96*inputData[,5]/10),na.rm=T)),las=2, legend.text=c("Observed","Random"))
error.bar(inputData[,4],barx[2,], 1.96*inputData[,5]/10)
dev.off()

