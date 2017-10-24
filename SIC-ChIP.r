#!/usr/bin/Rscript

# input arguments
args=commandArgs()
scriptPath=dirname(sub('--file=',"",grep("^--file=",args,value=TRUE)))
if('--help' %in% args){
  writeLines('
SIC-ChIP Shape Index Clustering for ChIP-seq peaks

Usage: /path/to/SIC-ChIP.r --bw=coverage_test.bw --bed=peaks_test.bed --out=results [options]

Alternative: Rscript /path/to/SIC-ChIP.r --bw=coverage_test.bw --bed=peaks_test.bed --out=results [options]

Mandatory options:
  --bw        BigWig file with the coverage function
  --bed       BED file with the coordinates of the peaks
  --out       name of the results folder/files

Options:
  --help      print short help message and exit
  --N         the distance between spline knots for computing the number of local maxima (default --N=20) 
  --toll      the minimum distance allowed between local maxima (default --toll=50)

IMPORTANT: the parameters N and toll should be carefully chosen depending on the peaks dataset
',stdout())
  q(save='no')
}
if((length(grep("^--bw=",args))==0)||(length(grep("^--bed=",args))==0)||(length(grep("^--out=",args))==0)){
  writeLines('
Input files required.

To get help type: /path/to/SIC-ChIP.r --help

or alternatively: Rscript /path/to/SIC-ChIP.r --help
',stdout())
  q(save='no')
}
args=setdiff(args[-(1:grep("^--args",args))],'--help')
args_names=unlist(lapply(strsplit(args,'='),function(arg) arg[1]))
not_valid=!(args_names %in% c('--bw','--bed','--out','--N','--toll'))
if(TRUE %in% not_valid)
  writeLines('
')
for(arg in args_names[not_valid]){
  writeLines(paste0(arg,' is not a valid argument
'),stdout())
}
args=strsplit(sub('--','',args)[!not_valid],'=')
names(args)=lapply(args,function(arg) arg[1])
args=lapply(args,function(arg) arg[2])
bw=args$bw
bed=args$bed
out=args$out
if(!is.null(args$N))
  N=as.numeric(args$N)
if(!is.null(args$toll))
  toll=as.numeric(args$toll)

writeLines('
Loading packages and functions...',stdout())

# load packages
suppressMessages(require(methods))
suppressMessages(require(rtracklayer))
suppressMessages(require(Rcpp))

# load functions
load(paste0(scriptPath,"/SIC-ChIP_functions.RData"))

# compile C++ code
Rcpp::sourceCpp(paste0(scriptPath,"/picco.cpp"))

# create results folder
suppressWarnings(dir.create(paste0("SIC-ChIP_",out)))


#######################################################
#### import coverage function for each peak in bed ####
#######################################################
writeLines('
Importing coverage function for each peak...',stdout())
peaks <- importPeaks.bw(bed,bw)
save(peaks,file=paste0("SIC-ChIP_",out,"/",out,"_peaks.RData"))


########################
#### create indices ####
########################
writeLines('
Computing shape indices...',stdout())
peaks.df.index=data.frame(chr=factor(as.vector(chrom(peaks)),levels=c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX')),
                          start=start(peaks),
                          end=end(peaks),
                          width=width(peaks))
if(!is.null(score(peaks)))
  peaks.df.index$score=score(peaks)
# maximum height
peaks.df.index$height=unlist(lapply(peaks$count,max))
# area
peaks.df.index$area=unlist(lapply(peaks$count,sum))
# full width at half maximum
peaks.df.index$width.height.2=mapply(function(count,height){max(which(count>=height/2))-min(which(count>=height/2))},peaks$count,peaks.df.index$height)
# number of local peaks
if(!exists('N'))
  N=min(20,peaks.df.index$width) # distance knots
if(N>min(peaks.df.index$width)){
  warning('\'N\' is too big. Setting it to the minimum peak width.')
  N=min(peaks.df.index$width)
}
if(N<4){
  if(min(peaks.df.index$width)<4)
    stop('There are very short peaks (<4 nucleotides).')
  warning('\'N\' is too low. Setting to the default value.')
  N=min(20,peaks.df.index$width)
}
if(1 %in% unlist(lapply(peaks$count,function(count) length(unique(count)))))
  warning('There are peaks with constant coverage.')
if(!exists('toll'))
  toll <- min(50,peaks.df.index$width) # minimum distance between maxima
if(toll>min(peaks.df.index$width)){
  warning('\'toll\' is too big. Setting it to the minimum peak width.')
  toll=min(peaks.df.index$width)
}
toll2=0.2 # minimum high (in percentage) of maxima
peaks.df.index$local.peaks=get_local_peaks(peaks$count,N,toll,toll2)
# shape index M divided by the maximum height
peaks.df.index$M.height=unlist(lapply(peaks$count,function(count) get_M(length(count)+1,c(0,count))))/peaks.df.index$height
save(peaks.df.index,file=paste0("SIC-ChIP_",out,"/",out,"_peaks_index.RData"))
# scatter plot
pdf(paste0("SIC-ChIP_",out,"/",out,"_scatterplot.pdf"),10,10)
pairs(peaks.df.index[,c("height","area","width.height.2","local.peaks","M.height")],
      col="black",main='Shape indices',pch=3,
      labels=expression(h,A,w [h/2],p[local],frac(M,h)),cex.labels=1.8)
invisible(dev.off())


####################
#### clustering ####
####################
writeLines('
Clustering peaks...',stdout())
# standardization
index.std=peaks.df.index[,c("height","area","width.height.2","local.peaks","M.height")]
index.std=apply(index.std,2,function(index){(index-mean(index))/sd(index)})
# k-means
K=lapply(1:10,function(k){kmeans(index.std,k,iter.max=100,nstart=10)})
within=unlist(lapply(K,function(K){K$tot.withinss}))
save(K,within,file=paste0("SIC-ChIP_",out,"/",out,"_kmeans.RData"))
# total within sum of squares plot
pdf(paste0("SIC-ChIP_",out,"/",out,"_tot_within_ss.pdf"),10,10)
plot(1:length(within),within/within[1],xlim=c(1,10),ylim=c(0,1),xlab='Number of clusters k',ylab='Total within SS (%)',main='K-means on shape indices',las=1,type='b')
axis(side=1,at=1:10)
invisible(dev.off())

writeLines('
DONE
',stdout())

q(save='no')
