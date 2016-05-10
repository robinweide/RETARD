readcovdata <- function(files, col=6, names=NULL) {
	require(BiocParallel)
	bp.param = MulticoreParam(workers=12)
	
	d <- bplapply(files, function(gzf) {
	  s <- read.table(gzfile(gzf))
	  md <- regmatches(gzf, regexec("(DRR0168[35]x)(\\d+)M-(\\d)-(rRNA|gcbins-.)", gzf))[[1]]
	  ag <- aggregate(s[,col], by=list(regionid=paste(s[,1], s[,2],s[,3], sep=":"), binwidth=s[,3]-s[,2]), mean)
	  ag$sample <- md[2]
	  ag$Mreads <- md[3]
	  ag$rep <- md[4]
	  ag$type <- md[5]
	  ag
    }, BPPARAM=bp.param)
	if(!is.null(names)) names(d) <- names else names(d) <- files
	melt(d,measure.vars="x", value.name="meancov")
}


library(gtools)
library(reshape2)
library(dplyr)
library(ggplot2)


#we have two samples: 
samples <- c("DRR01683x","DRR01685x")
#we have 5 tiers with 1-10 random sammples:
n <- c(350,rep(100,2), rep(50,4), rep(20,8), rep(10,10), rep(5,10))
m <- c(0,0:1, 0:3, 0:7,0:9, 0:9)
s <- paste0(n,"M-",m)

#the reference has 3 random sets
ref <- c("gcbins-1","gcbins-2","gcbins-3")
gcdata <- lapply(ref, function(n) {
  gcd <- read.delim(paste0(n,"-gc.txt"))
  gc <- gcd$X5_pct_gc
  names(gc) <- paste(gcd[,1], gcd[,2], gcd[,3], sep=":")
  gc
})
names(gcdata) <- ref

#ribosomal gc data
rgc <-  { gcd <- read.delim("rrna-w500-gc.bed");
  gc <- gcd$X6_pct_gc
  names(gc) <- paste(gcd[,1], gcd[,2], gcd[,3], sep=":")
  gc
}


#read the ribosomal coverage data
rcovfiles <- list.files(pattern=".*rRNA-rrna-w500.txt.gz$")
rcov <- readcovdata(rcovfiles, col=6, sub("-rrna-w500.txt.gz$","", rcovfiles))
rcov$gc <- rgc[as.character(rcov$regionid)]


#read the reference data
refcov <- NULL
for (re in ref) {
  reffiles <- list.files(pattern=paste0(".*-", re, ".txt.gz$"))
  rfcov <- readcovdata(reffiles, col=5, sub("-gcbins-..txt.gz$","", reffiles))
  #bind gc data
  rfcov$gc <- gcdata[[1]][as.character(rfcov$regionid)]
  rfcov$refname <- re
  refcov <- rbind(refcov, rfcov)
}

save(refcov, rcov, file="covdata.Rda");

#fit a lowess per sample on the refdata
fits <- with(refcov, by(refcov, list(L1, refname), function(x) loess(meancov ~ gc, data = x)))

for (r in ref) rcov[,r] <- sapply(1:nrow(rcov), function(i) {
  predict(fits[[sub("-rRNA$","", rcov$L1[i]), r]], data.frame(gc=rcov$gc[i]))
})
rcov$ref <- rowMeans(rcov[,c("gcbins-1","gcbins-2","gcbins-3")])
rcov$ratio <- rcov$meancov / rcov$ref

ggplot(subset(rcov, rcov$gc < 0.82), aes(x=rep, y=ratio, fill=sample)) + geom_boxplot() + facet_wrap(~Mreads) + ggsave("bw-perregions-3ref-rRDNA.png", width=240,height=120, unit="mm", dpi=96)

summarise(group_by(subset(rcov, rcov$gc < 0.82), sample, Mreads, rep), wmean=weighted.mean(ratio, binwidth)) -> aap
sumres <- as.data.frame(aap)
sumres$Mreads <- factor(paste(sumres$Mreads, "M reads"), levels=paste(unique(sort(as.numeric(sumres$Mreads))),"M reads"), ordered=T)
ggplot(sumres, aes(x=sample, y=wmean, color=sample)) + geom_point(size=2) + facet_grid(. ~Mreads) + theme(axis.text.x = element_blank()) + ggsave("wmean-3ref-rRDNA.png", width=180,height=90, unit="mm")

