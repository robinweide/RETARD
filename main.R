dataDir <- "~/MEGA/Work/PhD/projects/rDNA/data/rDNA_bedfiles/"
# obtain all files in a list-item  per bedFile: sample/origin/bed
gevondenBedFiles <- list.files(dataDir,pattern = "*.bed",recursive = T)
gevondenBedFiles.list <- strsplit(gevondenBedFiles ,split = "/")
df <- t(as.data.frame(gevondenBedFiles.list)) 
df <-  cbind(df, as.data.frame(unlist(paste0(dataDir,gevondenBedFiles))))
colnames(df) <- c('Donor', 'Sample', 'Ref','Path'); rownames(df) <- 1:nrow(df)

df$Ref <-  gsub(df$Ref, pattern = "_mapped.sorted.bed",replacement = '')
df$Ref <-  gsub(df$Ref, pattern = ".sorted.bed",replacement = '')
df$Ref <-  gsub(df$Ref, pattern = "_0_13314",replacement = '')
df$Ref <-  gsub(df$Ref, pattern = ".*_",replacement = '')
metadata <- read.csv("~/MEGA/Work/PhD/projects/rDNA/data/metadata.csv", sep=";")

sampleInfoDF <- merge(metadata, df, by = "Donor")

###
## Functions
###

getStats <- function(path){
  vec <- NULL
  df <- data.table::fread(as.character(path), data.table = F,verbose = F)
  N <- 1:nrow(df)
  #cat(max(N)," \n")
  if(max(N) > 1e5){
    N <- sample(N, size = 1e4,replace = F)
  }
  for(i in N){
    width <- df[i,3] - df[i,2]
    if(width >= 0){
      vec <- c(vec, rep(df[i,4], width))
    }
  }
  return(list(MED = median(vec), MAD = mad(vec), AVG = mean(vec), STD = sd(vec)))
  rm(vec)
  gc()
}

runRETARD <- function(sampleInfoDF){
  DF <- sampleInfoDF
  DF$MED <- 0
  DF$MAD <- 0
  DF$AVG <- 0
  DF$STD <- 0
  N <- nrow(DF)
  for(i in 1:N){
    cat(paste0("Sample ", i , " of ", N, "\n"))
    l <- getStats(DF[i,]$Path)
    DF[i,]$MED <- l$MED
    DF[i,]$MAD <- l$MAD
    DF[i,]$AVG <- l$AVG
    DF[i,]$STD <- l$STD
  }
  return(DF)
}

### 
## Run!
###

res <- runRETARD(sampleInfoDF)

samples <- unique(res$Sample)
dfr <- NULL
for(Si in 1:length(samples)){
  s <- samples[Si]
  df <- res[res$Sample == s,]
  don <- unique(df$Donor)
  age <- unique(df$Age)
  R <- df[df$Ref == 'RDNAalpha',]
  U <- df[df$Ref == 'U13369',]
  B <- df[df$Ref == 'blacklist',]
  A <- df[df$Ref == 'WAalpha',]
  G <- df[df$Ref == 'WGalpha',]
  dfr <- rbind(dfr,data.frame(Donor = don, Age = age, S = s, RB = R$MED/B$MED, UB = U$MED/B$MED, RA = R$MED/A$MED, UA = U$MED/A$MED, RG = R$MED/G$MED, UG = U$MED/G$MED))
}
dfr <- as.data.frame(dfr)
mdfr <- reshape2::melt(dfr, id.vars = c('Donor', 'S', 'Age'))
ggplot(mdfr , aes(y = value, x = variable, col = Donor)) + geom_point() + geom_step() + facet_grid(. ~ Age) + labs(x = "Method", y = "rDNA dosage", facet = "Age")
dfr[!grepl(dfr$S, pattern = "blood", ignore.case = T),]
plot(dfr[grepl(dfr$S, pattern = "blood", ignore.case = T),]$Age,dfr[grepl(dfr$S, pattern = "blood", ignore.case = T),]$RB, type = 'l')

### Per-age plot
dat <- dfr %>% group_by(Age) %>% summarise(Donor = unique(Donor),m = mean(UB), s = sd(UB))

newDat <- NULL
for(i in 1:nrow(dat)){
  newDat <- rbind(newDat, cbind(dat[i,1:4], as.data.frame(rnorm(100, dat[i,]$m, dat[i,]$s))))
}
colnames(newDat)[5] <- "raw"
plot(loess.smooth(span = 1, newDat$Age, newDat$raw), type = 'l', ylim = c(70,175), main = ""  , xlab = "Age", ylab = "rDNA dosage")
points(dfr$Age, dfr$UB, pch = 20,col= ifelse(grepl(dfr$S, ignore.case=T, pattern = "blood"), 'red','black'))
legend(x = 3, y = 177, legend = c('intestine', 'blood', "fit"), col = c(1,2,1), pch = c(20,20,NA), lwd = c(NA,NA,1), bty = "n", seg.len=.75)

