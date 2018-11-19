#!/usr/bin/env Rscript
library(randomFunctions)
id <- "~/bsu/peaky-sens-spec"
args <- getArgs(defaults=list(patt="N0.N1",taskid=1,file=file.path(id,"totest.txt")),numeric=c("taskid"))
d2 <- "~/share/Data/chic"

library(data.table)
totest <- scan(args$file,what="")
bait <- as.integer(totest[ args$taskid ])

od <- file.path(id,args$patt)
of <- file.path(od,paste0(bait,".csv"))
if(file.exists(of))
    stop("file exists: ",of)
if(!file.exists(od))
    dir.create(od)

library(peaky)
## https://htmlpreview.github.io/?https://github.com/cqgd/pky/blob/master/vignettes/introduction.html

d <- "/mrc-bsu/scratch/cew54/peaky/summary/derived"
## ## fragments
## (load(file.path(d,"hind-position.RData")))
## fragments <- hind
## setnames(fragments,c("chrom","chromStart","chromEnd","ID"))
## fragments$ID <- as.integer(fragments$ID)
## fragments$chrom <- as.integer(sub("chr","",fragments$chrom))
## fragments <- fragments[!is.na(chrom),]

bfile <- paste0(id,"/BTS-",args$patt,".RData")
message("loading ",bfile)
(load(bfile))

message("looking for bait ",bait)
relevant_bait = BTS[baitID==bait]
with(BTS,table(baitID==bait))

zoom = relevant_bait[abs(relevant_bait$dist)<5e6, ] # 10mb window - 5mb either side
zoom

plot(x=zoom$dist, 
     y=zoom$residual, 
     xlab="Distance from bait (bp)", 
     ylab="Adjusted readcount",
     main=paste("Bait",unique(zoom$baitID)))
abline(v=0,col="red",lty=2)

omega_power = -5
PKS = peaky(zoom, omega_power, iterations=1e7) # 10 million iterations. won't save
P = interpret_peaky(zoom, PKS, omega_power)[,.(baitID,preyID,mppc=rjmcmc_pos)]

fwrite(P,file=of)
