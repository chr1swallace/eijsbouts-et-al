#!/usr/bin/env Rscript
library(randomFunctions)
library(data.table)
args <- getArgs(defaults=list(patt="aCD4pchic",ifile="aCD4pchic-rep4.txt",taskid=1,rep=4,mult=10,bait=99653),
                numeric=c("taskid","mult","bait"))
d2 <- "/home/cew54/bsu/peaky/datasets"
print(args)

## input file set up, identify baits todo
od <- file.path("output",args$patt)
bfile <- file.path(od,"BTS.RData")
if(args$bait>0) {
    baits <- args$bait
} else {
    if(args$ifile!="") {
        totest <- scan(file.path("~/Projects/eijsbouts-et-al/cd4",args$ifile))
    } else {
        if(args$rep<3) {
            totest <- scan(file.path(od,"totest.txt"),what="")
        } else {
            totest <- scan(file.path(od,"totest2.txt"),what="")
        }
    }
    idx <- args$taskid*args$mult - 1:args$mult + 1
    baits <- as.integer(totest[idx])
}

if(length(baits)>1)
    baits <- sample(baits)## randomise, so long running baits don't always stop the rest of the pack
## output files, remove baits from list if output already exists
od <- file.path(od,paste0("rep-",args$rep))
of <- file.path(od,paste0(baits,".csv")) 
drop <- which(file.exists(of))
if(length(drop)) {
    baits <- baits[-drop]
    of <- of[-drop]
}
if(!length(baits))
    stop("all output files exist")
if(!file.exists(od))
    dir.create(od)
message("baits to run: ",length(baits))

library(peaky)
## https://htmlpreview.github.io/?https://github.com/cqgd/pky/blob/master/vignettes/introduction.html

message("loading ",bfile)
(load(bfile))
## totest <- scan("~/scratch/peaky/totest.txt",what="")

message("args: "); print(args)
message("baits: "); print(baits)
for(i in seq_along(baits)) {
    message("bait: ",baits[i])
	if(file.exists(of[i])) ## just in case another job has created the file in meantime
next
    zoom = BTS[baitID==baits[i]]
    message("nfrags: ",nrow(zoom))
## zoom = relevant_bait[abs(relevant_bait$dist)<5e6, ] # 10mb window - 5mb either side
if(interactive()) {
    plot(x=zoom$dist, 
     y=zoom$residual, 
     xlab="Distance from bait (bp)", 
     ylab="Adjusted readcount",
     main=paste("Bait",unique(zoom$baitID)))
abline(v=0,col="red",lty=2)
}

omega_power = -5
PKS = peaky(zoom, omega_power, iterations=1e7) # 10 million iterations. won't save (?) - if correlation between mppc from two runs is high, the combined inference would be the pointwise average
P = interpret_peaky(zoom, PKS, omega_power)[,.(baitID,preyID,residual,beta_mean,mppc=rjmcmc_pos)]

    message("saving to: ",of[i])
    fwrite(P,file=of[i])
}
