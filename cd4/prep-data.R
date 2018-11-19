#!/usr/bin/env Rscript
library(randomFunctions)
args <- getArgs(defaults=list(patt="aCD4val")) # aCD4pchic / aCD4val / nCD4pchic / nCD4val
d2 <- "/home/cew54/bsu/peaky/datasets"
library(data.table)

od <- file.path("~/bsu/peaky/pky",args$patt)
bfile <- file.path(od,"BTS.RData")
if(file.exists(bfile))
    stop("output file exists: ",bfile)

message(date())
message("patt: ",args$patt)
## library(peaky)
library(devtools)
load_all("~/RP/pky")
## https://htmlpreview.github.io/?https://github.com/cqgd/pky/blob/master/vignettes/introduction.html

d <- "/mrc-bsu/scratch/cew54/peaky/annot"
## fragments
(load(file.path(d,"hind-position.RData")))
fragments <- hind
setnames(fragments,c("chrom","chromStart","chromEnd","ID"))
fragments$ID <- as.integer(fragments$ID)
fragments$chrom <- as.integer(sub("chr","",fragments$chrom))
fragments <- fragments[!is.na(chrom),]

if(!file.exists(od))
    dir.create(od)
ifile <- file.path(od,paste0("fragint.RData"))
if(file.exists(ifile)) {
    message("loading ",ifile)
    load(ifile)
} else {
    x <- fread(paste0(d2,"/",args$patt,"_mini.csv"))
    ## add distance
    head(x)
    maxscore <- x[,.(mx=max(score)),by="baitID"]
    usebaits <- maxscore[mx>=5,]$baitID
    interactions <- x[baitID %in% usebaits &
                      baitID %in% fragments$ID &
                      preyID %in% fragments$ID,.(baitID,preyID,N)]
summary(interactions)
rm(x); gc()
## save(interactions,file=ifile)
}

## load("fragint.RData")
mfile <- file.path(od,paste0("models.RData"))
if(file.exists(mfile)) {
    message("loading ",mfile)
    (load(mfile))
} else {
    BI = bin_interactions(interactions, fragments, bins=5, max_dist=5e+6)
    BI$interactions <- BI$interactions[!is.na(dist),]
    BI$bins <- BI$bins[1:10,]
    BI$interactions <- copy_bait_covar(BI$interactions,"b.trans_res")
    BI$interactions[is.na(p.trans_res),p.trans_res:=0]
    BI$interactions[,b2b:=p.trans_res!=0]
    models = by(BI$interactions, BI$interactions$dist.bin, model_bin, subsample_size=10000,
                formula_add="b2b + p.trans_res")
    ## save(BI,models, file=mfile)
}

## if(file.exists(bfile)) {
##     message("loading ",bfile)
##     (load(bfile))
## } else {
residuals = lapply(models, "[[", "residuals")
bins = split(BI$interactions, BI$interactions$dist.bin)
BTS = split_baits(bins, residuals)
save(BTS,file=bfile)
BTS <- BTS[b.chr==1,] ## about 10% of size
save(BTS,file=sub("BTS","BTSchr1",bfile))

## }

## then save BTS, or only that part on chr 1?
## BTS1 <- BTS[b.chr==1,] ## about 10% of size
## totest <- BTS1[


## TODO



## relevant_bait = BTS[baitID==418]
## zoom = relevant_bait[abs(relevant_bait$dist)<1e6]
## plot(x=zoom$dist, 
##      y=zoom$residual, 
##      xlab="Distance from bait (bp)", 
##      ylab="Adjusted readcount",
##      main=paste("Bait",unique(zoom$baitID)))
## abline(v=0,col="red",lty=2)



## relevant_bait = BTS[baitID==418] 
## omega_power = -5
## PKS = peaky(zoom, omega_power, iterations=1e6)
## P = interpret_peaky(zoom, PKS, omega_power)[,.(baitID,preyID,mppc=rjmcmc_pos)]

## par(mfrow=c(3,1))
## plot(x=P$dist, xlab="Distance from bait (bp)",
##      y=P$residual, ylab="Adjusted readcount")

## plot(x=P$dist, xlab="Distance from bait (bp)",
##      y=P$beta_mean, ylab="Mean contact strength",
##      col="green")

## plot(x=P$dist, xlab="Distance from bait (bp)",
##      y=P$rjmcmc_pos, ylab="MPPC",
##      col="blue")
