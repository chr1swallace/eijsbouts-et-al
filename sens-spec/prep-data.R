#!/usr/bin/env Rscript
library(randomFunctions)
args <- getArgs(defaults=list(patt="N1|N2"))
d2 <- "~/share/Data/chic"
list.files(d2)
library(data.table)
od <- "~/bsu/peaky-sens-spec"

message(date())
message("patt: ",args$patt)
## library(peaky)
devtools::load_all("~/RP/pky")
## https://htmlpreview.github.io/?https://github.com/cqgd/pky/blob/master/vignettes/introduction.html

d <- file.path(od, "../peaky/annot")

## output file
bfile <- file.path(od,paste0("BTS-",make.names(args$patt),".RData"))
if(file.exists(bfile))
    stop("output already exists: ",bfile)


## fragments
(load(file.path(d,"hind-position.RData")))
fragments <- hind
setnames(fragments,c("chrom","chromStart","chromEnd","ID"))
fragments$ID <- as.integer(fragments$ID)
fragments$chrom <- as.integer(sub("chr","",fragments$chrom))
fragments <- fragments[!is.na(chrom),]


ifile <- file.path(od,paste0("fragint-",make.names(args$patt),".RData"))
if(!file.exists(ifile)) {
## system(paste0("zcat ",d2,"/Mac_for_Chris.txt.gz | awk -F',' '$3<=10 && $4<=10' | gzip -c > ",d2,"/Mac_for_Chris_1_10.txt.gz"))
## x <- fread(paste0("zcat ",d2,"/Mac_for_Chris.txt.gz | grep -v '[XYMT]' | awk '$3<=10 && $4<=10}'"),sep=",")
    ## x <- fread(paste0("zcat ",d2,"/Mac_for_Chris.txt.gz"))
    x <- fread(paste0("zcat ",d2,"/Mac_for_Chris_1_10.txt.gz"))
    nx <- fread(paste0("zcat ",d2,"/Mac_for_Chris.txt.gz | head -n 10"))
    setnames(x,names(nx))
##     print(colSums(x,na.rm=TRUE))
##       baitID   otherEndID      baitChr  otherEndChr         N0.1         N0.2 
## 1.700843e+14 1.747876e+14 4.102584e+09 4.034804e+09 2.288546e+08 2.314978e+08 
##         N0.3         N1.1         N1.2         N1.3         N2.1         N2.2 
## 1.731265e+08 1.802169e+08 1.688516e+08 1.191081e+08 2.066847e+08 1.837591e+08 
##         N2.3 
##     1.016280e+08
head(x)
setnames(x,"otherEndID","preyID")
cols <- grep(args$patt,names(x),value=TRUE)
    x[,N:=0]
    for(v in cols)
        x$N <- x$N + ifelse(is.na(x[[v]]),0,as.numeric(x[[v]]))
## x$N <- apply(x[,cols,with=FALSE],1,sum)
interactions <- x[N>0,.(baitID,preyID,N)]
summary(interactions)
rm(x); gc()
save(interactions,file=ifile)
}

## load("fragint.RData")
mfile <- file.path(od,paste0("models-",make.names(args$patt),".RData"))
if(file.exists(mfile)) {
    message("loading ",mfile)
    (load(mfile))
} else {
    message("loading ",ifile)
    load(ifile)
    BI = bin_interactions(interactions, fragments, bins=5, max_dist=5e+6)
    rm(interactions); gc()
    BI$interactions <- BI$interactions[!is.na(dist),]
    BI$bins <- BI$bins[1:10,]
    BI$interactions <- copy_bait_covar(BI$interactions,"b.trans_res")
    BI$interactions[is.na(p.trans_res),p.trans_res:=0]
    BI$interactions[,b2b:=p.trans_res!=0]
    models = by(BI$interactions, BI$interactions$dist.bin, model_bin, subsample_size=10000,
                formula_add="b2b + p.trans_res")
    ## models = by(BI$interactions, BI$interactions$dist.bin, model_bin, subsample_size=10000)
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
