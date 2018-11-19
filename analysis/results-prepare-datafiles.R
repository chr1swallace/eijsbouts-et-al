#!/usr/bin/env Rscript

library(data.table)
#ap <- readRDS("/scratch/wallace/peaky/tables/DT_aCD4pchic_double5e6.rds")
#av <- readRDS("/scratch/wallace/peaky/tables/DT_aCD4val_double5e6.rds")
#np <- readRDS("/scratch/wallace/peaky/tables/DT_nCD4pchic_double5e6.rds")
#x <- readRDS("/scratch/wallace/peaky/tables/DT_nCD4val_quadruple5e6.rds")
#head(x)

## new data location
##x <- readRDS("/mrc-bsu/scratch/cew54/peaky/summary/tables/DT_aCD4val_20e6_corr.rds")

dd <- "/mrc-bsu/scratch/cew54/peaky/derived" #"/scratch/wallace/peaky/tables"
files <- c("aCD4pchic-mppc.csv",
           "aCD4val-mppc.csv",
           "nCD4pchic-mppc.csv",
           "nCD4val-mppc.csv")
nm <- sub("-mppc.csv","",files)
files <- file.path(dd,files)
names(files) <- nm

for(i in seq_along(files)) {
    if(!file.exists(files[i]))
        stop("file not found: ",files[i])
}

## CORR <- vector("list",length(files))
## names(CORR) <- names(files)
## d2 <- "/home/cew54/bsu/peaky/raw_copy_subliminal/datasets"
ad <- "~/bsu/peaky/annot"
(load(file.path(ad,"hind-position.RData"))); hind[,hindID:=as.integer(hindID)]
hind[,mid:=(start+end)/2]
addpositions <- function(x) {
    x <- merge(x,hind[,.(hindID,mid)],by.x="baitID",by.y="hindID")
    x <- merge(x,hind[,.(hindID,mid)],by.x="preyID",by.y="hindID",suffixes = c(".bait",".prey"))
    x[,dist :=  mid.prey - mid.bait]
    return(x)
}
## (load(file.path(d3,"cao.RData"))) ## cao et al enh-promoter interactions from total CD4
## source("~/Projects/peaky/common.R")

do <- "~/bsu/peaky/datasets"
for(i in 1:3) {# seq_along(files)) {
    f <- files[i]
    nm <- names(files)[i]
    message("\n\nreading ",nm," - ",f)

    raw <- fread(f) 
    counts <- fread(paste0(do,"/",nm,"_mini.csv"))[,.(baitID,preyID,N,score)]
    data <- merge(raw,counts,by=c("baitID","preyID"))
    
    ## message("subsetting columns")
    ## setnames(x,
    ##      c("rjmcmc_pos_g-4.69897","rjmcmc_g-4.69897",
    ##        "beta_mean_g-4.69897", "predicted_g-4.69897"),
    ##      c("mppc","mppi","beta.post","peaky.pred"))
    ## x[,maxcor:=pmax(corr_5v5,corr_10v10)]
    ## use <- which(x$maxcor>0.75)

    ## message("saving column subsets")
    ## setkey(raw,baitID,preyID)
    ## save(raw,file=file.path(d,"../derived", paste0(nm,"-raw-locorr.RData")))
    setkey(data,baitID,preyID)
    data <- addpositions(data)
    setnames(data,"score","chicago")
    save(data,file=file.path(dd,"../analysis", paste0(nm,".RData")))

    ## peaky <- x[-use,.(baitID,preyID,mppc,mppi,beta.post,peaky.pred)]
    ## setkey(peaky,baitID,preyID)
    ## save(peaky,file=file.path(d,"../derived", paste0(nm,"-peaky-lowcorr.RData")))
    ## peaky <- x[use,.(baitID,preyID,mppc,mppi,beta.post,peaky.pred)]
    ## setkey(peaky,baitID,preyID)
    ## save(peaky,file=file.path(d,"../derived", paste0(nm,"-peaky.RData")))


    ## chicago score
    ## if("score" %in% names(x)) {
    ##     message("Chicago score found")
    ##     setnames(x,"score","chicago") 
    ##     chic <- x[-use,.(baitID,preyID,chicago,B,T)]
    ##     save(chic,file=file.path(d,"../derived", paste0(nm,"-chicago-lowcorr.RData")))
    ##     chic <- x[use,.(baitID,preyID,chicago,B,T)]
    ##     save(chic,file=file.path(d,"../derived", paste0(nm,"-chicago.RData")))
    ## } else {
        
    ##     message("attempting to fill in act-validation for chicago")
    ##     library(Chicago)
    ##     ##val <- readRDS(file.path(d,"../../raw/aCD4val_merge.Rds"))
    ##     val <- readRDS(file.path(d,"../../raw/aCD4val_merge_complete.Rd"))
    ##     y <- val@x[,.(baitID,otherEndID,score)]
    ##     setnames(y,"otherEndID","preyID")
    ##     setnames(y,"score","chicago")
    ##     val <- fread(file.path(CD4CHIC.DATA,"validation_peakmatrix","peakMatrix_validation_cutoff0.txt"))
    ##     setnames(val,"oeID","preyID")
    ##     x1 <- merge(x[,.(baitID,preyID,maxcor)],y,by=c("baitID","preyID"))
    ##     x2 <- merge(x[,.(baitID,preyID,maxcor)],y, by.x=c("baitID","preyID"),
    ##                 by.y=c("preyID","baitID"))
    ##     dim(x1)
    ##     dim(x2)
    ##     x <- unique(rbind(x1,x2))
    ##     dim(x)
    ##     use <- which(x$maxcor>0.75)
    ##     chic <- x[-use,.(baitID,preyID,chicago)]
    ##     save(chic,file=file.path(d,"../derived", paste0(nm,"-chicago-lowcorr.RData")))
    ##     chic <- x[use,.(baitID,preyID,chicago)]
    ##     save(chic,file=file.path(d,"../derived", paste0(nm,"-chicago.RData")))
    ## }

    ## cr <- x[,.(maxcor=max(maxcor),chicago.mx=max(chicago),chicago.n=sum(chicago>5)),by="baitID"]
    ## summary(cr$maxcor)
    ## cr[,experiment:=nm]
    ## CORR[[i]] <- cr

}

## save(CORR,file=file.path(d,"../derived/corr.RData"))
## cr <- rbindlist(CORR)
## library(stargazer)
## tt <- with(cr,table(experiment=experiment, "corr > 0.75"=ifelse(maxcor>0.75,">=0.75","<0.75")))
## stargazer(tt)
