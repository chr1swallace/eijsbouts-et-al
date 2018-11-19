#!/usr/bin/env Rscript

library(randomFunctions)
library(magrittr) 
args <- getArgs(defaults=list(patt="aCD4val")) # aCD4pchic / aCD4val / nCD4pchic / nCD4val
library(data.table)
setwd("~/Projects/eijsbouts-et-al/cd4")
od <- file.path("output",args$patt)

files.1 <- list.files(file.path(od,"rep-1"))
files.2 <- list.files(file.path(od,"rep-2"))
files.3 <- list.files(file.path(od,"rep-3"))
files.4 <- list.files(file.path(od,"rep-4"))

## first 20 mill
if(!length(files.1) || !length(files.2)) 
    stop("no files found with 2 reps")
files <- intersect(files.1,files.2)
message("files found in both reps: ",length(files))
getfiles <- function(files,rep) {
    x <- lapply(files, function(f) fread(file.path(od,paste0("rep-",rep),f))[,.(baitID,preyID,mppc)])
    names(x) <- sub(".csv","",files)
    x
}
d1 <- getfiles(files,1)
d2 <- getfiles(files,2)

merger <- function(d1,d2,suf1=1,suf2=2) {
    mapply(merge,d1,d2,
           MoreArgs=list(by=c("baitID","preyID"),suffixes=paste0(".",c(suf1,suf2))),
           SIMPLIFY=FALSE)  %>% rbindlist()
}
m <- merger(d1,d2,1,2)
mc <- m[,.(cr=cor(mppc.1,mppc.2)),by="baitID"]
if(interactive()) {
    hist(mc$cr)
    summary(mc$cr)
}
totest2 <- mc[cr<0.75,]$baitID

dim(mc)
sum(mc$cr>0.75)
mean(mc$cr>0.75)

## repeats if they exist
files <- intersect(files.3,files.4)
totest2 <- setdiff(totest2,as.numeric(sub(".csv","",files)))

message("further reps needed: ",length(totest2)," / ",nrow(mc)," baits")

cat(totest2,file=file.path(od,"totest2.txt"))

d3 <- getfiles(files,3)
d4 <- getfiles(files,4)

## add to d1,d2
if(length(d3) && length(d4)) {
    fi <- intersect(names(d3),names(d4))
	if(!length(fi))
stop("no files found with 2 reps")
    m1 <- merger(d3[fi],d4[fi],3,4)
    mm <- merge(m,m1,by=c("baitID","preyID"))
    mm[,mppc.1:=(mppc.1+mppc.2)/2]
    mm[,mppc.2:=(mppc.3+mppc.4)/2]
    m <- rbind(m[!(baitID %in% mm$baitID),],
               mm[,.(baitID,preyID,mppc.1,mppc.2)])
    mc <- m[,.(cr=cor(mppc.1,mppc.2)),by="baitID"]
}

dim(mc)
sum(mc$cr>0.75)
mean(mc$cr>0.75)

message("baits ok to go: ",sum(mc$cr>=0.75)," / ",nrow(mc)," baits")
fwrite(mc,file.path(od,"mppc-cor.txt"))

m[,mppc:=(mppc.1 + mppc.2)/2]
baits.ok <- mc[cr > 0.75,]$baitID
m <- m[baitID %in% baits.ok,.(baitID,preyID,mppc)]

fwrite(m,file.path(od,"mppc.csv"))
