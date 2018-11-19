#!/usr/bin/env Rscript

library(randomFunctions)
library(magrittr) 
args <- getArgs(defaults=list(patt="aCD4val")) # aCD4pchic / aCD4val / nCD4pchic / nCD4val
library(data.table)
od <- file.path("~/bsu/peaky-orig/chrisw",args$patt)

files.1 <- list.files(file.path(od,"rep-1"))
files.2 <- list.files(file.path(od,"rep-2"))
files.3 <- list.files(file.path(od,"rep-3"))
files.4 <- list.files(file.path(od,"rep-4"))

## first 20 mill
if(!length(files.1) || !length(files.2)) 
    stop("no files found with 2 reps")
files <- intersect(files.1,files.2)
message("files found in both reps: ",length(files))

corr <- fread(file.path(od,"mppc-cor.txt"))
files <- intersect(files,paste0(corr[cr>=0.75,]$baitID,".csv"))

getfiles <- function(files) {
    x <- lapply(files, function(f) {
        z <- fread(paste0("cat ",od,"/rep-*/",f))
        z <- z[mppc!="mppc",]
        z[,.(mppc=mean(as.numeric(mppc)),
             gamma_mean=mean(as.numeric(gamma_mean))),by=c("baitID","preyID","residual")]
    })  %>% rbindlist()
}
x <- getfiles(files)

fwrite(x,file.path(od,"mppc.csv"))

