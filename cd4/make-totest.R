#!/usr/bin/env Rscript
library(randomFunctions)
args <- getArgs(defaults=list(patt="aCD4val"))
d2 <- "~/share/Data/chic"
list.files(d2)
library(data.table)
od <- file.path("~/bsu/peaky/chrisw",args$patt)

library(peaky)
## https://htmlpreview.github.io/?https://github.com/cqgd/pky/blob/master/vignettes/introduction.html

files <- file.path(od,"BTS.RData") #list.files(od,pattern="BTS-",full=TRUE)
## message("files found: ")
## print(files)

(load(files))

## add chicago score
d2 <- "/home/cew54/bsu/peaky/raw_copy_subliminal/datasets"
x <- fread(paste0(d2,"/",args$patt,"_mini.csv"))
m <- merge(BTS,x[,.(baitID,preyID,score)],by=c("baitID","preyID"))
m <- m[abs(dist) < 5e+6, ]
m <- m[,.(score=max(score)), by="baitID"]

## then save BTS, or only that part on chr 1?
totest <- unique(m$baitID)

cat(totest,file=file.path(od,"totest.txt"))
