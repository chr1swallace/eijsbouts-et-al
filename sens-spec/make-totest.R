#!/usr/bin/env Rscript
library(randomFunctions)
library(magrittr)
args <- getArgs(defaults=list(patt="N2"))
d2 <- "~/share/Data/chic"
list.files(d2)
library(data.table)
od <- "~/bsu/peaky-sens-spec"

library(peaky)
## https://htmlpreview.github.io/?https://github.com/cqgd/pky/blob/master/vignettes/introduction.html

files <- list.files("~/bsu/peaky-sens-spec",pattern="BTS-N",full=TRUE)
message("files found: ")
print(files)

totest <- lapply(files, function(f) {
    message(f)
    (load(f))
    ## then save BTS, or only that part on chr 1?
    unique(BTS[abs(dist)<5e+6 & fdr.res<0.1,.(b.chr,baitID)])
})

length(totest)
lapply(totest,nrow)

m <- unique(rbindlist(totest))
nrow(m)

table(m$b.chr)

length(unique(unlist(totest)))

totest <- m[b.chr==1,]$bait
cat(totest,file="~/bsu/peaky-sens-spec/totest.txt",sep="\n")


## check where we are
dirs <- lapply(files, function(f) gsub("BTS-|.RData","",f))
done <- lapply(dirs, list.files)  %>% lapply(., function(f) sub(".csv","",f)  %>% as.numeric())
names(done) <- sapply(dirs,basename)
sapply(done,length)

done1 <- lapply(done, intersect, totest)
length(totest)
sapply(done1,length)
