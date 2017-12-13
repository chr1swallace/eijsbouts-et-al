library(data.table)
library(ggplot2)
library(magrittr)
library(cowplot)

f <- function(x) {
    k <- x[,.(maxd=max(abs(dist)),nfrag=.N),by="baitID"]
    lapply(k[,.(maxd,nfrag)],summary)
}

##d <- "/scratch/wallace/peaky/derived"
d <- "/mrc-bsu/scratch/cew54/peaky/summary/derived"
list.files(d)
peaky.files <- list.files(d,pattern="peaky.RData")
chic.files <- list.files(d,pattern="chicago")
(load(file.path(d,"cao.RData"))) ## cao et al enh-promoter interactions from total CD4
source("~/Projects/peaky/common.R")
##x <- readRDS("/mrc-bsu/scratch/cew54/peaky/summary/tables/DT_aCD4val_20e6_corr.rds")
names(title)
nm <- "actval"
for(nm in names(title)) {
   message(" !!! ",nm) 
   fnm <- file.path(d,paste0(nm,"-joined.RData"))
   ## if(file.exists(fnm)) {
   ## (load(fnm))
   ## } else {
   ## load raw files
   (load(file.path(d,paste0(nm,"-raw.RData"))));
   (load(file.path(d,paste0(nm,"-peaky.RData"))));
   (load(file.path(d,paste0(nm,"-chicago.RData"))));
   (load(file.path(d,"hind-position.RData"))); hind[,hindID:=as.integer(hindID)]
   ## merge - every thing paired should be keyed so that by not needed
   setkey(chic,baitID,preyID)
   x <- merge(raw,chic,all.x=TRUE,by=c("baitID","preyID"))
   x <- merge(x,peaky,by=c("baitID","preyID"))
   x <- addpositions(x)
   
   if(grepl("val",nm)) {
       (load(file.path(d,"b2g.RData")))
       b2g <- b2g[biotype %in% c("protein_coding"),]
       (load(file.path(d,"rnaseq.RData")))
       if(grepl("act",nm)) {
           rnaseq[,trans:=act_1 + act_2]
       } else {
           rnaseq[,trans:=non_1 + non_2]
       }
       b2g <- merge(b2g[,.(id,gene,baitID)],rnaseq[,.(id,trans)],by="id")
       setnames(b2g,"baitID","preyID")
       b2g <- b2g[,.(z=mean(log2(trans+1))),by="preyID"]
       x <- merge(x,b2g,by="preyID",all.x=TRUE)
       x[,y:=as.integer(!is.na(z))] # z only not NA if corresponds to a bait in the promcap expt
   }
   
   if(grepl("prom",nm)) {
       o <- (load(file.path(d,paste0("chromhmm-",substr(nm,1,3),".RData"))))
       chromhmm <- eval(as.symbol(o))
       rm(list=o)
       chromhmm[,enhprom:=pmin(1,E4+E5+E6+E7+E8+E9+E10+E11)]
       x <- merge(x,chromhmm[,.(hindID,enhprom,E14)],
                  by.x="preyID",by.y="hindID",all.x=TRUE)
       x[,z:=enhprom]
       x <- merge(x,cao,by=c("baitID","preyID"),all.x=TRUE)
       x[abs(dist)<1e+6,y:=ifelse(is.na(cao),0,1)]
    }
   
   x[is.na(chicago),chicago:=0]
   x[,maxchic:=max(chicago),by="baitID"]
   x[,maxres:=max(residual),by="baitID"]
   ux <- unique(x,by="baitID")
   message("Number of baits passing 5/2 thresholds for max chicago score/NB residual")
   print(with(ux,table(r=maxres>=2,c=maxchic>=5)))
   save(x,file=fnm)


   ## fnm <- file.path(d,paste0(nm,"-joined.RData"))
   ## (load(fnm))
   ## print(dim(x))
   fnm <- file.path(d,paste0(nm,"-joined-fm.RData"))
   x <- x[maxchic>=5,]
   save(x,file=fnm)
 }

message("all is well")
print(date())
