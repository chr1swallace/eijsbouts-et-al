source("~/Projects/peaky/common-v2.R")
source("~/Projects/peaky/viz.R")
                #plotone(xb,zoom=TODO[[g]]$zooms[[i]])
## source("~/Projects/peaky/extract-mppc-credsets.R")
library(magrittr)
d <- "/mrc-bsu/scratch/cew54/peaky"

(load(file.path(d,"annot","b2g.RData")))
head(b2g)
baits <- b2g$baitID

(load(file.path(d,"annot","hind-position.RData")))
hind[,hindID:=as.integer(hindID)]
## (load(file.path(d,"rnaseq.RData")))
## head(rnaseq)
## b2g <- merge(b2g,rnaseq,by="id")

## select some examples
## nm <- "actval"
## for(nm in names(title)) {
##    message(" !!! ",nm) 
##    fnm <- file.path(d,paste0(nm,"-joined-fm.RData"))
##    (load(fnm))
##    y <- x[,.(chicago=sum(chicago>5)),by="baitID"]
##    y <- y[chicago>quantile(chicago,0.99),]
##    y <- merge(y,b2g,by="baitID",all.x=TRUE)

addstartend <- function(x) {
    x <- merge(x,hind[,.(hindID,start,end)],by.x="baitID",by.y="hindID")
    x <- merge(x,hind[,.(hindID,start,end)],by.x="preyID",by.y="hindID",suffixes=c(".bait",".prey"))
    x
}

   ## AHR
b2g[gene=="AHR",]
b2g[gene=="PTGER4",]

## IRF8
b2g[gene=="IRF8",]
b2g[gene=="STAT4",]
b2g[gene=="RORC",]
b2g[gene=="INPP4B",]
b2g[gene=="BACH2",]
b2g[gene=="IL6ST",]
b2g[,off:=(baitStart+baitEnd)/2]
b2g[gene %in% c("IL6ST","ANKRD55"),]

oia <- 55530266 - 55289478
oia

TODO <- list(## "AHR"= list(b=666461,NAME="AHR",cell="act",
              ##            zooms= list(prom=c(-1e+6,1e+6) #,
              ##                           #farleft=c(-5.1e+6,-4.1e+6),
              ##                           #farright= c(3.5e+6,4.5e+6)
              ##                        )),
             ##             "BACH2" = list(b=636904,NAME="BACH2",cell="act",
             ##                             zooms=list(all=c(-5e+5,5e+5))),
             "ETS1"=list(b=139778,NAME="ETS1",cell="act",
                         zooms=list(all=c(-1e+6,1e+6))),
             ##"GAPDH"=list(b=143494,NAME="GAPDH",cell="non",
             ##zooms=list(all=c(-1e+6,1e+6))),
             ##            "PTGER4" = list(b=567688,NAME="PTGER4",cell="act",
             ##                            zooms=list(all=c(-8e+5,2e+5))),
             ## "IL6ST" = list(b=571941,NAME="IL6ST",cell="act",
             ##                zooms=list(prox=c(-2e+5 - oia,2e+6 - oia))),
             "AC009505.2" = list(b=361840,NAME="AC0095052",cell="non",
                                  zooms=list(all=c(-6e+5,8e+5)))
             ## TNFAIP8 = list(b=592091,NAME="TNFAIP8",cell="non",zooms=list(all=c(-5e+5,5e+5))),
             ## CLVS1 = list(b=724747,NAME="CLVS1",cell="non"),
             ## FYN=list(b=643601,NAME="FYN",cell="non"),
             ## ## CD6=list(b=120463,NAME="CD6",cell="act",zooms=list(all=c(-2e+5,6e+5))),
             ## PTPN4=list(b=365957,NAME="PTPN4",cell="act",zooms=list(all=c(-2e+5,7e+5))),
             ## ## ARRDC3=list(b=583213,NAME="ARRDC3",cell="act"),
             ## PCF11=list(b=125742,NAME="PCF11",cell="act",zooms=list(all=c(-5e+5,5e+5))),
             ## FAM53B=list(b=100122,NAME="FAM53B",cell="act"),
             ## JADE2=list(b=596929,NAME="JADE2",cell="act"),
             ## BCL2=list(b=315353,NAME="BCL2",cell="act"),
             ## BOD1=list(b=608571,NAME="BOD1",cell="act")
             ## "ANKRD55" = list(b=572013,NAME="ANKRD55",cell="act",
             ##               zooms=list(prox=c(-2e+5,1e+6)))
             ## "INPP4B" = list(b=540390,NAME="INPP4B",
             ##                 zooms=list(proxwide=c(-2.2e+6,2e+5),
             ##                            proxnarrow=c(-1e+6,2e+5))),
             ##"STAT4" = list(b=388501,NAME="STAT4",
             ##               zooms=list(prox=c(-5e+4,2e+5))),
             ##             "RORC" = list(b=35775,NAME="RORC",
             ##                           zooms=list(prox=c(-1e+6,1e+6))),
             ##             "IRF8"=list(b=277910, NAME="IRF8", cell="act",
             ##                         zooms=list(prox=c(0,200000)))
)

##' NUMBERS
##' 
(load(file.path(d,"analysis","numbers.RData")))
z <- dblockdetail
uz <- z[,.(length=unique(block_length),sum.mppc=sum(mppc)),by=c("baitID","block","expt")]
uz <- merge(uz,b2g[,.(baitID,gene,biotype)],by="baitID",allow.cartesian = TRUE)
head(uz)
summary(uz$length)
uz <- uz[order(length,decreasing = TRUE),]
uz[length>40,]
uz[gene=="ETS1" & expt=="Promoter, act", .(sum(length), sum(sum.mppc))]
uz[gene=="AC009505.2" & expt=="Promoter, non", .(sum(length), sum(sum.mppc))]
subset(uz,gene=="AC009505.2")

## number of fragments per cluster for ETS1
CELL <- "act"
g <- "ETS1"
## (load(file.path(d,paste0(CELL,"prom","-joined.RData"))))
(load(file.path(d,"analysis",paste0(substr(CELL,1,1),"CD4pchic.RData"))))
## add residuals

b <- TODO[[g]]$b
NAME <- TODO[[g]]$NAME
             
OFFSET <- offset <- unique((b2g[baitID==b,]$baitStart + b2g[baitID==b,]$baitEnd)/2) 
CHR <- chr <- as.character(unique(b2g[baitID==b,]$baitChr))
xb <- xbi <- data[baitID==b,]
##xb <- merge(xb,cs[expt==paste0("Promoter, ",CELL),],by=c("baitID","preyID"),all.x=TRUE)
xb <- addstartend(xb)
REF <- gr <- GRanges(seq=Rle(rep(chr,nrow(xb))),
                     IRanges(start=xb$dist+offset,width=4000))
FROM=min(xb$dist)+offset
TO=max(xb$dist)+offset

#genomics <- load.genomics(CELL)
#tracks <- lapply(genomics,make.dt)
## rna <- load.stranded(CELL,CHR)
## peaks <- load.peaks(CELL,CHR,FROM,TO)
genes <- getgenes(CHR,FROM,TO,OFFSET,AC009505.2 = g=="AC009505.2")
## chromhmm <- load.chromhmm(CELL,CHR,FROM,TO)

xb[,cumc:=cumsum(mppc)]
ggplot(xb,aes(x=mid.prey,y=cumc)) + geom_path() + geom_path(aes(y=sqrt(mppc)*7),col="red") + xlim(1.275e+8,1.29e+8) + background_grid() + geom_hline(yintercept=1:7)
xb[mid.prey>1.28e+8 & mid.prey<1.285e+08,]$cumc
4.35305 - 2.43495

## if(!file.exists(f)) {
##     plotone(xb)
##     ggsave(file.path(d,"../figures",paste0("ex-",NAME,"-",CELL,".pdf")),height=h,width=w)
## }

library(data.table)

## GLOBAL VARIABLES THAT WILL BE USED BY FUNCTIONS IN viz.R
CELL="non"; g <- "AC009505.2"; i <- 1
h <- 12; w <- 10
for(CELL in c("non","act")) {
    message("CELL: ",CELL)
    (load(file.path(d,"analysis",paste0(substr(CELL,1,1),"CD4pchic.RData"))))
    res <- fread(file.path("/home/cew54/bsu/peaky/supp",paste0("suppdata-cd4-",substr(CELL,1,1),"CD4pchic.csv.gz")))
    ## (load(file.path(d,paste0(CELL,"prom","-joined.RData"))))
    data <- merge(data,res[,.(baitID,preyID,residual)],by=c("baitID","preyID"))
    for(g in names(TODO)) {
        message(CELL,"\t",g)
        if(TODO[[g]]$cell!=CELL)
            next
        b <- TODO[[g]]$b
        NAME <- TODO[[g]]$NAME
        f <- file.path(d,"figures",paste0("ex-",NAME,"-",CELL,".pdf"))
        fz <- if("zooms" %in% names(TODO[[g]])) {
                  file.path(d,"figures",paste0("ex-",NAME,"-",CELL,"-zoom",names(TODO[[g]]$zoom),".pdf"))
              } else {
                  NULL
              }
        ## if(all(file.exists(c(f,fz))))
        ##     next
        
        OFFSET <- offset <- unique((b2g[baitID==b,]$baitStart + b2g[baitID==b,]$baitEnd)/2) 
        CHR <- chr <- as.character(unique(b2g[baitID==b,]$baitChr))
        xb <- xbi <- data[baitID==b,]
        ##xb <- merge(xb,cs[expt==paste0("Promoter, ",CELL),],by=c("baitID","preyID"),all.x=TRUE)
        xb <- addstartend(xb)
        REF <- gr <- GRanges(seq=Rle(rep(chr,nrow(xb))),
                             IRanges(start=xb$dist+offset,width=4000))
        FROM=min(xb$dist)+offset
        TO=max(xb$dist)+offset
        
        genomics <- load.genomics(CELL)
        tracks <- lapply(genomics,make.dt)
        ## rna <- load.stranded(CELL,CHR)
        ## peaks <- load.peaks(CELL,CHR,FROM,TO)
        genes <- getgenes(CHR,FROM,TO,OFFSET,AC009505.2 = g=="AC009505.2")
        chromhmm <- load.chromhmm(CELL,CHR,FROM,TO)
        
        ## if(!file.exists(f)) {
        ##     plotone(xb)
        ##     ggsave(file.path(d,"../figures",paste0("ex-",NAME,"-",CELL,".pdf")),height=h,width=w)
        ## }
        
        if("zooms" %in% names(TODO[[g]])) {
            for(i in seq_along(TODO[[g]]$zooms)) {
                ## if(file.exists(fz[i]))
                ##     next
                plotone(xb,zoom=TODO[[g]]$zooms[[i]])
                ggsave(fz[i], height=h,width=w)
            }
        } else {
            ## if(!file.exists(f)) {
                plotone(xb)
                ggsave(file.path(d,"figures",paste0("ex-",NAME,"-",CELL,".pdf")),height=h,width=w)
            ## }
        }
    }
}




   ##                  ## ac27
    ## ggplot(as.data.frame(binned$ac27), aes(xmin=start-offset,xmax=end-offset,ymin=0,ymax=(V1))) + geom_rect(col="cornflowerblue",fill="cornflowerblue") + ggtitle("H3K27ac") + xlim(from,to),
    ##                  ## me1
    ## ggplot(as.data.frame(binned$me1), aes(xmin=start-offset,xmax=end-offset,ymin=0,ymax=(V1))) + geom_rect(col="cornflowerblue",fill="cornflowerblue") + ggtitle("H3K4me1") + xlim(from,to),
    ##                  ## me3
    ## ggplot(as.data.frame(binned$me3), aes(xmin=start-offset,xmax=end-offset,ymin=0,ymax=(V1))) + geom_rect(col="cornflowerblue",fill="cornflowerblue") + ggtitle("H3K4me3") + xlim(from,to),
    ## rna
                                        #ggplot(rna.this,aes(xmin=start-offset,xmax=end-offset,ymax=y)) + geom_rect(ymin=0,col="darkblue",fill="darkblue") + xlim(from,to) + ggtitle("RNAseq (log2")
                  ## ggplot(as.data.frame(binned$rna), aes(xmin=start-offset,xmax=end-offset)) +
                  ##  geom_rect(aes(ymax=log2(V1)),ymin=0,col="darkblue",fill="darkblue")+ xlim(from,to) + ggtitle("RNAseq (log2)")
                  ## geom_rect(aes(ymax=(plus)),ymin=0,col="red",fill="red") +
                  ## geom_rect(aes(ymin=(minus)),ymax=0,col="darkblue",fill="darkblue") + xlim(from,to)
               

## strand(ac27.this) <- "*"
## strand(rna.this) <- "*"
## ac27track <- DataTrack(ac27binned,genome="hg19",chromosome=chr,name="H3K27ac",
##                       type="l",window=-1)
## rnatrack <- DataTrack(rnabinned,genome="hg19",chromosome=chr,name="RNAseq",
##                       type="l",window=-1,col=c("cornflowerblue","purple"),groups=c("+","-"))

## tracks <- list(a=axis,
##                g=genes,
##                n=DataTrack(gr,data=log(xb$N+1),name="Log (Counts + 1)"),
##                r=DataTrack(gr,data=xb$residual,name="residuals",baseline=c(-2,0,2),col.baseline="gray"),
##                c=DataTrack(gr,data=xb$chicago,name="CHiCAGO",baseline=5,col.baseline="gray"),
##                p=DataTrack(gr,data=xb$mppc,name="MPPC",type="l"),
##                rna=rnatrack 
##               ,ac27=ac27track
##                )
## from=min(xb$dist)+offset;to=max(xb$dist)+offset
## plotTracks(tracks,chromosome=chr,from=from,to=to)



## rlim=c(-12, 12)
## clim=220
## ac27 <- chpt(file.path(seqpath,"H3K27ac_PoolQTW_Act.bam"), "H3K27ac", ylim=clim, "green4")
## plotTracks(ac27,chromosome="chr7",from=1e+6,to=10e+6)

## from,to=to)


## rtrack
## plotTracks(ntrack)


## ?DataTrack
## data(twoGroups)
## head(twoGroups)
## dTrack <- DataTrack(twoGroups, name = "uniform")
## plotTracks(dTrack)

## rnat <- function(file, name, ylim) {
##     logt <- function(x) sign(x) * log2(abs(x)+0.1)
## #    kk <- strandedBamImport("data_srt.2_non.bam")
##     d <- DataTrack(range=file, importFunction=strandedBamImport, stream=T, baseline=0, col.baseline="grey",
##               col=c("#e31a1c","#1f78b4"), groups=c("Forward", "Reverse"), type=c("histogram", "g"), col.histogram=NA,
##               name=name, ylim=ylim, transformation=logt)
##     displayPars(d)$col.grid <- adjustcolor("grey", alpha.f = 0.2)
##     displayPars(d)$v <- 0
##     d
## }
## chpt <- function(file, name, ylim, col) {
##                                         #    AlignmentsTrack(file, isPaired=F, type=c("mountain"), name=name, ylim=c(0, ylim))
##     d <-
##     ## AlignmentsTrack(file, isPaired=F, type=c("coverage"), name=name, ylim=c(0, ylim),
##     DataTrack(file, isPaired=F, type=c("polygon"), window=-1, windowSize=1000, baseline=0, name=name, ylim=c(0, ylim)
##              ## ,transformation=function(x){x[x<10] <- NA; x}
##               )
##     displayPars(d)$col=col
##     d
## }
## regt <- function(file, name, chromosome) {
##     d <- read.table(file)
##     colnames(d) <- c("chromosome", "start", "end", "state")
##     d$start <- d$start+1
##     s <- d[d$chromosome==chromosome,]
##     txn <- NA
##     enhancer <- "#cab2d6"
##     promoter <- "#cab2d6"
##     rep <- NA
##     AnnotationTrack(chromosome=s$chromosome, start=s$start, end=s$end, feature=s$state,
##                     name=name, genome="ENSEMBL", stacking="dense", stackHeight=0.15,
##                     E1=txn, E2=txn, E3=txn, E4=enhancer, E5=enhancer,
##                     E6=enhancer, E7=promoter, E8=promoter, E9=promoter, E10=promoter,
##                     E11=enhancer, E12=rep, E13=rep, E14=rep, E15=rep,
##                     col=NA)
## }
## rgnt <- function(file, name, chromosome) {
##     d <- read.table(file)
##     colnames(d) <- c("chromosome", "start", "end")
##     d$start <- d$start+1
##     s <- d[d$chromosome==chromosome,]
##     enhancer <- "#cab2d6"
##     AnnotationTrack(chromosome=s$chromosome, start=s$start, end=s$end,
##                     name=name, genome="ENSEMBL", stacking="dense", stackHeight=0.15)
## }
## ant <- function(range) AnnotationTrack(range=range, genome="ENSEMBL", name="GWAS", stackHeight=0.3, col=NA)
                                        #genomics <- load.genomics("act")
                                        #binned <- lapply(genomics,binner,binsize=1000)
                                        #binned <- lapply(binned, trimzeros)
