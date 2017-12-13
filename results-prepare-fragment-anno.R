library(data.table)

d <- "/mrc-bsu/scratch/cew54/peaky/summary/tables"

## positions
require(rtracklayer)
f <- file.path(CD4CHIC.DATA, "Digest_Human_HindIII.bed")
hind <- import(f)
seqlevels(hind) <- paste0("chr",seqlevels(hind))
hind <- keepStandardChromosomes(hind)
hind <- as.data.table(hind)
hind[,strand:=NULL]
hind[,width:=NULL]
setnames(hind,"name","hindID")
setkey(hind,hindID)
save(hind,file=file.path(d,"../derived","hind-position.RData"))


## bait2gene
int.genes <- fread(file.path(CD4CHIC.DATA,"merged_samples_12Apr2015_full_denorm_bait2baits_e75.tab"),select=c(1:8))
setnames(int.genes, names(int.genes), c("id","gene","biotype","strand","baitChr","baitStart","baitEnd","baitID"))
setkey(int.genes,baitID,id)
b2g <- unique( int.genes )
save(b2g,file=file.path(d,"../derived/b2g.RData"))

## chromatin data
require(GenomicRanges)
require(S4Vectors)
require(magrittr)
(load(file.path(CD4CHIC.OUT,"chipseq","hind-our-peaks.RData")))
    ##     hind <- get.hind()
    hind <- mcols(hind)[,-2] %>% as.data.frame() %>% as.data.table()
hind$name %<>% as.integer()
cat(names(hind),sep=",")

chromatin <- hind[,.(name,H3K27ac_Act,H3K27ac_NAct,H3K4me3_Act,H3K4me3_NAct,H3K27me3_NAct,H3K27me3_Act,H3K36me3_Act,H3K36me3_NAct,H3K4me1_Act,H3K4me1_NAct,H3K9me3_Act,H3K9me3_NAct)]
setnames(chromatin,"name","hindID")
setkey(chromatin,hindID)
save(chromatin,file=file.path(d,"../derived/chromatin-peaks.RData"))


## chromhmm
f <- function(n) {
 (load(file.path(CD4CHIC.OUT,"chipseq",paste0("hind-",n,".RData"))))
        hind.dt <- mcols(hind) %>% as.data.table()
        hind.dt[,bait:=NULL]
        hind.dt$name %<>% as.integer()
 setnames(hind.dt,"name","hindID")
 setkey(hind.dt,hindID)
        return(hind.dt)
     }
    chromhmm.act <- f("Act_15")
chromhmm.non <- f("NAct_15")
save(chromhmm.act,file=file.path(d,"../derived/chromhmm-act.RData"))
save(chromhmm.non,file=file.path(d,"../derived/chromhmm-non.RData"))

## exprs
erna <- fread(file.path(CD4CHIC.DATA,"rna-with-regrna-diff-v2.csv"))
setnames(erna,c("logFC","adj.P.Val","1_act","2_act","1_non","2_non"),
         c("logFC.erna","FDR.erna","act_1","act_2","non_1","non_2"))
##    erna <- erna[type %in% c("regulatory","lincRNA","protein_coding","pseudogene"), ]
    intergenic.ids <- scan(file.path(CD4CHIC.DATA, "distant-erna.csv"),what="")
    erna[,intergenic:=id %in% intergenic.ids]
    erna[id %in% intergenic.ids,type:="intergenic.reg"]
rnaseq <- erna

save(rnaseq,file=file.path(d,"../derived/rnaseq.RData"))

## interactions
f <- function(x) {
    setnames(x,"oeID","preyID")
    dbl <- intersect(x$baitID,x$preyID)
    x2 <- x[preyID %in% baitID,]
    setnames(x2,c("baitID","preyID"),c("preyID","baitID"))
    rbind(x,x2[,names(x),with=FALSE])
}


## not good enough - need full chicago results.
## see files in /home/cew54/scratch/cd4chic/non-bybait
## column 5 =  bait, column 3 = prey (or vice versa)
## last column (29) is chicago score
## should extract B as well if possible
## for other column ids, see https://bioconductor.org/packages/release/bioc/vignettes/Chicago/inst/doc/Chicago.html

## or, for promoter, get from the peaky objects??

prom <- fread(file.path(CD4CHIC.DATA,"merged_samples_12Apr2015_full.txt"))
prom <- prom[,.(baitID,oeID,Total_CD4_Activated,Total_CD4_NonActivated)]
setnames(prom,c("Total_CD4_Activated","Total_CD4_NonActivated"),c("chic.act","chic.non"))
chicago <- f(prom)
save(chicago,file=file.path(d,"../derived/chicago-prom.RData"))

val <- fread(file.path(CD4CHIC.DATA,"validation_peakmatrix","peakMatrix_validation_cutoff0.txt"))
head(val)
val <- val[,.(baitID,oeID,aCD4val_merge,uCD4val_merge)]
setnames(val,c("aCD4val_merge","uCD4val_merge"),c("chic.act","chic.non"))
chicago <- f(val)
save(chicago,file=file.path(d,"../derived/chicago-val.RData"))

length(unique(val$baitID))
length(unique(val$oeID))
length(unique(prom$baitID))
length(unique(prom$oeID))

source("~/Projects/cd4chic/activation-analyses/R/common.R")
int.genes <- get.b2gene()
head(int.genes)
save(int.genes,file=file.path(d,"b2gene.RData"))

## setkey(int,baitID)
##     setkey(int.genes,baitID)
##     int <- merge(int,int.genes[,.(id,gene,biotype,strand,baitID)],allow.cartesian=TRUE)
##     ## ## subset to protein coding only
##     ## int <- int[biotype=="protein_coding",]
##     ## int[ oeID==120852 & baitID==120843, .(baitChr,baitStart,baitID,oeChr,oeStart,oeID,id,gene,Total_CD4_Activated,Total_CD4_NonActivated)]
##     int <- int[,b2b:= oeID %in% int$baitID]
##     int[,Background:=apply(int[,.(Erythroblasts,Megakaryocytes)],1,max)]
##     if(!is.null(threshold))
##         int <- int[pmax(Background,Total_CD4_Activated,Total_CD4_NonActivated)>=threshold, ]
##     int <- int[, c("baitChr","oeChr") := list(paste0("chr",baitChr), paste0("chr",oeChr)) ]
##     int <- int[, c("baitLength","oeLength"):= list(abs(baitStart-baitEnd), abs(oeStart-oeEnd)) ]
##     int <- int[,.(baitChr,baitStart,baitID,oeChr,oeStart,oeID,id,gene,Total_CD4_Activated,Total_CD4_NonActivated,Background,biotype,baitLength,oeLength,b2b)]
##     return(unique(int,by=NULL)) # sometimes one bait maps to multiple promoters for the same gene
## }
