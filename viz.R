library(Gviz)
library(GenomicRanges)
#library(GenomicInteractions)
d2 <- (file.path(CD4CHIC.ROOT,"seq/track"))
#source(file.path(d2,"stranded.r"))
options(ucscChromosomeNames=FALSE)

library(biomaRt)
e75.genemart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",  host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

axis <- GenomeAxisTrack(col="black", fontcolor="black")
seqpath <- "/scratch/cew54/cd4chic/chris"
library(Repitools)
library(Rsamtools)
library(GenomicAlignments)


bamreader <- function(path,ref,strand="*") {
    bam_fields <- c("pos", "qwidth","strand") #I believe these are the minimal parameters required
    wh <- ScanBamParam(what=bam_fields,which=ref,
                       flag=scanBamFlag(isNotPassingQualityControls=FALSE,
                                        isUnmappedQuery=FALSE,isDuplicate=FALSE ))
    aligns <- readGAlignments(path, param = wh)
    aligns.info <- values(aligns)
    aligns.gr <- as(aligns, "GRanges")
    values(aligns.gr) <- aligns.info
    if(strand=="*") {
        gr1 <- bindAsGRanges(coverage(aligns.gr))
        return(gr1[with(mcols(gr1), V1>=1)])
    } else {
        grp <- bindAsGRanges(coverage(aligns.gr[ strand(aligns.gr)=="+" ]))
        grn <- bindAsGRanges(coverage(aligns.gr[ strand(aligns.gr)=="-" ]))
        return(list(pos=grp[with(mcols(grp), V1>=1)],
                    neg=grn[with(mcols(grn), V1>=1)]))
    }
}

load.genomics <- function(what=c("act","non")) {
    what <- match.arg(what)
    if(what=="act") {
        atacfile <- "/scratch/wallace/jonny/Sample_7_Cutler_S1_L001_001-trimmed-pair1_mapd_nomt_sort_nodup.bam"
        rnafile <- file.path(seqpath,"data_srt.2_act.bam")
        ac27file <- file.path(seqpath,"H3K27ac_PoolQTW_Act.bam")
        me1file <- file.path(seqpath,"H3K4me1_PoolQTW_Act_1.bam")
        me3file <- file.path(seqpath,"H3K4me3_PoolQTW_Act_1.bam")
    } else {
        #atacfile <- "/scratch/wallace/jonny/Sample_8_Cutler_S2_L001_001-trimmed-pair1_mapd_nomt_sort_nodup.bam"
        rnafile <- file.path(seqpath,"data_srt.2_non.bam")
        ac27file <- file.path(seqpath,"H3K27ac_PoolQTW_NAct.bam")
        me1file <- file.path(seqpath,"H3K4me1_PoolQTW_NAct_1.bam")
        me3file <- file.path(seqpath,"H3K4me3_PoolQTW_NAct_1.bam")
    }
    #atac <- bamreader(atacfile,ref=REF)
    ac27gr <- bamreader(ac27file,ref=REF)
    me1gr <- bamreader(me1file,ref=REF)
    me3gr <- bamreader(me3file,ref=REF)
    ## srnagr <- bamreader(rnafile,strand="+-")
    rnagr <- bamreader(rnafile,ref=REF)
    return(list(ac27=ac27gr,me1=me1gr,me3=me3gr,rna=rnagr))
}

load.stranded <- function(what=c("act","non"),chr) {
    what <- match.arg(what)
    if(what=="act") {
        rnafile <- file.path(seqpath,"data_srt.2_act.bam")
    } else {
        rnafile <- file.path(seqpath,"data_srt.2_non.bam")
    }
    ## srnagr <- bamreader(rnafile,strand="+-")
    ## rna <- bamreader(rnafile,ref=REF,strand="+-")
    ## z <- coverage(rna$pos,weight="V1")[[chr]]
    ## n <- length(z@values)
    ## p <- data.table(start=c(1,cumsum(z@lengths)-1)[1:n],
    ##                 end=cumsum(z@lengths),
    ##                 y=z@values)
    ## z <- coverage(rna$neg,weight="V1")[[chr]]
    ## n <- length(z@values)
    ## m <- data.table(start=c(1,cumsum(z@lengths)-1)[1:n],
    ##                 end=cumsum(z@lengths),
    ##                 y=-z@values)
    ## rna <- rbind(p,m)
    rna <- bamreader(rnafile,ref=REF)
    z <- coverage(rna,weight="V1")[[chr]]
    n <- length(z@values)
    rna <- data.table(start=c(1,cumsum(z@lengths)-1)[1:n],
                      end=cumsum(z@lengths),
                      y=z@values)
    rna[y!=0,]
}
load.chromhmm <- function(what=c("act","non"),CHR,FROM,TO) {
    what <- match.arg(what)
    if(what=="act") {
        f <- file.path(seqpath,"../chromhmm","Act_15_segments.bed")
    } else {
        f <- file.path(seqpath,"../chromhmm","NAct_15_segments.bed")
    }
    d <- fread(f)
    setnames(d,c("chr","start","end","state"))
    use <- d$chr==CHR & d$start < TO & d$end > FROM &
    d$state %in% c("E4","E5","E6","E7","E8","E9","E10","E11")
    d[use,]
}
             
make.rnadt <- function(gr1,CHR,FROM,TO) {
    z <- coverage(gr1[ strand(gr1)=="+" ],weight="V1")[[CHR]]
    n <- length(z@values)
    p <- data.table(start=c(1,cumsum(z@lengths)-1)[1:n],
                    end=cumsum(z@lengths),
                    y=z@values)
    z <- coverage(gr1[ strand(gr1)=="-" ],weight="V1")[[CHR]]
    n <- length(z@values)
    m <- data.table(start=c(1,cumsum(z@lengths)-1)[1:n],
                    end=cumsum(z@lengths),
                    y=-z@values)
    rna <- rbind(p,m)
    rna[end > FROM & start < TO,]
}
make.dt <- function(gr1) {
    z <- coverage(gr1,weight="V1")[[CHR]]
    n <- length(z@values)
    p <- data.table(start=c(1,cumsum(z@lengths)-1)[1:n],
                    end=cumsum(z@lengths),
                    y=z@values)
    p[end > FROM & start < TO,]
}


## todo binnedaverage
## http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/GenomicRanges/html/genomicvars.html

binner <- function(gr1, binrange=NULL, nbins=1000,
                   binsize=(max(end(binrange)) - min(start(binrange)))/nbins) {
    sumPerBin <- function(x, binsize, mcolnames=NULL) {
        if (!is(x, "GenomicRanges"))
            stop("'x' must be a GenomicRanges object")
        if (any(is.na(seqlengths(x))))
            stop("'seqlengths(x)' contains NAs")
        bins <- IRangesList(lapply(seqlengths(x),
                                   function(seqlen)
                                       IRanges(breakInChunks(seqlen, binsize))))
        ans <- as(bins, "GRanges")
        seqinfo(ans) <- seqinfo(x)
        if (is.null(mcolnames))
            return(ans)
        sumMCol <- function(colname) {
            cvg <- coverage(x, weight=colname)
            views_list <- RleViewsList(
                           lapply(names(cvg),
                                  function(seqname)
                                      Views(cvg[[seqname]], bins[[seqname]])))
            unlist(viewMeans(views_list)*binsize, use.names=FALSE)
        }
        mcols(ans) <- DataFrame(lapply(mcols(x)[mcolnames], sumMCol))
        #use <- mcols(ans)[[mcolnames]] > 0
                                        #ans[ use ]
        ans
    }
    sumPerBin(gr1,round(binsize),"V1")
}

trimzeros <- function(x,cols) {
    use <- rowSums(abs(as.data.frame(mcols(x)))) > 0
    x[use ]
}

##ac27binned <- trimzeros(binner(genomics$ac27,binsize=1000))
## tmp <- lapply(rnagr,binner,binsize=1000)
## rnabinned <- tmp[[1]]
## names(mcols(rnabinned)) <- "plus"
## rnabinned$minus <- -mcols(tmp[[2]])[["V1"]]
## rnabinned <- trimzeros(rnabinned)

## tracks for region with gviz
gviz <- function(what=c("act","non"),chr,from,to) {
    

    what <- match.arg(what)
    if(what=="act") {
        ac27file <- file.path(seqpath,"../chip","H3K27ac_PoolQTW_Act_1_macs_peaks.narrowPeak")
        me1file <- file.path(seqpath,"../chip","H3K4me1_PoolRG_Act_2_macs_peaks.gappedPeak")
        me3file <- file.path(seqpath,"../chip","H3K4me3_PoolQTW_Act_1_macs_peaks.narrowPeak")
    } else {
        ac27file <- file.path(seqpath,"../chip","H3K27ac_PoolQTW_NAct_1_macs_peaks.narrowPeak")
        me1file <- file.path(seqpath,"../chip","H3K4me1_PoolRG_NAct_2_macs_peaks.gappedPeak")
        me3file <- file.path(seqpath,"../chip","H3K4me3_PoolQTW_NAct_1_macs_peaks.narrowPeak")
    }
}


## plot with ggplot

getgenes <- function(chr,from,to,offset,AC009505.2=FALSE) {
    genes1 <- BiomartGeneRegionTrack(genome="hg19", name="Genes", transcriptAnnotation="symbol",
                                mart=e75.genemart,
                                collapseTranscripts="meta",
                                filters=list(with_ox_refseq_mrna=T),
                                chromosome=chr,from=from,to=to)
    genes2 <- BiomartGeneRegionTrack(genome="hg19", name="Genes", transcriptAnnotation="symbol",
                                mart=e75.genemart,
                                collapseTranscripts="meta",
                                filters=list(with_ox_refseq_ncrna=TRUE),
                                chromosome=chr,from=from,to=to)
     ## genes3 <- BiomartGeneRegionTrack(genome="hg19", name="Genes", transcriptAnnotation="symbol",
     ##                            mart=e75.genemart,
     ##                            collapseTranscripts="meta",
     ##                            filters=list(symbol="AC009505.2"),
     ##                            chromosome=chr,from=from,to=to)
     genes.df <- rbind(as.data.frame(genes1@range),as.data.frame(genes2@range))
    exons.dt <- as.data.table(genes.df)
    exons.dt <- exons.dt[ start > from & end < to,]
    exons.dt[,start:=start-offset]
    exons.dt[,end:=end-offset]
    exons.dt[,y:=ifelse(strand=="+",0.2,-0.2)]
    genes.dt <- exons.dt[,.(start=min(start),end=max(end),strand=unique(strand),y=unique(y)),
                         by=c("gene","symbol")]
    genes.dt[strand=="-",c("start","end"):=list(end,start)]
    if(AC009505.2)
        genes.dt <- rbind(genes.dt,
                          data.table(gene="ENSG00000235522",
                                     symbol = "AC009505.2", 
                                     start = 773.5,
                                     end = -9869.5, strand = "-",
                                     y=-0.2))
    return(list(exons=exons.dt,genes=genes.dt))
}
## genes <- getgenes(CHR,FROM,TO,OFFSET)


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @param what "act" or "non"
##' @param chr
##' @param from
##' @param to
##' @return 
##' @author Chris Wallace
load.peaks <- function(what=c("act","non"),chr,from,to) {
    what <- match.arg(what)
      if(what=="act") {
        ac27file <- file.path(seqpath,"../chip","H3K27ac_PoolQTW_Act_1_macs_peaks.narrowPeak")
        me1file <- file.path(seqpath,"../chip","H3K4me1_PoolRG_Act_2_macs_peaks.gappedPeak")
        me3file <- file.path(seqpath,"../chip","H3K4me3_PoolQTW_Act_1_macs_peaks.narrowPeak")
    } else {
        ac27file <- file.path(seqpath,"../chip","H3K27ac_PoolQTW_NAct_1_macs_peaks.narrowPeak")
        me1file <- file.path(seqpath,"../chip","H3K4me1_PoolRG_NAct_2_macs_peaks.gappedPeak")
        me3file <- file.path(seqpath,"../chip","H3K4me3_PoolQTW_NAct_1_macs_peaks.narrowPeak")
    }
    ac27 <- fread(ac27file,skip=1,select=c(1:3,5))
    setnames(ac27,c("chr","start","end","score"))
    use <- ac27$chr==chr & ac27$start <= to & ac27$end >= from
    ac27 <- ac27[use,]
   me3 <- fread(me3file,skip=1,select=c(1:3,5))
    setnames(me3,c("chr","start","end","score"))
    use <- me3$chr==chr & me3$start <= to & me3$end >= from
    me3 <- me3[use,]
    me1 <- fread(me1file,skip=1,select=c(1,5,7,8))
    setnames(me1,c("chr","score","start","end"))
    use <- me1$chr==chr & me1$start <= to & me1$end >= from
    me1 <- me1[use,]
    ac27$mark <- "H3K27ac"
    me1$mark <- "H3K4me1"
    me3$mark <- "H3K4me3"
    rbind(ac27,me3,me1[,names(ac27),with=FALSE])
}

asinh_trans <- function() { scales:::trans_new("asinh",asinh,sinh,domain=c(0,Inf)) }
kb_trans <- function() {
    trans <- x /1000
    inv <- x*1000
    scales:::trans_new("kb",trans, inv,domain=c(-Inf,Inf))
}

    

                                        #    rnan <- rnat(file.path(seqpath,"data_srt.2_non.bam"), "RNA-seq (n)", rlim)
#    rnaa <- rnat(file.path(seqpath,"data_srt.2_act.bam"), "RNA-seq (a)", rlim)
#    ac27n <- chpt(file.path(seqpath,"H3K27ac_PoolQTW_NAct.bam"), "H3K27ac (n)", clim, "green4")
#    ac27a <- chpt(file.path(seqpath,"H3K27ac_PoolQTW_Act.bam"), "H3K27ac (a)", clim, "mediumpurple4")
#tracks$reg <- rgnt(file.path(seqpath,"reg.bed"), name="regRNA", chromosome)

#source("~/Projects/peaky/strandedBamImport.R")
#libType="fr-firststrand"
#rnatrack <- DataTrack(rnafile, genome="hg19", chromosome=chr, importFunction=strandedBamImport, stream=TRUE, legend=TRUE, col=c("cornflowerblue","purple"), groups=c("Forward","Reverse"))

#plotTracks(ac27track,chromosome="7",from=from,to=to)


## 1mm = 2.83465pt
pt <- 1/2.83465
plotone <- function(xb,zoom=NULL,pt=0.1) {
    if(!is.null(zoom))
        xb <- xb[dist>=zoom[[1]] & dist <= zoom[[2]], ]
    from <- min(xb$dist) 
    to <- max(xb$dist)
    xb[c(1,nrow(xb)),residual:=0]
    ## peaks.this <- peaks[end-offset > from & start-offset < to,]
    ## rna.this <- rna[end-offset > from & start-offset < to,]
    ## rna.this[,RNA.strand:=ifelse(y>0,"+","-")]
    tracks.this <- lapply(tracks, function(x) x[end-offset > from & start-offset < to, ])
    ## raw 
    yv <- with(genes$genes[ start>from & end<to,],
               rep(seq(0.03,0.23,by=0.05),length.out=length(y)))
    plots <- list(ggplot(xb,aes(x=dist,y=N+1)) + geom_point(size=pt,col="black",fill="grey30") +
                  ##geom_path(aes(y=B+T)) +
                  scale_y_log10() + ggtitle("Counts (log10 scale)") + ylim(0,max(xb$N+1)),
### resid
                  ggplot(xb,aes(x=(start.prey-offset),xend=(end.prey-offset),yend=residual,#x=dist,
                                y=residual)) + geom_point(col="black",fill="grey30",size=pt) + ggtitle("NB residual") + geom_hline(yintercept=0,col="grey"),
### CHiCAGO
                  ggplot(xb,aes(x=dist,
                                        #x=(start.prey-offset),xend=(end.prey-offset),yend=chicago,
                                ymin=0,ymax=chicago,
                                y=chicago)) +
                                        #geom_ribbon() +
                  geom_area(col="black",fill="grey30") +
                  geom_area(aes(y=pmin(5,chicago)),col="grey",fill="grey") +
                                        #                                        geom_point(size=pt) +
                  ggtitle("CHiCAGO score (asinh scale)") + geom_hline(yintercept=5,col="grey") +
                  scale_y_continuous(trans=asinh_trans(),breaks=c(0,5,20,40),limits=c(0,max(xb$chicago))),
### Credset
                  ## ggplot(xb,aes(xmin=(start.prey-offset),xmax=(end.prey-offset),ymin=0,ymax=h)) + geom_rect(fill="grey20",col="grey20") + ggtitle("MPPC finemap"),
### MPPC

                  ggplot() +
                  geom_rect(aes(xmin=(start-offset),xmax=(end-offset)),
                            ymin=0,ymax=max(xb$mppc)^0.5,fill="violet",col="violet",
                            data=chromhmm) +
                  ## geom_point(aes(x=(mid.prey-offset), y=mppc),
                  ##            size=pt,data=xb) +
                  geom_area(aes(x=(mid.prey-offset), y=mppc),
                            col="black",fill="grey30",
                            data=xb) +
                  ggtitle("MPPC (sqrt scale)") +
                  scale_y_sqrt()
                 ,
### bamfiles
                  ac27=ggplot(tracks.this$ac27,aes(xmin=start-offset,xmax=end-offset,ymax=y,ymin=0)) + geom_rect(col="black",fill="black") + ggtitle("K3K27ac"),
                  me1=ggplot(tracks.this$me1,aes(xmin=start-offset,xmax=end-offset,ymax=y,ymin=0)) + geom_rect(col="black",fill="black") + ggtitle("H3K4me1"),
                  me3=ggplot(tracks.this$me3,aes(xmin=start-offset,xmax=end-offset,ymax=y,ymin=0)) + geom_rect(col="black",fill="black") + ggtitle("H3K4me3"),
### CHROMHMM - B&W
                  ##ggplot(chromhmm, aes(xmin=(start-offset),xmax=(end-offset),ymin=0,ymax=1,fill=state,col=state)) + geom_rect()  +
                  ## chromhmm=ggplot(chromhmm, aes(xmin=(start-offset),xmax=(end-offset),ymin=0,ymax=1,fill=Chromatin)) + geom_rect(fill="grey20",col="grey20")  + ggtitle("active chromatin"),
### RNAseq

                  ## p.rna=ggplot(rna.this,aes(xmin=start-offset,xmax=end-offset,
                  ##                           ymax=sign(y) * log2(abs(y)),ymin=0)) +
                  ##                           ## fill=RNA.strand,col=RNA.strand)) +
                  ##   geom_rect() +
                  ##   theme(legend.position="none") +
                  ##   ggtitle("RNAseq (log2 counts)"), #+
                    ## scale_colour_manual("Genes (strand)",values=c("grey20","darkslateblue")) +
                    ## guides(colour=FALSE) +
                    ## scale_fill_manual("Genes (strand)",values=c("grey20","darkslateblue")) +
                    ## guides(fill=FALSE)
### genes
g=ggplot(genes$genes[ start>from & end<to,],
         aes(x=start,xend=end,col=strand))+
  geom_segment(aes(y=sign(y) * yv, #rep(seq(0,0.15,by=0.05),length.out=length(y)),
                   yend=sign(y) * yv), #rep(seq(0,0.15,by=0.05),length.out=length(y))),
               arrow = arrow(length = unit(0.1,"cm")),size=1) +
  geom_text(aes(x=start,y=sign(y) * yv + sign(y) * 0.02, #rep(seq(0.05,0.2,by=0.05),length.out=length(y)),
                label=symbol),hjust=1,size=12 * mm2pt) +
  ## geom_rect(aes(xmin=start,xmax=end,ymin=y-0.03, ymax=y+0.03,col=strand),
  ##           data=genes$exons[ start>from & end<to,]) +
  scale_y_continuous(limits=c(-0.3,0.3)) +
  scale_colour_manual("Genes (strand)",values=rev(c("grey20","darkslateblue"))) + guides(colour=FALSE)
)
    ## add bait locations
    bloc <- xb[preyID %in% baits,]$dist %>% unique()%>%setdiff(.,0)
    plots <- lapply(plots, function(p) {
        p + #geom_vline(xintercept=bloc,col="grey",linetype="dashed",alpha=0.5) +
                                        #background_grid() +
          geom_vline(xintercept=0,col="Firebrick1",linetype="dashed") +
          xlim(from,to) + 
          theme(#legend.position="top",
                axis.text.x=element_blank(),
                axis.title.x=element_blank(),
                axis.title.y=element_blank()
                )
    })
    plots[-c(1:4)] %<>% lapply(., function(p) {
        p + theme(axis.text.y=element_blank(), #scale_y_continuous(breaks=NULL) 
                  axis.line.y = element_blank(), 
                  axis.ticks.y = element_blank())
    })
    ## turn x axis back on for genes
    plots[[length(plots)]] <- plots[[length(plots)]] +
      scale_x_continuous("Distance from bait (bp)"## , trans="kb"
                        ,limits=c(from,to)) + theme(axis.text.x = element_text(),axis.title.x=element_text()) 
    ## make plot
    plot_grid(plotlist=plots,ncol=1,align="v",
              rel_heights=c(rep(1,length(plots)-1),2))
}
## if(interactive())
##     plotone(xb,zoom=TODO[[g]]$zoom[[i]])



