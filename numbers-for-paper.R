#' load stuff
library(data.table)
library(ggplot2)
library(magrittr)
library(cowplot)
library(stargazer)

## d <- "/scratch/wallace/peaky/derived"
d <- "/mrc-bsu/scratch/cew54/peaky/summary/derived"
source("~/Projects/peaky/common.R")

##' make data
if(FALSE) {
    (load(file.path(d,"corr.RData")))
    xc <- rbindlist(CORR)
    xc[,expt:=factor(title.ab[experiment],levels=title.ab[o])]
    ##' add posterior expectation of gamma
    DATA <- structure(vector("list",length(title)),names=names(title))
    for(nm in names(title)) {
        message(" !!! ",nm) 
        fnm <- file.path(d,paste0(nm,"-joined-fm.RData"))
        (load(fnm))
        DATA[[nm]] <- x
        DATA[[nm]]$experiment <- nm
    }
    intnames <- intersect(names(DATA[[1]]),names(DATA[[2]]))
    data <- lapply(DATA,"[",j=intnames,with=FALSE) %>% rbindlist()
    data <- data[order(experiment,baitID,dist),]
    setkey(data,baitID,experiment)
    gamma <- unique(data[maxchic>5,.(gamma_mean.c5=sum(ifelse(chicago>=5,mppc,0)),
                                     gamma_mean=sum(mppc),
                                     chicago.n5=sum(chicago>=5),
                                     chicago.n3=sum(chicago>=3 & chicago<5),
                                     chicago.n2=sum(chicago>=2 & chicago<5)
                                     ),by=c("baitID","experiment")])
    x <- merge(xc,gamma,by=c("baitID","experiment"))
    head(x)
    
    ##' add blocks
    data[,csig:=as.numeric(chicago>5),by=c("baitID","experiment")]
    data[,start:=c(0,diff(csig))==1 | (csig==1 & c(0,diff(preyID))!=1),by=c("baitID","experiment")]
                                        #data[,end:=c(diff(csig),0)==-1,by=c("baitID","experiment")]
                                        #data[,mid:=csig==1 & start==FALSE & end==FALSE,by=c("baitID","experiment")]
    data[,block:=cumsum(start) * (csig==1),by=c("baitID","experiment")]
    
    ## summary(data$block)
    ## data[block>80,]
    ## y <- data[baitID==624482 & experiment=="nonval",]
    ## y[block %in% c(40:45) ,]
    ## head(y)
    ## summary(y)
    ## ggplot(y,aes(x=preyID,y=block)) + geom_path() + geom_point(aes(col=chicago>5),size=5)
    
    data[block>0,block_length:=.N,by=c("experiment","baitID","block")]
    data[block>0,exp_gamma:=sum(mppc),by=c("experiment","baitID","block")]
    
    dblock <- unique(data[block>0,],by=c("experiment","baitID","block"))
    head(dblock)
    summary(dblock)
                                        #dblock[,expt:=factor(title.ab[experiment],levels=title.ab[o])]
    data[,expt:=factor(title.ab[experiment],levels=title.ab[o])]
    dblock[,expt:=factor(title.ab[experiment],levels=title.ab[o])]
    data[,expt2:=factor(title.pv[experiment])]
    dblockdetail <- data[block>0, .(baitID,preyID,N,residual,chicago,mppc,dist,block,block_length,exp_gamma,expt,expt2,y,z)]
    dblock[,expt2:=factor(title.pv[experiment])]
    x[,expt2:=factor(title.pv[experiment])]

    save(dblock,dblockdetail,x,xc,file=file.path(d,"numbers.RData"))
    save(data,file=file.path(d,"alldata.RData"))
} else {
    (load(file.path(d,"numbers.RData"))) 
    (load(file.path(d,"alldata.RData")))
}

##' TEXT: approx value of Fb
(load(file.path(d,"hind-position.RData")))
hind[,mid:=(start+end)/2]
s <- split(hind,hind$seqnames)
s1 <- s[[1]]
s1
s1[,n:=0]
for(i in as.integer(1:nrow(s1)))
    set(s1,i,"n",sum(abs(s1$mid - s1$mid[i])<5e+6))
s1
summary(s1$n)

##' TABLE: numbers of baits with corr > 0.75
tl <- xc[,.(cr75=sum(maxcor>0.75), tot=.N),by="experiment"]
tab <- xc[chicago.mx>=5,.(cr75=sum(maxcor>0.75), tot.c5=.N),by="experiment"]
tab[,pc75:=100*cr75/tot.c5]
tab[,tot:=tl$tot]
tab

latex(tab[,.(tot.c5,cr75,pc75)],
      file=file.path(d, "../ptables/tab-peaky-corr-btwn-chains.tex"),
      ## rgroup=c("All baits","CHiCAGO signif. baits"), n.rgroup=c(4,4),
      colhead=c("Total baits","n. $\\rho > 0.75$","\\% $\\rho > 0.75$"),
      rowname=title.ab,rowlabel="Experiment",
      where="!hbp",
      cdec=c(0,0,1),
      label="tab:corr-mppc",
      booktabs=TRUE,
      caption="Number and \\% of baits for which correlation between MPPC between two parallel runs exceeded 0.75.  The total baits is the number for which at least one prey fragment has a CHiCAGO score \\textgreater 5.")

##' FIGURE: Dn of cor btwn chains
xc5 <- xc[chicago.mx>5,]
ggplot(xc5, aes(x=maxcor,fill=expt)) + geom_histogram(col="grey") +
facet_wrap(~expt,scales="free_y") +
geom_vline(xintercept=0.75,linetype="dashed") +
theme(legend.position="none") +
labs(x="Correlation",y="Number of baits")
ggsave(file=file.path(d,"../figures/cor-btwn-chains.pdf"),height=6,width=6)
xc[,median(maxcor),by="expt"]## >=0.85


##' FIGURE: Distribution of n>5 per bait
library(ggplot2)
ggplot(x[chicago.n5>0,],aes(x=chicago.n5,fill=expt2,col=expt2)) +
geom_histogram(bins=150) +
facet_wrap(~expt2,scales="free_y") + scale_x_log10("Number of preys with CHiCAGO score > 5 per bait (log10 scale)") + scale_y_continuous("Number of baits") + theme(legend.position="none")
ggsave(file.path(d,"../figures/nchic5-per-bait.pdf"),height=4,width=8)

##' TABLE: Dn of n>5 per bait
res <- with(x[chicago.n5>0, ], tapply(chicago.n5, expt,summary))
latex(do.call("rbind",res),file=file.path(d,"../ptables/nchic5-per-bait.tex"),
      rowlabel="Experiment",
      booktabs=TRUE,
      where="!hbp",
      label="tab:nchic5-per-bait",
      caption="Summary of distribution of number of preys with CHiCAGO score \\textgreater 5 per bait, amongst those baits with at least one such prey")

################################################################################

##' bait-prey pair comparisons

sqrt_trans <- function() {
    trans <- function(x) sqrt(x)
    inv <- function(x) x^2
    scales:::trans_new("sqrt",trans,inv,domain=c(0,Inf))
}

library(scales)

dsub <- split(data,data$expt) %>%
lapply(., function(x) x[sample(1:nrow(x),min(100000,nrow(x))),]) %>%
rbindlist()
## dsub0 <- dsub[mppc>0.01 | chicago>0.01,]

dim(data)
dim(dsub)

cr <- dsub[,.(rho=cor(mppc,chicago,method="spearman"),
              p=cor.test(mppc,chicago,method="spearman")$p.value),
           by="expt"]
cr

max(dsub$mppc)
(max(data$chicago))
asinh(max(data$chicago))
lx <- c(0,1,3,5,10,30); asinh(lx)
ly <- c(0,0.2,0.5,1); sqrt(ly)
ggplot(dsub, aes(x=asinh(chicago),y=(mppc),col=expt)) +
geom_point(alpha=0.5,data=dsub) +
## geom_hex(bins=50) +
geom_smooth(colour="black") +
geom_text(aes(label=sprintf("rho==%4.2f",rho)),data=cr,x=0.3,y=0.8,parse=TRUE,col="black") +
guides(colour=FALSE) +
## geom_rect(xmin=0,ymin=0,xmax=asinh(0.01),ymax=sqrt(0.01),aes(fill=expt,col=expt)) +
scale_y_continuous("MPPC (sqrt scale)",trans=sqrt_trans(),
                   breaks=sqrt(ly),labels=ly) +
scale_x_continuous("CHiCAGO score (asinh scale)",trans=trans_new("asinh","asinh","sinh"),
                   breaks=asinh(lx),labels=lx) +
facet_wrap(~expt)  

ggsave(file.path(d,"../figures/chic-mppc.png"),height=6,width=6) #,res=300)
capt <- paste("\\caption{Higher MPPC correlates with higher CHiCAGO scores. Figure shows random sample of ",
          nrow(dsub)," of ",nrow(data)," bait-prey pairs.  Spearman's rank correlation ($\\rho$) is shown.",
         "}")
capt
cat(capt, file=file.path(d,"../figures/chic-mppc.tex"))

##' distance profiles

data[,cdist:=cut(abs(dist),quantile(abs(dist),seq(0,1,by=0.02)),include.lowest=TRUE),by="expt"]

q <- function(x,nm) {
    ret <- as.list(quantile(x,c(0.05,0.25,0.5,0.75,0.95)));
    names(ret) %<>% sub("%","",.) %>% paste(nm,.,sep=".")
    ret
}
ddata <- data[,c(d=median(abs(dist)),
                 q(chicago,"chic"),
                 q(mppc,"mppc")),
              by=c("expt2","cdist")]
ddata

mdata <- melt(ddata,id.vars=c("expt2","cdist","d"),
              measure.vars=list(c("chic.5","mppc.5"),
                                c("chic.25","mppc.25"), 
                                c("chic.50","mppc.50"), 
                                c("chic.75","mppc.75"), 
                                c("chic.95","mppc.95")),
              value.name=c("q5","q25","q50","q75","q95"))
mdata[,variable:=c("CHiCAGO","MPPC")[variable]]
head(mdata)

ggplot(mdata,aes(x=d/1e+6,y=q50,ymin=q5,ymax=q95,col=expt2)) + geom_linerange() + facet_grid(variable~expt2,scales="free_y") + geom_smooth(se=FALSE,col="grey30") + geom_smooth(aes(y=q25),col="grey30",se=FALSE,linetype="dashed") + geom_smooth(aes(y=q75),se=FALSE,linetype="dashed",col="grey30") +
  labs(x="Distance from bait (mb)",y="Statistic") +
  theme(legend.position="none")
ggsave(file.path(d,"../figures/chicago-mppc-decay-dist.pdf"),height=6,width=8)


################################################################################

##' mppc shapes in runs

## head(dblockdetail)
## y <- dblockdetail
## y[,x:=NULL]
## y[,i:=1:.N]
## y[,medi:=runif(.N)]
## y[,medi:=as.numeric(median(i)),by=c("expt","baitID","block")]
## y[,x:=i-medi,by=c("expt","baitID","block")]
## y[,ay:=mppc/sum(mppc),by=c("expt","baitID","block")]
## head(y)
## y[,group:=paste(expt,baitID,block)]

## y5 <- y[block_length==5,.(ay=mean(ay)),by=c("x")]
## ggplot(y[block_length==5,],aes(x=x,y=ay)) + geom_path(aes(group=group),alpha=0.1) + geom_path(data=y5,col="Firebrick1")

################################################################################

##' comparison in runs/regions

##' FIGURE: contacts per bait, c vs mppc
xc <- x[,cor(chicago.n5,gamma_mean,method="spearman"),by="expt2"]
xc[,rho:=paste0("rho=='",sprintf("%3.2f",V1),"'")]
xc
p1 <- ggplot(x,aes(x=chicago.n5,y=gamma_mean,col=expt2)) + geom_point(alpha=0.5) + geom_abline(linetype="dashed") + facet_wrap(~expt2) + geom_smooth(method="lm",se=FALSE) +
geom_text(aes(label=rho),data=xc,x=2,y=-1,col="black",parse=TRUE) +
scale_x_log10("Number of contacts called by CHiCAGO (log scale)") + scale_y_log10("Expected number of direct contacts\n(log scale)") + theme(legend.position="none")## +
## geom_density_2d()
p1
ggsave(file.path(d,"../figures/contacts-per-bait-chic-mppc.pdf"),height=4,width=8)

##' examine long runs
ss <- split(dblock,dblock$expt2)
slope <- lapply(ss, function(x) lm(exp_gamma ~ block_length - 1,data=x) %>%
                                summary() %>%
                                "[["("coefficients")) %>%
do.call("rbind",.) %>%
as.data.frame()
rownames(slope) <- slope$expt2 <- names(ss)
slope

##' FIGURE: c vs mppc in runs
p2 <- ggplot(dblock,aes(x=block_length,y=exp_gamma,col=expt2)) + geom_point() +
geom_smooth(se=FALSE,method="lm",formula="y~x-1") +
geom_text(aes(x=80,y=100*Estimate,label=sprintf("slope = %4.2f%%",100*Estimate)),data=slope) +
geom_abline(linetype="dashed") +
labs(x="Length of run",y="Expected number of direct contacts") +
facet_wrap(~expt2) +
theme(legend.position="none")
p2
ggsave(file.path(d,"../figures/contacts-per-run-chic-mppc.pdf"),height=4,width=8)

#plot_grid(p1,p2,nrow=2)

## ggplot(dblock,aes(x=block_length,y=exp_gamma/block_length,col=expt)) + geom_point(alpha=0.5) +
## facet_wrap(~expt) + scale_x_log10("CHiCAGO run length") + scale_y_continuous("posterior mean proportion of direct contacts within run")

## dblock[block_length>10,median(exp_gamma/block_length),by="expt"]
## dblock[block_length>10,mean(exp_gamma/block_length),by="expt"] ## ~ 7%
## dblock[block_length==10,median(exp_gamma),by="expt"]
## dblock[block_length==10,mean(exp_gamma),by="expt"] ## ~ 7%

#' above is expected number in any block, what number of baits are needed to capture signal (credible set style)?

y <- dblockdetail
head(y)
y <- y[order(expt,baitID,block,1-mppc),]
head(y[block_length>1,])
y[baitID==254 & block==5,]
y[,mppc.cov:=cumsum(mppc)/sum(mppc),by=c("expt2","expt","baitID","block")]
y[,mppc.el:=1:.N,by=c("expt2","expt","baitID","block")]
y[,group:=factor(paste(expt,baitID,block))]
head(y,2)

ggplot(y,aes(x=mppc.el/block_length,y=mppc.cov,col=expt2,group=group)) +
geom_path(alpha=0.02) +
geom_smooth(se=FALSE) +
facet_wrap(~expt2) +
theme(legend.position="none")

##' FIGURE: what propn of fragments are required to cover {thr} proportion of total mppc?
f <- function(thr) {
    yt <- unique(y[mppc.cov>=thr & block_length>2,],by=c("expt2","expt","baitID","block"))
    yt[,bl:=cut(block_length,c(2,10,50,100,150))]
    yt[,thr:=thr]
    yt[,.(expt2,expt,bl,thr,mppc.prop=mppc.el/block_length)]
}
propcov <- lapply(c(seq(0.1,0.9,0.1)),f) %>% rbindlist()
propcov[,thr:=thr*100]

ggplot(propcov,aes(x=bl,y=mppc.prop*100,col=expt2,fill=factor(thr))) +
geom_boxplot(outlier.size=0.5) +
scale_fill_grey("% summed MPPC recovered",start=0.9,end=0.1) +
facet_grid(expt2~.) +
labs(x="Length of run (binned)",y="% fragments included") +
guides(fill = guide_legend(nrow = 1),#fill = guide_legend(reverse=TRUE),
colour=FALSE) +
theme(legend.position="bottom") + background_grid()

ggsave(file.path(d,"../figures/run-coverage.pdf"),width=8,height=8)

f <- function(thr) {
    yt <- unique(y[mppc.el>=block_length * thr & block_length>2,],by=c("expt2","expt","baitID","block"))
    yt[,bl:=cut(block_length,c(2,10,50,100,150))]
    yt[,thr:=thr]
    yt[,.(expt2,expt,bl,thr,mppc.cov)]
}
propcov <- lapply(c(seq(0.1,0.9,0.1)),f) %>% rbindlist()
propcov[,thr:=thr*100]

ggplot(propcov,aes(x=bl,y=100*mppc.cov,col=expt2,fill=factor(thr))) +
geom_boxplot(outlier.size=0.5) +
scale_fill_grey("% fragments included",start=0.9,end=0.1) +
facet_grid(expt2~.) +
labs(x="Length of run (binned)",y="% summed MPPC recovered") +
guides(fill = guide_legend(nrow = 1),#fill = guide_legend(reverse=TRUE),
colour=FALSE) +
theme(legend.position="bottom") + background_grid()

ggsave(file.path(d,"../figures/run-coverage-sw.pdf"),width=8,height=8)















################################################################################

## junk below here


ggplot(propcov,aes(x=factor(thr),y=mppc.prop,fill=bl)) + geom_boxplot() 


ggplot(yt,aes(x=bl,y=mppc.el/block_length)) + geom_boxplot() 

geom_abline() +
geom_smooth(method="lm",se=FALSE)
ggplot(yt,aes(x=block_length,y=mppc.el)) + geom_hex() + geom_abline() +
geom_smooth(method="lm",se=FALSE)

dblock[block_length>100,]



## start of a run has cdiff=1
## mid run cdiff=0
## last in a run has cdiff=-1

summary(x$chicago.n5)

x[,chic.cut:=cut(chicago.n5,c(0,3,7,18,400))]
thr <- unique(x[,.(chic.cut)])
thr$min <- c(0,3,7,18)
thr$max <- c(3,7,18,350)

ggplot(x,aes(x=gamma_mean,col=chic.cut)) + geom_density() +
geom_vline(aes(xintercept=min,col=chic.cut),data=thr) +
facet_wrap( ~ expt,scales="free_y")

x[,chicago.n35:=chicago.n3>chicago.n5]
ggplot(x,aes(x=log10(chicago.n5),y=log10(gamma_mean))) + geom_jitter(alpha=0.5) + geom_abline() + facet_wrap(chicago.n35~expt)

ggplot(x,aes(x=chicago.n3,y=gamma_mean,col=expt)) + geom_point() + geom_abline() + facet_wrap(~expt)
ggplot(x,aes(x=chicago.n5,y=gamma_mean.c5,col=expt)) + geom_point() + geom_abline() + facet_wrap(~expt)

x[gamma_mean>40 & chicago.n5==1,]

y <- data[baitID==330965 & experiment=="actprom",]
head(y)

y[,lN:=log10(N)]
y2 <- melt(y[,.(dist,lN,mppc,mppi,residual,chicago,z,y)],c("dist","y"))
head(y2)
ggplot(y2,aes(x=dist,y=value)) + geom_path() + facet_grid(variable~.,scales="free_y") + geom_vline(aes(xintercept=dist),data=y2[y==1,],col="red") + xlim(-5e+5,5e+5)
sum(y$mppc)
gamma[baitID==330965,]

x2 <- x[,.(gamma_mean=sum(MPPC)),by=c("expt","baitID")]
d <- "/mrc-bsu/scratch/cew54/peaky/summary/tables" #"/scratch/wallace/peaky/tables"
files <- c(actprom="DT_aCD4pchic_20e6_corr.rds",
           actval="DT_aCD4val_20e6_corr.rds",
           nonprom="DT_nCD4pchic_20e6_corr.rds",
           nonval="DT_nCD4val_20e6_corr.rds")
NPER <- vector("list",4)
names(NPER) <- names(files)
for(nm in names(files)) {
    message(nm)
    x <- readRDS(file.path(d,files[nm]))
    head(x)
    setnames(x,"rjmcmc_pos_g-4.69897","gpos")
    x2 <- x[,.(gamma_mean=sum(gpos)),by="baitID"]
    x2[,experiment:=nm]
    NPER[[nm]] <- copy(x2)
}
nper <- rbindlist(NPER)
chic <- rbindlist(CORR[o])
nper <- merge(nper,chic,by=c("baitID","experiment"))
head(nper)
nper[,expt:=factor(title.ab[experiment],levels=title.ab[o])]
ggplot(nper,aes(x=chicago.n5,y=gamma.mean)) + geom_point() + facet_wrap(~expt) + geom_abline()
res1 <- with(x[chicago.n5>0, ], tapply(chicago.n5, expt,summary))
res2 <- with(x[chicago.n5>0, ], tapply(gamma.mean, expt,summary))


    dim(x2)
    head(x2)
ggplot(x2,aes(x=gamma_mean)) + geom_histogram()


(load(file.path(d,peaky.files[1])))


DATA <- structure(vector("list",length(title)),names=names(title))
for(nm in names(title)) {
   message(" !!! ",nm) 
   fnm <- file.path(d,paste0(nm,"-joined.RData"))
   (load(fnm))
#   if(grepl("prom",nm))
                                        #       x$z <- lgit(x$z)
 #  if(grepl("val",nm))
 #      x$y <- x$y * x$z >= log2(15)
   x$int <- runif(nrow(x),0.99,1.01)
   x$expt <- nm
   DATA[[nm]] <- x[,.(nm,baitID,preyID,N,residual,chicago,mppc,beta.post,mid.bait,mid.prey,dist,z,maxchic,maxres)]
}

lapply(DATA,head,3)
class(DATA[[1]])
DATA <- rbindlist(DATA)

#' distance window
dw <- DATA[,.(maxd=max(abs(dist)),nfrag=.N),by=c("nm","baitID")]
summary(dw)
ggplot(dw,aes(x=maxd)) + geom_histogram() + facet_wrap(~nm,scales="free_y")
ggplot(dw,aes(x=nfrag)) + geom_histogram() + facet_wrap(~nm,scales="free_y")
dw[,.(summary(maxd)),by="nm"]

summary(DATA[,.(max(abs(dist))),by=c("nm","baitID")]$V1)

#' Number of baits with CHICAGO > 5 for at least 1 prey in window
tmp <- DATA[,.(over5=maxchic>5),by=c("nm","baitID")]
(tmp[,.(over5=sum(over5),under5=sum(!over5),total=.N),by=c("nm")])
