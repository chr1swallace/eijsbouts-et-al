library(data.table)
library(ggplot2)
library(magrittr)
library(cowplot)
library(stargazer)


## d <- "/scratch/wallace/peaky/derived"
d <- "/mrc-bsu/scratch/cew54/peaky/summary/derived"
list.files(d)
(load(file.path(d,"alldata.RData")))
## NB max dist for cao is 1e+6, 95% at 2e+5, need to set obs outside this to NA rather than 0 to avoid distance being a confounder
data[grepl("Prom",expt) & abs(dist)>2e+5, y:=NA]

## tmp <- data[expt=="Promoter, non",]$z
## hist(tmp,breaks=100)
## tt <- table(tmp==0,tmp==1)
## tt[1,1]/sum(tt)

## cuts 2/16m rows, but allows logistic model
data[expt %in% c("Promoter, non","Promoter, act") & z!=0 & z!=1, z:=NA]
data <- data[,.(baitID,preyID,chicago,mppc,z,y,expt2,expt,block)]

source("~/Projects/peaky/common.R")

DATA <- split(data,data$expt)
    ## val:
    ## z - transcription log2(count+1)
    ## y - is a bait in the prom expt
    ## prom:
    ## z - propn of frag in E4-E11 states
    ## y - !is.na(cao.score) for fragment - NB, defined only for dist < 1MB

names(DATA)
################################################################################

## per fragment
fitter <- function(din,Y,n=100000,f=mod) {
    cat(".")
    d <- din[!is.na(din[[Y]]) & !is.nan(din[[Y]]),]
    if(nrow(d)>n)
        d <- d[sample(1:nrow(d),n),]
    d$int <- runif(nrow(d),0.99,1.01)
    list(M0= paste(Y,"~int-1") %>% as.formula %>% f(., d),
         Mp= paste(Y,"~sqrt(mppc)") %>% as.formula %>% f(., d),
         Mc= paste(Y,"~asinh(chicago)") %>% as.formula %>% f(., d),
         Mpc= paste(Y,"~asinh(chicago) + sqrt(mppc)") %>% as.formula %>% f(., d))
         ## Mpci= paste(Y,"~asinh(chicago) * sqrt(mppc)") %>% as.formula %>% f(., d))
}

## logistic models of y
MY <- lapply(DATA, fitter, Y="y",f=lmod)
promidx <- c("Promoter, non", "Promoter, act")
validx <- c("Validation, non","Validation, act")
## linear models of z
MZ <- c(lapply(DATA[validx],fitter,Y="z"),
        lapply(DATA[promidx],fitter,Y="z",f=lmod))
names(MZ)

#rm(DATA); gc()

domat <- function(x) {
    names(x) %<>% sub("\\.d\\.f\\.","",.)
    matrix(x,4,length(x)/4,
           dimnames=list(sub(".*\\.","",names(x))[1:4],
                         sub("\\..*","",names(x))[seq(1,length(x),by=4)]))
}
getbic <- function(M)
    do.call("c",M)%>% sapply(.,BIC)
sumbic <- function(M) {
    bic <- getbic(M) %>% domat(.)
    t(t(bic[-1,]) - bic[1,])
}

tab <-rbind(promy=sumbic(MY[promidx]),
            promz=sumbic(MZ[promidx]),
            valy=sumbic(MY[validx]),
            valz=sumbic(MZ[validx]))
tab
tabf <- format.df(tab,dec=1,dollar=TRUE)
## for(i in c(1,5,9,13)) {
for(i in c(1,4,7,10)) {
    for(j in c(1,2)) {
        wh <- which.min(tab[i:(i+2),j])
        tabf[ i+wh-1, j] %<>%paste0("\\textbf{* }",.)
    }
}
tabf

latex(tabf,
      colheads=c("Non-activated","Activated"),
      dec=1,
      col.just=c("r","r"),
      booktabs=TRUE,
      rowlabel="Model",
      rgroup=c("a: Promoter: match to Cao et al.","b: Promoter: link to active chromatin",
               "c: Validation: overlap baited promoter","d: Validation: expression at linked promoter"),
      rowname=rep(c("MPPC","CHiCAGO","MPPC + CHiCAGO"),3),
      caption="$\\Delta$BIC from the intercept only model for four measures of biological plausibility of contacts.  The best fitting model (lowest $\\Delta$BIC) is highlighted by \\textbf{*}.
\\textbf{a}--\\textbf{d} are defined in full in the Methods.  Briefly,
\\textbf{a} whether the bait-prey pair corresponds to published CD4$^+$ T cell promoter-enhancer networks\\cite{Cao2017-qx};
\\textbf{b} whether the prey fragment overlaps active chromatin states defined by \\cite{burren2017chromosome};
\\textbf{c} whether the prey from overlaps a gene promoter;
\\textbf{d} the level of expression of a gene associated with the prey fragment.
In all cases, a robust clustered model was used to account for repeated observations at the prey fragment.",
      label="tab:dbic-alldata", 
      file=file.path(d,"../ptables/dbic-alldata.tex"))

## plot fits
plotter <- function(M) {
    df <- expand.grid(chicago=c(2,3,4,5,10,50),mppc=seq(0,1,by=0.01))
    df <- cbind(df, M=(predict(M,newdata=df,conf.int=0.95)))
    if("lrm" %in% class(M)) {
        df[,-c(1,2)] <- exp(df[,-c(1,2)])
        yl <- 1.2
        yt <- "Predicted odds of overlap"
    } else {
        yl <- 0.5
        yt <- "Predicted expression"
    }
    ggplot(df,aes(x=mppc,y=M.linear.predictors,#ymin=M.lower,ymax=M.upper,
                  fill=as.factor(chicago))) + geom_path(aes(col=as.factor(chicago))) +
    #geom_ribbon(alpha=0.1) +
    #geom_path(aes(y=one.linear.predictors),col="black") +
    #geom_ribbon(aes(ymin=one.lower,ymax=pmin(yl,one.upper)),fill="black",alpha=0.1) + 
    scale_colour_discrete("CHiCAGO score") +
    scale_fill_discrete("CHiCAGO score") +
    scale_x_continuous("MPPC") + scale_y_continuous(yt) +
      guides(col = guide_legend(nrow = 1)) +
      theme( legend.position = #"bottom",
            c(.95, .05),
       legend.justification = c("right", "bottom"))
}
plotter(MY[[1]][[4]])

titles <- outer(c(#"a: published T cell promoter-enhancer network",
"promoter-enhancer network",
            "active chromatin states",
            "gene promoter",
            "gene expression"),
            c("(non)", "(act)"),
            "paste") %>% t() %>% as.vector()
titles

plots <- lapply(c(MY[promidx],MZ[promidx],MY[validx],MZ[validx]), function(M) plotter(M[[4]]))
for(i in seq_along(plots)) {
    plots[[i]] <- plots[[i]] + ggtitle(titles[[i]]) +theme(plot.title = element_text(size = rel(1), hjust = 0,
                                                                   face="bold"))
    if(i < length(plots))
        plots[[i]] <- plots[[i]] + theme(legend.position="none")
}
plot_grid(plotlist=plots,hjust=0,nrow=4,align="hv")#,labels=letters[1:8])

ggsave(file.path(d,"../figures/fits-alldata.pdf"),height=10,width=8)

rm(DATA,MY,MZ); gc()

################################################################################

##' now do same within runs

head(data)
data[,bl:=.N,by=c("expt","baitID","block")]
data <- data[block!=0,]
## data[,h:=NULL]
## data <- merge(data,cs,by=c("expt","baitID","preyID"))
## data[,sc:=sqrt(mppc)/sum(sqrt(mppc)),by=c("expt","baitID","block")]

de <- data[grepl("Promoter",expt) & bl>=2,]
summary(de$bl)
quantile(de$bl,seq(0,1,0.1))
de[,blq:=cut(bl,c(2,4,10,20,200),include.lowest = TRUE)]
DATA <- split(de,paste(de$expt,de$blq))
names(DATA)

fitter <- function(din, X, Y, n=100000, f=mod) {
    d <- din[!is.na(din[[Y]]) & !is.nan(din[[Y]]),]
    if(nrow(d)>n)
        d <- d[sample(1:nrow(d),n),]
    d$int <- runif(nrow(d),0.99,1.01)
    list(M0= paste(Y,"~int-1") %>% as.formula %>% f(., d),
         Mp= paste(Y,"~sqrt(mppc)") %>% as.formula %>% f(., d),
         Mc= paste(Y,"~asinh(chicago)") %>% as.formula %>% f(., d),
         Mpc= paste(Y,"~asinh(chicago) + sqrt(mppc)") %>% as.formula %>% f(., d))
         ## Mpci= paste(Y,"~asinh(chicago) * sqrt(mppc)") %>% as.formula %>% f(., d))
}


## logistic models of y
MY <- lapply(DATA, fitter, Y="y",f=lmod)
MZ <- lapply(DATA, fitter, Y="z",f=lmod)

MM <- MZ
cf <- sapply(MM, function(M) coefficients(M[[4]]))
se <- sapply(MM, function(M) vcov(M[[4]])  %>%  diag())
df <- merge(melt(cf),melt(se),by=c("Var1","Var2"),suffixes=c(".cf",".se"))
df$x <- sub(".*, ","",df$Var2)
df$cell <- sub(" .*","",df$x)
df$x <- sub(".* ","",df$x)
    
ggplot(df[df$Var1!="Intercept",], aes(x=x,y=value.cf,ymin=value.cf-1.96*value.se,ymax=value.cf+1.96*value.se)) + geom_pointrange() + facet_grid(cell~Var1) + geom_hline(yintercept = 0,col="red")


## promidx <- c("Promoter, non", "Promoter, act")
## validx <- c("Validation, non","Validation, act")
## ## linear models of z
## MZ <- c(lapply(DATA[validx],fitter,Y="z"),
##         lapply(DATA[promidx],fitter,Y="z",f=lmod))
## names(MZ)

#rm(DATA); gc()

domat <- function(x) {
    names(x) %<>% sub("\\.d\\.f\\.","",.)
    matrix(x,4,length(x)/4,
           dimnames=list(sub(".*\\.","",names(x))[1:4],
                         sub("\\..*","",names(x))[seq(1,length(x),by=4)]))
}
getbic <- function(M)
    do.call("c",M)%>% sapply(.,BIC)
sumbic <- function(M) {
    bic <- getbic(M) %>% domat(.)
    t(t(bic[-1,]) - bic[1,])
}

## tab <-rbind(promy=sumbic(MY[promidx]),
##             promz=sumbic(MZ[promidx]),
##             valy=sumbic(MY[validx]),
##             valz=sumbic(MZ[validx]))
tab <-rbind(promy=sumbic(MY),
            promz=sumbic(MZ))
tab
tabf <- format.df(tab,dec=1,dollar=TRUE)
## for(i in c(1,5,9,13)) {
for(i in c(1,4)) {
    for(j in 1:ncol(tab)) {
        wh <- which.min(tab[i:(i+2),j])
        tabf[ i+wh-1, j] %<>%paste0("\\textbf{* }",.)
    }
}
colnames(tabf)
tabf <- rbind(tabf[,c(8,7,5,6)],tabf[,c(4,3,1,2)])
tabf
latex(tabf,
      colheads=c("2-4","5-10","11-20","21+"),
      dec=1,
      col.just=c("r","r","r","r"),
      booktabs=TRUE,
      rowlabel="Model (cell)",
      collabel="Stretch",
      rgroup=c("a: match to Cao et al. (non)",
               "b: link to active chromatin (non)",
               "c: match to Cao et al. (act)",
               "d: link to active chromatin (act)"),
      rowname=rep(c("MPPC","CHiCAGO","MPPC + CHiCAGO"),3),
      caption="$\\Delta$BIC from the intercept only model for four measures of biological plausibility of contacts.  The best fitting model (lowest $\\Delta$BIC) is highlighted by \\textbf{*}.
\\textbf{a}--\\textbf{d} are defined in full in the Methods.  Briefly,
\\textbf{a} whether the bait-prey pair corresponds to published CD4$^+$ T cell promoter-enhancer networks\\cite{Cao2017-qx};
\\textbf{b} whether the prey fragment overlaps active chromatin states defined by \\cite{burren2017chromosome};
\\textbf{c} whether the prey from overlaps a gene promoter;
\\textbf{d} the level of expression of a gene associated with the prey fragment.
In all cases, a robust clustered model was used to account for repeated observations at the prey fragment.",
      label="tab:dbic-alldata", 
      file=file.path(d,"../ptables/dbic-alldata.tex"))

## plot fits
plotter <- function(M) {
    df <- expand.grid(chicago=c(2,3,4,5,10,50),mppc=seq(0,1,by=0.01))
    df <- cbind(df, M=(predict(M,newdata=df,conf.int=0.95)))
    if("lrm" %in% class(M)) {
        df[,-c(1,2)] <- exp(df[,-c(1,2)])
        yl <- 1.2
        yt <- "Predicted odds of overlap"
    } else {
        yl <- 0.5
        yt <- "Predicted expression"
    }
    ggplot(df,aes(x=mppc,y=M.linear.predictors,#ymin=M.lower,ymax=M.upper,
                  fill=as.factor(chicago))) + geom_path(aes(col=as.factor(chicago))) +
    #geom_ribbon(alpha=0.1) +
    #geom_path(aes(y=one.linear.predictors),col="black") +
    #geom_ribbon(aes(ymin=one.lower,ymax=pmin(yl,one.upper)),fill="black",alpha=0.1) + 
    scale_colour_discrete("CHiCAGO score") +
    scale_fill_discrete("CHiCAGO score") +
    scale_x_continuous("MPPC") + scale_y_continuous(yt) +
      guides(col = guide_legend(nrow = 1)) +
      theme( legend.position = #"bottom",
            c(.95, .05),
       legend.justification = c("right", "bottom"))
}
plotter(MY[[1]][[4]])

titles <- outer(c(#"a: published T cell promoter-enhancer network",
"promoter-enhancer network",
            "active chromatin states",
            "gene promoter",
            "gene expression"),
            c("(non)", "(act)"),
            "paste") %>% t() %>% as.vector()
titles

plots <- lapply(c(MY[promidx],MZ[promidx],MY[validx],MZ[validx]), function(M) plotter(M[[4]]))
plots <- lapply(MY, function(M) plotter(M[[4]]))
for(i in seq_along(plots)) {
    plots[[i]] <- plots[[i]] + #ggtitle(titles[[i]]) +
        theme(plot.title = element_text(size = rel(1), hjust = 0,
                                                                   face="bold"))
    if(i < length(plots))
        plots[[i]] <- plots[[i]] + theme(legend.position="none")
}
plot_grid(plotlist=plots,hjust=0,nrow=4,align="hv")#,labels=letters[1:8])

ggsave(file.path(d,"../figures/fits-alldata.pdf"),height=10,width=8)

rm(DATA,MY,MZ); gc()


quantile(DATA[[2]]$bl)
args <- expand.grid(j=c(1,2),
                    Y=c("y","z"),
                    X=c("sqrt(mppc)","h","sc"),
                    blmin=c(2,6,11,21,51))
args$blmax <- c("2"=3,"4"=5,"6"=10,"11"=20,"21"=50,"51"=400)[as.character(args$blmin)]
args
M <- lapply(1:nrow(args),function(i) 
    fitter(DATA[[args[i,"j"]]], Y=args[i,"Y"], X=args[i,"X"], blmin=args[i,"blmin"], blmax=args[i,"blmax"])
    )

M <- do.call("rbind",M) %>% as.data.frame()
M <- cbind(M,args)
M
M$cell <- c("non-activated","activated")[M$j]
M <- subset(M,!(M$Y=="y" & M$blmin==51))
head(M)
M$yz <- ifelse(M$Y=="y", "a: Match to Cao et al.", "b: Overlap open chromatin")
M$bl <- paste0(M$blmin,"-",M$blmax) %>% sub("51-400","51+",.) %>% factor(.,levels=c("2-3","4-5","6-10","11-20","21-50","51+"))

M[M$X=="sqrt(mppc)","bic"]
M[M$X=="sc","bic"]
M[M$X=="h","bic"]
M[M$X=="sqrt(mppc)","bic"] - M[M$X=="h","bic"]
M[M$X=="sqrt(mppc)","bic"] - M[M$X=="sc","bic"]

ggplot(M[M$X=="sqrt(mppc)",],aes(x=bl,y=cf.mppc,ymin=l.mppc,ymax=u.mppc)) + geom_pointrange() + facet_grid(Y~cell) +
  geom_hline(yintercept=0) +
  labs(x="Run length", y="log odds ratio + 95% CI")
ggsave(file.path(d,"../figures/mppc-enrichment-within-blocks.pdf"),height=6,width=6)



fitter <- function(din,Y,n=100000,f=mod,blmin,blmax) {
    d <- din[bl>=blmin & bl<=blmax & !is.na(din[[Y]]) & !is.nan(din[[Y]]),]
    if(nrow(d)>n)
        d <- d[sample(1:nrow(d),n),]
    lapply(v, function(thr) {
        cat(".")
        d[,hthr:=h>=thr]
        tt <- table(d$hthr)
        if(length(tt)!=2 | any(tt<100))
            return(NULL)
        names(tt) <- c("h0","h1")
        t2 <- table(d$hthr,d[[Y]])
        if(length(t2)!=4)
            return(NULL)
        y1 <- t2[2,2]/sum(t2[2,])
        y0 <- t2[1,2]/sum(t2[1,])
        m <- paste(Y,"~hthr") %>% as.formula %>% f(., d)
        cf <- coefficients(m)[2]
        s <- sqrt(vcov(m)[2,2])
        c(thr=thr,blmin=blmin,blmax=blmax,tt,prop.y1.retained=y1,prop.y1.dropped=y0,cf=cf,s=s,l=cf-s,u=cf+s,p=2*pnorm(abs(cf/s),lower.tail=FALSE))
    })
}

library(parallel)
options(mc.cores=2)
DATA <- split(data,data$expt)
#MY <- mclapply(DATA[promidx], function(d) fitter(d,Y="y",f=lmod) %>% do.call("rbind",.)) %>%
#  lapply(., as.data.frame)
MZ1 <- mclapply(DATA[promidx], function(d) fitter(d,Y="z",f=lmod,blmin=6,blmax=10) %>% do.call("rbind",.)) %>%
  lapply(., as.data.frame)
MZ2 <- mclapply(DATA[promidx], function(d) fitter(d,Y="z",f=lmod,blmin=11,blmax=50) %>% do.call("rbind",.)) %>%
  lapply(., as.data.frame)
MZ3 <- mclapply(DATA[promidx], function(d) fitter(d,Y="z",f=lmod,blmin=51,blmax=400) %>% do.call("rbind",.)) %>%
  lapply(., as.data.frame)
MZ4 <- mclapply(DATA[promidx], function(d) fitter(d,Y="z",f=lmod,blmin=2,blmax=5) %>% do.call("rbind",.)) %>%
  lapply(., as.data.frame)

M <- list(MZ1,MZ2,MZ3,MZ4)
for(i in seq_along(M)) {
    M[[i]][[1]]$cell <- "non-activated"
    M[[i]][[2]]$cell <- "activated"
}
M <- do.call("c",M) %>% do.call("rbind",.)

## MY[[2]]$cell <- "act"
## MZ[[1]]$cell <- NA
## MZ[[2]]$cell <- NA
## MY[[1]]$ti <- "a: Match to Cao et al (non)"
## MY[[2]]$ti <- "a: Match to Cao et al (act)"
## MZ[[1]]$ti <- "b: Overlap active chromatin (non)"
## MZ[[2]]$ti <- "b: Overlap active chromatin (act)"
## M <- c(MY,MZ) %>% do.call("rbind",.)
## M
head(M)
M$bl <- paste0(M$blmin,"-",M$blmax) %>% sub("51-400","51+",.) %>% factor(.,levels=c("2-5","6-10","11-50","51+"))

p3 <- ggplot(M,aes(x=1-thr,col=bl)) +
  geom_point(aes(y=prop.y1.retained),pch="+") +
  geom_point(aes(y=prop.y1.dropped),pch="o") +
  labs(y="Propn overlap",x="Proportion of posterior covered") +
  ## geom_hline(yintercept=1,col="grey") +
  facet_wrap(~cell)+ theme(legend.position="none",
        strip.background=element_rect(fill="white"),
        strip.text=element_text(face="bold"))
p2 <- ggplot(M,aes(x=1-thr,y=exp(cf.hthr),ymin=exp(l.hthr),ymax=exp(u.hthr),col=bl)) + geom_pointrange() +
  labs(y="Odds ratio",x="Proportion of posterior covered") +
  geom_hline(yintercept=1,col="grey") +
  facet_wrap(~cell)+ theme(legend.position="none",
        strip.background=element_rect(fill="white"),
        strip.text=element_text(face="bold"))
p1 <- ggplot(M[!is.na(M$cell),],aes(x=1-thr,y=100*h1/(h1+h0),col=bl)) + geom_point() + geom_path() + 
  labs(y="% of fragments kept",x="Proportion of posterior covered") +
  scale_colour_discrete("Run length") +
  theme(legend.position=c(0.05,0.95),legend.justification=c("left","top"),
        strip.background=element_rect(fill="white"),
        strip.text=element_text(face="bold")) +
  facet_wrap(~cell) 
plot_grid(p1,p2,p3,ncol=1)

ggsave(file.path(d,"../figures/credset-enrichment.pdf"),height=8,width=8)

## val:
    ## z - transcription log2(count+1)
    ## y - is a bait in the prom expt
    ## prom:
    ## z - propn of frag in E4-E11 states
    ## y - !is.na(cao.score) for fragment - NB, defined only for dist < 1MB

do.call("rbind",M)

lapply(M, coefficients)
lapply(M,vcov)

ggplot(y1,aes(x=sqrt(mppc),y=z)) + facet_wrap(~expt,scales="free") + geom_point(alpha=0.1) + geom_smooth()  

DATA <- list(C0=y0[expt %in% c("Promoter, non","Promoter, act") & !is.na(z),] ,
             C1=y1[expt %in% c("Promoter, non","Promoter, act") & !is.na(z),] ,
             C5=y[expt %in% c("Promoter, non","Promoter, act") & !is.na(z) & block_length==5,] ,
             C10=y[expt %in% c("Promoter, non","Promoter, act") & !is.na(z) & block_length>5 & block_length<=10,] ,
             C20=y[expt %in% c("Promoter, non","Promoter, act") & !is.na(z) & block_length>10 & block_length<=20,])
lapply(DATA, function(x) table(x$z))

DATA <- lapply(DATA, function(x) {
    split(x, x$expt)[1:2]
}) %>% do.call("c",.)

## use max 1 million obs
maxn <- 1e+6
n <- sapply(DATA,nrow)
for(i in which(n>maxn))
    DATA[[i]] <- DATA[[i]][ sample(1:n[i],maxn), ]

## model
MODS <- lapply(DATA, function(x) lmod(z ~ sqrt(mppc), x))
lapply(MODS,coefficients)

stargazer(MODS)


MODS <- lapply(DATA, function(d) {
    split(d,d$expt)[1:2] %>% lapply(.,function(x) lmod(z ~ sqrt(mppc),x))
})
MODS

y1p <-
M1 <- split(y1p,y1p$expt)[1:2] %>% lapply(.,function(x) lmod(z ~ sqrt(mppc),x))
M1

y5p <-
M1 <- split(y1p,y1p$expt)[1:2] %>% lapply(.,function(x) lmod(z ~ sqrt(mppc),x))
M1




MZ
sumbic(MZ)

sd <- data.table(p,bl,gapbefore,gapafter,beg,end,drop)
    sd[,drop:=any(drop),by="bl"]
    sd$drop
    
y[,gapbefore:=c(diff(preyID),0),by=c("expt","baitID")]
y[,gapafter:=c(0,diff(preyID)),by=c("expt","baitID")]
y
y[,beg:=1:.N==1,by=c("expt","baitID")]
y[,end:=1:.N==.N,by=c("expt","baitID")]
dropblocks <- y[beg==TRUE & gapbeforeends==TRUE & (abs(gapbefore)<3 | abs(gapafter)<3),]
dropblocks$drop <- TRUE
y <- merge(y,

head(y)
y <- y[order(expt,baitID,block,1-mppc),]
head(y[block_length>1,])
y[baitID==254 & block==5,]
y[,mppc.cov:=cumsum(mppc)/sum(mppc),by=c("expt","baitID","block")]
y[,z.cov:=cumsum(z)/sum(z,na.rm=TRUE),by=c("expt","baitID","block")]
y[,y.cov:=cumsum(y)/sum(y,na.rm=TRUE),by=c("expt","baitID","block")]
y[,mppc.el:=1:.N,by=c("expt","baitID","block")]
y[,group:=factor(paste(expt,baitID,block))]
thr <- 0.5
y[is.na(y) | is.nan(y), y:=0]
y[is.na(z) | is.nan(z), z:=0]
f <- function(thr) {
    y[is.nan
    yt <- unique(y[mppc.cov>=thr & block_length>2,],
                 by=c("expt","baitID","block"))
    yt[,bl:=cut(block_length,c(2,10,50,100,150))]
    yt[,thr:=thr]
    yt[,.(expt,baitID,bl,block,thr,mppc.cov,y.cov,z.cov,mppc.prop=mppc.el/block_length)]
}
propcov <- lapply(c(seq(0.1,0.9,0.1)),f) %>% rbindlist()

propcov

ggplot(propcov,aes(x=mppc.cov,y=z.cov,col=expt))+
geom_hex() + geom_smooth(col="black") + 
facet_grid(bl~expt) +
labs(x="MPPC cov",y="y cov") +
theme(legend.position="bottom") + background_grid()

    ggplot(propcov,aes(x=mppc.cov,y=y.cov,col=expt))+
geom_point() + geom_smooth() + 
facet_grid(bl~expt) +
labs(x="MPPC cov",y="y cov") +
theme(legend.position="bottom") + background_grid()






nm <- "actval"

DATA <- MODELSPLUS <- CMODELSPLUS <- MODELS <- CMODELS <- PLOTS <- structure(vector("list",length(title)),names=names(title))
for(nm in names(title)) {
   message(" !!! ",nm) 
   fnm <- file.path(d,paste0(nm,"-joined-fm.RData"))
   (load(fnm))
   if(grepl("prom",nm)) {
       x[abs(dist) > 1e+6, y:=NA]
       ##     x <- merge(x,cao,by=c("baitID","preyID"),all.x=TRUE)
       ##     x[,y:=ifelse(is.na(cao),0,1)]
   }
   ##       x$z <- lgit(x$z)
                                        #  if(grepl("val",nm))
                                        #      x$y <- x$y * x$z >= log2(15)
   x$int <- runif(nrow(x),0.99,1.01)
   DATA[[nm]] <- x[!is.na(chicago) & chicago > 2,]
}

DATA5 <- lapply(DATA, function(x) x[chicago > 5,])

fitter <- function(DATA) {
    ## val:
    ## z - transcription log2(count+1)
    ## y - is a bait in the prom expt
    ## prom:
    ## z - propn of frag in E4-E11 states
    ## y - !is.na(cao.score) for fragment - NB, defined only for dist < 1MB

    ## model quantitative trait, z
    DNOMISS <- lapply(DATA, function(x) { x <- x[!is.na(z) & !is.nan(z),]
        x[,dist:=log(abs(dist))]
        x[,x1:=sqrt(mppc)]
        x[,x2:=asinh(chicago)]
        return(x) })
    MODELS0 <- lapply(DNOMISS, function(x) mod(z ~ int-1,x))
    MODELS <- lapply(DNOMISS, function(x) mod(z ~ sqrt(mppc),x))
    MODELSPLUS <- lapply(DNOMISS, function(x) mod(z ~ sqrt(mppc) + asinh(chicago) ,x))
    MODELSPLUS2 <- lapply(DNOMISS, function(x) mod(z ~ asinh(chicago) ,x))
    
    ## model binary trait, y
    DNOMISS <- lapply(DATA, function(x) { x <- x[!is.na(y),]; return(x) })
    CMODELS0 <- lapply(DNOMISS, function(x) lmod(y ~ int-1,x))
    CMODELS <- lapply(DNOMISS, function(x) lmod(y ~ sqrt(mppc),x))
    CMODELSPLUS <- lapply(DNOMISS, function(x) lmod(y ~ sqrt(mppc) + asinh(chicago),x))
    CMODELSPLUS2 <- lapply(DNOMISS, function(x) lmod(y ~ asinh(chicago),x))
    
    prom.L1 <- list(MODELSPLUS2[[3]],MODELS[[3]],MODELSPLUS[[3]],
                    MODELSPLUS2[[1]],MODELS[[1]],MODELSPLUS[[1]])
    prom.L0 <- list(MODELS0[[3]],
                    MODELS0[[1]])
    val.L1 <- list(MODELSPLUS2[[4]],MODELS[[4]],MODELSPLUS[[4]],
                   MODELSPLUS2[[2]],MODELS[[2]],MODELSPLUS[[2]])
    val.L0 <- list(MODELS0[[4]],
                   MODELS0[[2]])
    cprom.L1 <- list(CMODELSPLUS2[[3]],CMODELS[[3]],CMODELSPLUS[[3]],
                     CMODELSPLUS2[[1]],CMODELS[[1]],CMODELSPLUS[[1]])
    cprom.L0 <- list(CMODELS0[[3]],
                     CMODELS0[[1]])
    cval.L1 <- list(CMODELSPLUS2[[4]],CMODELS[[4]],CMODELSPLUS[[4]],
                    CMODELSPLUS2[[2]],CMODELS[[2]],CMODELSPLUS[[2]])
    cval.L0 <- list(CMODELS0[[4]],
                    CMODELS0[[2]])
    return(list(prom.L1=prom.L1,prom.L0=prom.L0,
                val.L1=val.L1,val.L0=val.L0,
                cprom.L1=cprom.L1,cprom.L0=cprom.L0,
                cval.L1=cval.L1,cval.L0=cval.L0))
}


mylatex <- function(L,L0,file,...) {
    bics <- sapply(L,BIC)
    b0 <- sapply(L0,BIC)
    n <- length(bics)/2
    bics[1:n] <- bics[1:n] - b0[1]
    bics[-c(1:n)] <- bics[-c(1:n)] - b0[2]
                                        #bics <- bics[-1] - bics[1]
    x <- stargazer(L,
            add.lines=list(BIC=c("$\\Delta$BIC",sprintf("%16.1f",bics))),
            keep.stat=c("rsq","n"),
            dep.var.caption="",
            dep.var.labels.include=FALSE,
            report="vcs",
            omit="Constant",
            omit.table.layout="n",
            column.labels=c("non-activated","activated"),
            column.separate=c(n,n),
            #ci=rep(TRUE,length(L)),
            ...
            )
    x %<>% sub("[-1.8ex]","",.,fixed=TRUE)
    hl <- grep("hline",x)
    ## first hline
    if(length(hl)) {
    x[ hl[[1]] ] <- "\\toprule"
    x[ hl[[1]]+1 ] <- ""
    x[ hl[3:(length(hl)-2)] ] <- "\\midrule"
    x[ hl[ length(hl)-1] ] <- ""
      x[ hl[length(hl)]] <- "\\bottomrule"
    }
        cat(x,file=file,sep="\n")
}
xtra <- "against CHiCAGO score and/or MPPC for non-activated (1-3) and activated (4-6) CD4$^+$ T cells for promoter-capture sets.  MPPC and CHiCAGO scores were sqrt and asinh transformed, respectively. $\\Delta$BIC is the difference in Bayesian Information Criterion vs the null model (intercept only)."


fit5 <- fitter(DATA5)

mylatex(fit5$prom.L1,fit5$prom.L0,        
        title=paste("Higher MPPC is independently predictive of greater fragment overlap with active chromatin states amongst fragments with CHiCAGO scores $>$5.  Entries show regression coefficients (std. error) from regressions of the proportion of fragments in active chromatin states ",xtra),
        label="tab:E4E11reg",
        file=file.path(d,"../ptables/tab-E4E11reg.tex"))

mylatex(fit5$cprom.L1,fit5$cprom.L0,        
        title=paste("Higher MPPC is independently predictive that a bait-prey pair correspond to an independent set of predicted enhancer-promoters in total CD4$^+$ T cells, amongst bait-prey pairs with CHiCAGO scores $>$5.  Entries show regression coefficients (std. error) from logistic regressions of pair overlap",xtra),
        label="tab:E4E11reg",
        file=file.path(d,"../ptables/tab-E4E11reg.tex"))

## mylatex(L1,L0,
##         title=paste("Higher MPPC is independently predictive of higher RNA transcription. Entries show regression coefficients (std. error) from regressions of log2 read counts in genes with promotersoverlapping fragments",xtra),
##         label="tab:RNAseqreg",
##         file=file.path(d,"../ptables/tab-RNAseqreg.tex"))

mylatex(fit5$cval.L1,fit5$cval.L0,
        title=paste("Higher MPPC is independently predictive that a fragment contains an annotated promoter amongst fragments with CHiCAGO scores $>$5.  Entries show regression coefficients (std. error) from logistic regressions of fragment overlap with annotated TSS",xtra),
        label="tab:TSSreg",
        file=file.path(d,"../ptables/tab-TSSreg.tex"))

## what about chic scores < 5

fitter <- function(DATA,vars=c("both")) {
    DNOMISS <- lapply(DATA, function(x) {
        x <- x[!is.na(z) & !is.nan(z),]
        x[,dist:=log(abs(dist))]
        x[,x1:=sqrt(mppc)]
        x[,x2:=asinh(chicago)]
        return(x) })
    if(vars=="both")
        lapply(DNOMISS, function(x) mod(z ~ sqrt(mppc) * asinh(chicago),x))
    else
        lapply(DNOMISS, function(x) mod(z ~ sqrt(mppc) , x))
}
yfitter <- function(DATA,vars=c("both")) {
    DNOMISS <- lapply(DATA, function(x) {
        x <- x[!is.na(y) & !is.nan(y),]
        x[,dist:=log(abs(dist))]
        x[,x1:=sqrt(mppc)]
        x[,x2:=asinh(chicago)]
        return(x) })
    if(vars=="both")
        lapply(DNOMISS, function(x) lmod(y ~ sqrt(mppc) * asinh(chicago),x))
    else
        lapply(DNOMISS, function(x) lmod(y ~ sqrt(mppc) , x))
}

fitboth <- fitone <- vector("list",length(DATA))
idx <- grep("Promoter",title)
fitboth[idx] <- fitter(DATA[idx],vars="both")
fitone[idx] <- fitter(DATA[idx],vars="one")
fitboth[-idx] <- yfitter(DATA[-idx],vars="both")
fitone[-idx] <- yfitter(DATA[-idx],vars="one")

plotter <- function(mboth,mone) {
    df <- expand.grid(chicago=c(2,3,4,5,10,50),mppc=seq(0,1,by=0.01))
    df <- cbind(df, both=(predict(mboth,newdata=df,conf.int=0.95)))
    df <- cbind(df, one=(predict(mone,newdata=df,conf.int=0.95)))
    if("lrm" %in% class(mone)) {
        df[,-c(1,2)] <- exp(df[,-c(1,2)])
        yl <- 1.2
        yt <- "Predicted odds of overlap with promoter bait"
    } else {
        yl <- 0.5
        yt <- "Predicted overlap with E4-E11"
    }
    df$zone <- exp(predict(mone,newdata=df))
    ggplot(df,aes(x=mppc,y=both.linear.predictors,ymin=both.lower,ymax=pmin(yl,both.upper),
                  fill=as.factor(chicago))) + geom_path(aes(col=as.factor(chicago))) +
    #geom_ribbon(alpha=0.1) +
    #geom_path(aes(y=one.linear.predictors),col="black") +
    #geom_ribbon(aes(ymin=one.lower,ymax=pmin(yl,one.upper)),fill="black",alpha=0.1) + 
    scale_colour_discrete("CHiCAGO score") +
    scale_fill_discrete("CHiCAGO score") +
    scale_x_continuous("MPPC") + scale_y_continuous(yt,limits=c(0,yl)) #+
    #theme(legend.position="bottom")
}

plots <- mapply(plotter,fitboth,fitone,SIMPLIFY=FALSE)
for(i in seq_along(title)) {
    plots[[i]] <- plots[[i]] + ggtitle(title[[i]])
    if(i==3)
        plots[[i]] <- plots[[i]] + theme(legend.position=c(0.25,0.1),legend.direction = "horizontal",legend.background = element_rect(color="black",size = 1, linetype = "solid"))
    else
        plots[[i]] <- plots[[i]] + theme(legend.position="none")
}
plot_grid(plotlist=plots[c(1,3,2,4)],labels="auto")
ggsave(file.path(d,"../figures","MPPC-vs-outcome.pdf"),height=10,width=10)


## fit4 <- fitter(DATA4)
## fit3 <- fitter(DATA3)
## fit2 <- fitter(DATA2)

## proms <- list(non.5=fit5$prom.L1[[3]],non.4=fit4$prom.L1[[3]],non.3=fit3$prom.L1[[3]],non.2=fit2$prom.L1[[3]],
##               act.5=fit5$prom.L1[[6]],act.4=fit4$prom.L1[[6]],act.3=fit3$prom.L1[[6]],act.2=fit2$prom.L1[[6]]) 
## lapply(proms,coefficients)
## vals <- list(non.5=fit5$val.L1[[3]],non.4=fit4$val.L1[[3]],non.3=fit3$val.L1[[3]],non.2=fit2$val.L1[[3]],
##               act.5=fit5$val.L1[[6]],act.4=fit4$val.L1[[6]],act.3=fit3$val.L1[[6]],act.2=fit2$val.L1[[6]]) 
## lapply(vals,coefficients)

## makeplotdata <- function(m,nm) {
##     cf <- coef(m)[-1]
##     v <- diag(vcov(m))[-1]
##     lo <- cf - 1.96 * sqrt(v)
##     hi <- cf + 1.96 * sqrt(v)
##     data.frame(var=names(cf),cf=cf,v=v,lo=lo,hi=hi,nm=nm)
## }
## library(ggplot2)
## library(magrittr)
## plotter <- function(L) {
## df <- mapply(makeplotdata,L,names(L),SIMPLIFY=FALSE)%>%do.call("rbind",.)
## df$chicago <- factor(gsub("non.|act.","",df$nm),levels=2:5)
## levels(df$chicago) <- c("(2,3]","(3,4]","(4,5]",">5")
## df$cell <- substr(df$nm,1,3)
## w=0.1
## df$x <- as.numeric(df$chicago) + ifelse(df$cell=="non",-1,1)*w

## ggplot(df,aes(x=x,y=cf,ymin=lo,ymax=hi,col=cell)) + geom_pointrange() + facet_wrap(~var) +
## scale_x_continuous("CHiCAGO score",breaks=seq_along(levels(df$chicago)),labels=levels(df$chicago))
## }

## plotter(proms)
## plotter(vals)

## library(texreg)
## plotreg(proms[[1]])
## plotreg(proms)






##    ## plots
## ##    if(nrow(tmp)>100000) {
## ##        tmp <- tmp[sample(1:nrow(tmp),100000),]
## ##    }

## ##    xti <- if(grepl("prom",nm)) { "Fraction covering E4-E11" } else { "Mean RNA-seq readcount (log2)" }
## ##    mlabels <- c(0,0.01,0.05,0.1,0.2,0.5,0.8)
## ##    mbreaks <- sqrt(mlabels)

## ##    PLOTS[[nm]] <- ggplot(tmp,aes(y=z,x=sqrt(mppc))) + geom_point() + geom_smooth(method="lm") +
## ##    scale_x_continuous("MPPC (sqrt scale)",breaks=mbreaks,labels=mlabels) +
## ##    scale_y_continuous(xti) + background_grid() #+ ggtitle(title[nm])
## ## }
## ## lapply(MODELS,BIC)

