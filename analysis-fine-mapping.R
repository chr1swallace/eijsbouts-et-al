library(data.table)
library(ggplot2)
library(magrittr)
library(cowplot)
library(stargazer)


## d <- "/scratch/wallace/peaky/derived"
d <- "/mrc-bsu/scratch/cew54/peaky/summary/derived"
list.files(d)
peaky.files <- list.files(d,pattern="peaky.RData")
chic.files <- list.files(d,pattern="chicago")
(load(file.path(d,"cao.RData")))
b=17508
b=17777
b=666461 # AHR
source("~/Projects/peaky/common.R")

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

