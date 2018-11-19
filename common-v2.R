d <- "/mrc-bsu/scratch/cew54/peaky/summary/derived"
library(data.table)
addpositions <- function(x) {
    (load(file.path(d,"hind-position.RData")))
  hind[,hindID:=as.integer(hindID)]
  hind[,mid:=(start+end)/2]
    x <- merge(x,hind[,.(hindID,mid)],by.x="baitID",by.y="hindID")
    x <- merge(x,hind[,.(hindID,mid)],by.x="preyID",by.y="hindID",suffixes = c(".bait",".prey"))
    x[,dist :=  mid.prey - mid.bait]
    return(x)
}

## robust cluster regression models
library(rms)
lmod <- function(f,x,maxn=NULL) {
    f <- as.formula(f)
    m <- lrm(f, x=TRUE, y=TRUE, data=x)
    robcov(m, cluster=x$preyID)
}

mod <- function(f, x, do.print=TRUE,dropfirst=TRUE) {
#    print(f)
    if(!("formula" %in% class(f)))
        f <- as.formula(f)
    mod <- ols(f, x=TRUE, y=TRUE, data=x)
    mod.r <- robcov(mod, cluster=x$preyID)
    return(mod.r)
}

modprint <- function(mod.r,dropfirst=TRUE) {
    cf <- coefficients(mod.r)
    se <- ( vcov(mod.r) %>% diag() %>% sqrt() )
    z <- cf/se
    p <- pnorm(abs(z),lower.tail=FALSE) * 2
    ret <- cbind(Coef=cf, Lower.95=cf-1.96*se, Upper.95=cf+1.96*se, P=p)
    if(dropfirst)
        return(ret[-1,,drop=FALSE])
    return(ret)
}

lgit <- function(x) ifelse(x==0,NA,log(x/(1-x)))

## plots
library(ggplot2)
library(cowplot)
cur <- theme_get()
theme_set(cur %+replace%
          theme(plot.title = element_text(size = rel(1.2), hjust = 0, vjust = 1, )#,
                ##                 panel.grid.major = element_line(colour = "grey90", size = 0.2))
                ))

title <- c("aCD4pchic"="Promoter capture, activated",
           aCD4val="Validation, activated",
           nCD4pchic="Promoter capture, non-activated",
           nCD4val="Validation, non-activated")
title.pv <- c(aCD4pchic="Promoter Capture",
              aCD4val="Validation",
              nCD4pchic="Promoter Capture",
              nCD4val="Validation")
title.ab <- c(aCD4pchic="Promoter, act",
              aCD4val="Validation, act",
              nCD4pchic="Promoter, non",
              nCD4val="Validation, non")
o <- c(3,1,4,2)

addtitles <- function(x) {
    x[,expt:=factor(title.ab[experiment],levels=title.ab[o])]
    x[,expt2:=factor(title.pv[experiment])]
    x
}
