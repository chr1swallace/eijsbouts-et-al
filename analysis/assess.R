library(data.table)
library(ggplot2)
library(magrittr)
library(cowplot)
library(stargazer)
##     baitID preyID   N   chicago    mppc mid.bait mid.prey     dist  z y
##  1:    248    236  15 1.6451756 0.01330  1116188  1041676 -74511.5  0 1
##  2:    248    237  14 1.4126643 0.01350  1116188  1045696 -70492.0  0 0
##  3:    248    238  38 5.3406739 0.01965  1116188  1048650 -67538.0  0 0
##  4:    248    239  16 0.6600332 0.02580  1116188  1055856 -60332.0 NA 1
##  5:    248    240  66 8.9079607 0.06170  1116188  1063452 -52735.5  0 0
##  6:    248    241  75 9.7373344 0.07780  1116188  1066646 -49542.5  0 0
##  7:    248    242  81 7.8683215 0.13420  1116188  1076502 -39685.5 NA 0
##  8:    248    243 107 7.4648761 0.17775  1116188  1087597 -28591.0 NA 0
##  9:    248    244 104 5.5184275 0.15920  1116188  1091624 -24564.0 NA 0
## 10:    248    245 151 7.9368161 0.15045  1116188  1094376 -21812.0 NA 0
## 11:    248    246 282 5.2219046 0.14055  1116188  1101928 -14260.5  0 0
source("~/Projects/peaky/common-v2.R")

## d <- "/scratch/wallace/peaky/derived"
## v2data <- "~/bsu/peaky-orig/chrisw"
d <- "~/bsu/peaky"
list.files(d)
list.files(file.path(d,"analysis"))

(data.files <- file.path(d,"analysis",paste0(names(title),".RData")))
data <- lapply(data.files, function(f) {
    tmp <- eval(as.symbol(load(f)))
    tmp[,experiment:=sub(".RData","",basename(f))]
    tmp })
names(data) <- sub(".RData","",basename(data.files))

## data used to evaluate biological plausibility
(load(file.path(d,"annot","cao.RData"))) ## cao et al enh-promoter interactions from total CD4
## NB max dist for cao is 1e+6, 95% at 2e+5, need to set obs outside this to NA rather than 0 to avoid distance being a confounder
(load(file.path(d,"annot","b2g.RData")))
b2g <- b2g[biotype %in% c("protein_coding"),]
(load(file.path(d,"annot","rnaseq.RData")))


## act transcribed
## z = average expression at promoters
## y = prey is a bait in the promcmp expt
rnaseq[,trans:=act_1 + act_2]
a <- merge(b2g[,.(id,gene,baitID)],rnaseq[,.(id,trans)],by="id")
setnames(a,"baitID","preyID")
a <- a[,.(z=mean(log2(trans+1))),by="preyID"]
data[["aCD4val"]]  %<>% merge(., a, by="preyID",all.x=TRUE)
data[["aCD4val"]][,y:=as.integer(!is.na(z))]# z only not NA if corresponds to a bait in the promcap expt

## non transcribed
rnaseq[,trans:=non_1 + non_2]
a <- merge(b2g[,.(id,gene,baitID)],rnaseq[,.(id,trans)],by="id")
setnames(a,"baitID","preyID")
a <- a[,.(z=mean(log2(trans+1))),by="preyID"]
data[["nCD4val"]]  %<>% merge(., a, by="preyID",all.x=TRUE)
data[["nCD4val"]][,y:=as.integer(!is.na(z))]# z only not NA if corresponds to a bait in the promcap expt

## act prom
## z = sum of pp of open/active chromatin states at prey
## y = contact detected by Cao
obj <- (load(file.path(d,"annot","chromhmm-act.RData")))
chromhmm <- eval(as.symbol(obj))
rm(list=obj)
chromhmm[,enhprom:=pmin(1,E4+E5+E6+E7+E8+E9+E10+E11)]
x <- merge(data[["aCD4pchic"]],chromhmm[,.(hindID,enhprom,E14)],
           by.x="preyID",by.y="hindID",all.x=TRUE)
       x[,z:=enhprom]## cuts 2/16m rows, but allows logistic model
       ## x[z=0 & z!=1,z:=NA]
       x[z>0,z:=1]
       x <- merge(x,cao,by=c("baitID","preyID"),all.x=TRUE)
       x[abs(dist)<2e+5,y:=ifelse(is.na(cao),0,1)]
data[["aCD4pchic"]] <- x

## non prom
obj <- (load(file.path(d,"annot","chromhmm-non.RData")))
chromhmm <- eval(as.symbol(obj))
rm(list=obj)
chromhmm[,enhprom:=pmin(1,E4+E5+E6+E7+E8+E9+E10+E11)]
x <- merge(data[["nCD4pchic"]],chromhmm[,.(hindID,enhprom,E14)],
           by.x="preyID",by.y="hindID",all.x=TRUE)
       x[,z:=enhprom]
       ## x[z=0 & z!=1,z:=NA]
       x[z>0,z:=1]
       x <- merge(x,cao,by=c("baitID","preyID"),all.x=TRUE)
       x[abs(dist)<2e+5,y:=ifelse(is.na(cao),0,1)]
data[["nCD4pchic"]] <- x

data <- lapply(data, function(x) x[,.(baitID,preyID,N,chicago,mppc,mid.bait,mid.prey,dist,z,y,experiment)])  %>% rbindlist()
data  <-  addtitles(data)

## source("~/Projects/peaky/common.R")

DATA <- split(data,data$expt)
    ## val:
    ## z - transcription log2(count+1)
    ## y - is a bait in the prom expt
    ## prom:
    ## z - propn of frag in E4-E11 states
    ## y - !is.na(cao.score) for fragment - NB, defined only for dist < 200kb

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
      file=file.path(d,"ptables/dbic-alldata.tex"))

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

ggsave(file.path(d,"figures/fits-alldata.pdf"),height=10,width=8)

rm(DATA,MY,MZ); gc()

################################################################################

##' now do same within runs
bak <- data

head(data)
##' add blocks
## data <- data[order(expt,baitID,preyID),]
data[,csig:=as.numeric(chicago>5),by=c("baitID","experiment")]
data[,start:=c(0,diff(csig))==1 | (csig==1 & c(0,diff(preyID))!=1),by=c("baitID","experiment")]
data[,block:=cumsum(start) * (csig==1),by=c("baitID","experiment")]
    
data[block>0,block_length:=.N,by=c("experiment","baitID","block")]
data[block>0,exp_gamma:=sum(mppc),by=c("experiment","baitID","block")]
data[,bl:=.N,by=c("expt","baitID","block")]
data <- data[block!=0,]
## data[,h:=NULL]
## data <- merge(data,cs,by=c("expt","baitID","preyID"))
## data[,sc:=sqrt(mppc)/sum(sqrt(mppc)),by=c("expt","baitID","block")]



de <- data[grepl("Promoter",expt) & bl>=2,]
summary(de$bl)
quantile(de$bl,seq(0,1,0.1))
de[,blq:=cut(bl,c(2,4,10,20,200),include.lowest = TRUE)]
## de[,blq:=cut(bl,c(2,20,200),include.lowest = TRUE)]
DATA <- split(de,paste(de$expt,de$blq))
names(DATA)

fitter <- function(din, X, Y, n=1000000, f=mod) {
    d <- din[!is.na(din[[Y]]) & !is.nan(din[[Y]]),]
    ## if(nrow(d)>n)
    ##     d <- d[sample(1:nrow(d),n),]
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
    
## mppc useful in longer blocks
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
tab <-rbind(#promy=sumbic(MY),
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
## tabf <- rbind(tabf[,c(4,3)],tabf[,c(2,1)])
tabf

latex(tabf,
      colheads=c("2-4","5-10","11-20","21+"),
      dec=1,
      col.just=c("r","r","r","r"),
      cgroup="Stretch length",
      booktabs=TRUE,
      rowlabel="Model",
      rgroup=c("non","act"),
      ## rgroup=c("a: Promoter: match to Cao et al.","b: Promoter: link to active chromatin",
      ##          "c: Validation: overlap baited promoter","d: Validation: expression at linked promoter"),
      rowname=rep(c("MPPC","CHiCAGO","MPPC + CHiCAGO"),4),
      caption="$\\Delta$BIC from the intercept only model for
whether the prey fragment overlaps active chromatin states defined by \\cite{burren2017chromosome} in non-activated and activated cells.  The best fitting model (lowest $\\Delta$BIC) is highlighted by \\textbf{*}.
In all cases, a robust clustered model was used to account for repeated observations at the prey fragment.",
      label="tab:dbic-stretch", 
      file=file.path(d,"ptables/dbic-stretch.tex"))


## kk <- latex(tabf,
##       colheads=c("2-4","5-10","11-20","21+"),
##       dec=1,
##       col.just=c("r","r","r","r"),
##       booktabs=TRUE,
##       rowlabel="Model (cell)",
##       cgroup="Stretch length",
##       n.cgroup = 4,
##       rgroupTexCmd = "",
##       rgroup=c("\\multicolumn{5}{l}{\\bfseries a: match to Cao et al. (non)}",
##                "\\multicolumn{5}{l}{\\bfseries b: link to active chromatin (non)}",
##                "\\multicolumn{5}{l}{\\bfseries a: match to Cao et al. (act)}",
##                "\\multicolumn{5}{l}{\\bfseries b: link to active chromatin (act)}"),
##       rowname=rep(c("MPPC","CHiCAGO","MPPC + CHiCAGO"),4),
##       caption="$\\Delta$BIC from the intercept only model for four measures of biological plausibility of contacts.  The best fitting model (lowest $\\Delta$BIC) is highlighted by \\textbf{*}.
## \\textbf{a}--\\textbf{b} are defined in full in the Methods.  Briefly,
## \\textbf{a} whether the bait-prey pair corresponds to published CD4$^+$ T cell promoter-enhancer networks\\cite{Cao2017-qx};
## \\textbf{b} whether the prey fragment overlaps active chromatin states defined by \\cite{burren2017chromosome};
## In all cases, a robust clustered model was used to account for repeated observations at the prey fragment.",
##       label="tab:dbic-stretch", 
##       file=file.path(d,"ptables/dbic-stretch-pre.tex"))


