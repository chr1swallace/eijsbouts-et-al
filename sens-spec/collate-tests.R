#!/usr/bin/env Rscript
library(magrittr)
library(data.table)
od <- "~/bsu/peaky-sens-spec"
(dirs <- list.files(od,full=TRUE,pattern="^N.*[012]$"))

d <- dirs[1]
head(list.files(d))
for(d in dirs) {
    message(d)
    z <- fread(paste0("cat ",d,"/*"))
    z <- z[baitID!="baitID",]
    z[,mppc:=as.numeric(mppc)]
    z[,baitID:=as.integer(baitID)]
    z[,preyID:=as.integer(preyID)]
    setnames(z,"mppc",basename(d))
    if(d==dirs[[1]]) {
        x <- z
    } else {
        x <- merge(x,z,by=c("baitID","preyID"))
    }
}

head(x)

m <- melt(x,c("baitID","preyID"))

library(ggplot2)
library(cowplot)
plot.mppc.dn <- ggplot(m,aes(x=value)) + geom_histogram(col="grey",fill="black") + scale_y_log10() + labs(x="MPPC",y="total count (log scale)")

## ggplot(m,aes(x=value,y=stat(count),col=variable)) +
##   geom_histogram(bins=10,position=position_dodge2(padding=0.1)) +
##   scale_y_log10() 

## ggplot(m,aes(x=value,y=stat(count),col=variable)) +
##   geom_density() +
##   scale_y_sqrt() 

## ggplot(m,aes(x=value,col=variable)) +
##   stat_ecdf(geom = "step") + scale_x_sqrt() + scale_y_sqrt()

## ks.test(m[variable=="N1.N2",]$value,m[variable=="N0",]$value)

## ggplot(m,aes(x=value,y=stat(density),col=variable)) + geom_freqpoly(binwidth) + scale_y_log10()


m[,mean(value<0.01),by="variable"] # 
m[,mean(value<0.1),by="variable"] # 
m[,mean(value<0.01)] # 83%
m[,mean(value<0.1)] # 99%


library(pROC)
## roc_obj <- roc(category, prediction)
## auc(roc_obj)
ss <- function(test,ref,lo=0.01,hi=0.05) {
    z <- copy(x)
    z <- z[ z[[ref]] < lo | z[[ref]] > hi, ]
    z[,category:=as.numeric(z[[ref]]>hi)]
    o <- roc(z$category, z[[test]])
    data.table(sens=o$sensitivities,spec=o$specificities,thr=o$thresholds,auc=o$auc,
               test=test,ref=ref,lo=lo,hi=hi)
    ## auc(o)
    ## lapply(c(seq(0,0.02,by=0.001),seq(0.03,1,by=0.01)), function(thr) {
    ##     sens=sum(z[[test]]>thr & z[[ref]]>hi)/sum(z[[ref]]>hi)
    ##     spec=sum(z[[test]]<=thr & z[[ref]]<lo)/sum(z[[ref]]<lo)
    ##     data.table(thr=thr,sens=sens,spec=spec,lo=lo,hi=hi,ref=ref,test=test)
    ## })  %>% rbindlist()
}

library(parallel)
options(mc.cores=4)
s0 <- mclapply(c(0.01,0.02,0.05,0.1), function(h) ss("N1.N2","N0",lo=0.01,hi=h))  %>%  rbindlist()
s1 <- mclapply(c(0.01,0.02,0.05,0.1), function(h) ss("N0.N2","N1",lo=0.01,hi=h))  %>%  rbindlist()
s2 <- mclapply(c(0.01,0.02,0.05,0.1), function(h) ss("N0.N1","N2",lo=0.01,hi=h))  %>%  rbindlist()
r0 <- mclapply(c(0.01,0.02,0.05,0.1), function(h) ss("N0","N1.N2",lo=0.01,hi=h))  %>%  rbindlist()
r1 <- mclapply(c(0.01,0.02,0.05,0.1), function(h) ss("N1","N0.N2",lo=0.01,hi=h))  %>%  rbindlist()
r2 <- mclapply(c(0.01,0.02,0.05,0.1), function(h) ss("N2","N0.N1",lo=0.01,hi=h))  %>%  rbindlist()
## r <- rbind(r0,r1,r2)
## ggplot(r,aes(x=1-spec,y=sens,col=factor(hi))) + geom_path() + facet_wrap(~test)

r01 <- mclapply(c(0.01,0.02,0.05,0.1), function(h) ss("N0","N1",lo=0.01,hi=h))  %>%  rbindlist()
r02 <- mclapply(c(0.01,0.02,0.05,0.1), function(h) ss("N0","N2",lo=0.01,hi=h))  %>%  rbindlist()
r10 <- mclapply(c(0.01,0.02,0.05,0.1), function(h) ss("N1","N0",lo=0.01,hi=h))  %>%  rbindlist()
r12 <- mclapply(c(0.01,0.02,0.05,0.1), function(h) ss("N1","N2",lo=0.01,hi=h))  %>%  rbindlist()
r20 <- mclapply(c(0.01,0.02,0.05,0.1), function(h) ss("N2","N0",lo=0.01,hi=h))  %>%  rbindlist()
r21 <- mclapply(c(0.01,0.02,0.05,0.1), function(h) ss("N2","N1",lo=0.01,hi=h))  %>%  rbindlist()

## plot faceted by test/ref, different cutoffs in same plot
l <- list(s0,s1,s2,r0,r1,r2,r01,r02,r10,r12,r20,r21)
r <- rbindlist(l)
rp <- r[thr==hi,]
r[,d:=(sens - (1-spec))/sqrt(2)]
r[,dmax:=max(d),by=c("ref","test","hi")]
rd <- r[dmax==d,]
ggplot(r,aes(x=1-spec,y=sens,col=factor(hi))) +
  geom_path() +
  ## geom_point(data=rp) +
  ## geom_text(aes(label=thr),data=rd) +
  facet_wrap(~test+ref,labeller=label_both)

rd[,.(mean(thr),mean(sens),mean(spec)),by=hi]
ggplot(rd,aes(x=hi,y=thr)) + geom_point() + geom_abline()

## plot faceted by cutoff, different test/ref in same plot
r[,nref:=ifelse(nchar(ref)==2,3,6)]
r[,ntest:=ifelse(nchar(test)==2,3,6)]
r[,lab:=paste0("R=",nref," T=",ntest)]
r[,grp:=paste(ref,test)]
r[,cutoff:=paste0("cutoff: ",hi)]
plot.roc <- ggplot(r,aes(x=1-spec,y=sens,col=lab,group=grp)) +
  geom_path() +
  ## geom_point(data=rp) + 
  facet_wrap(~cutoff) +
  background_grid() +
  scale_colour_discrete("Ref/Test sample size") +
  labs(x="1-specificity",y="sensitivity") +
  theme(legend.position=c(0.55,0.2))
plot.roc

auc <- unique(r,by=c("ref","test"))
plot.auc <- ggplot(auc,aes(x=lab,y=auc,col=lab)) + geom_boxplot() + geom_point() + facet_wrap(~cutoff) + labs(x="Ref/Test sample size",y="AUC") + theme(legend.position="none") + background_grid()
plot.auc

plot.left <- plot_grid(plot.mppc.dn,plot.auc,labels=c("a","c"),nrow=2)
plot_grid(plot.left,plot.roc,labels=c("","b"),ncol=2,rel_widths=c(0.4,0.6))
ggsave(file.path(od,"../peaky/figures/sensspec.pdf"),height=10,width=10)

