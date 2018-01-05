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
results <- vector("list",length(title))
names(results) <- names(title)
for(nm in names(title)) {
   message(" !!! ",nm) 
   fnm <- file.path(d,paste0(nm,"-joined.RData"))
   (load(fnm))
   setkey(x,baitID)
   x <- x[order(baitID,dist),]
   width=250
   windows = data.table(w_start = seq(1,max(1,nrow(x)-(width-1)),by=width))
   windows[,w_end:=w_start+width]
   windows[,bstart:=x$baitID[w_start]]
   windows[,bend:=x$baitID[w_end]]
   windows <- windows[bstart==bend,]
   dim(windows)
   stats = apply(windows,1,function(w){x[w['w_start']:w['w_end'],
                                         .(mn=mean(residual),v=var(residual))]})  %>% rbindlist()
   results[[nm]] <- cbind(windows,stats,nm=nm)
}

results <- rbindlist(results)

m <- melt(results,c("nm","w_start","w_end","bstart","bend"))

results$lab <- "mean"
pm <- ggplot(results,aes(x=mn,y=..density..,fill=nm)) + geom_histogram() + facet_grid(.~nm) + geom_vline(xintercept=0,col="black",linetype="dashed") + xlab("Mean") + theme(legend.position="none") + ylab("Density")
results$lab <- "variance"

f <- function(x) (x)
pv <- ggplot(results,aes(x=f(v),y=..density..,fill=nm)) + geom_histogram() + facet_grid(.~nm) + geom_vline(xintercept=f(1),col="black",linetype="dashed") + xlab("Variance") + theme(legend.position="none") + ylab("Density")
pv

plot_grid(pm,pv,nrow=2)
ggsave(file.path(d,"../figures","omega-residuals-allbaits.pdf"),height=6,width=8)

