library(data.table)
library(ggplot2)
library(magrittr)
library(cowplot)

d <- "/scratch/wallace/peaky/derived"
list.files(d)
peaky.files <- list.files(d,pattern="peaky.RData")
chic.files <- list.files(d,pattern="chicago")
b=17508
b=17777
b=666461 # AHR
source("~/Projects/peaky/common.R")

files <- list.files(d,pattern="join")
for(nm in names(title)) {
    f <- paste0(nm,"-joined.RData")
    message("!!! ",f)
    (load(file.path(d,f)))
    head(x)
    ncall <- sum(x$chicago>5)
    x[maxres<2,mppc:=0]
    x <- x[order(x$mppc,decreasing=TRUE),]
    x$peaky <- 1:nrow(x)<=ncall
    x$cat <- paste0(ifelse(x$chicago>5,"C","-"),
                    ifelse(x$peaky,"P","-"))
    x <- x[cat!="--",]
    Q <- tapply(abs(x$dist)/1000,x$cat,quantile,seq(0,1,by=0.001),SIMPLIFY=FALSE) %>%
    do.call("cbind",.) %>% as.data.frame() %>% melt(.,"C-")
    colnames(Q)[1] <- "base"

    tt <- table(x$cat)
    nt <- names(tt)
    tt <- sprintf("%4.1f",100*tt/sum(tt))
    x$cat <- factor(x$cat,levels=nt)
    levels(x$cat) <- paste0(nt," (",tt,"%)")
    print(tapply(abs(x$dist),x$cat,summary))
    ggplot(x[cat!="--",],aes(x=abs(dist),col=cat)) + geom_histogram(binwidth=5000) + facet_grid(cat~.,scales="free") + geom_vline(xintercept=5000,col="grey",linetype="dashed") + scale_x_continuous("Distance from bait") + theme(legend.position="none")
    
    #ggplot(Q,aes(x=base,y=value,col=variable)) + geom_path() + geom_abline() 

    ggsave(file=file.path(d,"../figures",paste0(nm,"-distance.pdf")),height=6,width=6)
}
