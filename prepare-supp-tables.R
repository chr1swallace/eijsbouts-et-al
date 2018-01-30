
library(data.table)
#ap <- readRDS("/scratch/wallace/peaky/tables/DT_aCD4pchic_double5e6.rds")
#av <- readRDS("/scratch/wallace/peaky/tables/DT_aCD4val_double5e6.rds")
#np <- readRDS("/scratch/wallace/peaky/tables/DT_nCD4pchic_double5e6.rds")
#x <- readRDS("/scratch/wallace/peaky/tables/DT_nCD4val_quadruple5e6.rds")
#head(x)

## new data location
##x <- readRDS("/mrc-bsu/scratch/cew54/peaky/summary/tables/DT_aCD4val_20e6_corr.rds")

d <- "/mrc-bsu/scratch/cew54/peaky/summary/tables" #"/scratch/wallace/peaky/tables"
files <- c(actprom="DT_aCD4pchic_20e6_corr.rds",
           actval="DT_aCD4val_20e6_corr.rds",
           nonprom="DT_nCD4pchic_20e6_corr.rds",
           nonval="DT_nCD4val_20e6_corr.rds")

for(i in seq_along(files)) {
    f <- files[i]
    nm <- names(files)[i]
    of <- file.path(d,"../derived", paste0("suppdata-",nm,".csv"))
    ofz <- file.path(d,"../derived", paste0("suppdata-",nm,".csv.gz"))
    if(file.exists(of)) unlink(of)
    if(file.exists(ofz)) unlink(ofz)
    if(file.exists(of) || file.exists(ofz)) {
        message("\n\nskipping ",nm,": already exists ",of)
        next
    }
    message("\n\nreading ",nm," from ",f)

    x <- readRDS(file.path(d,f))

    message("subsetting columns")
    setnames(x,
         c("rjmcmc_pos_g-4.69897","rjmcmc_g-4.69897",
           "beta_mean_g-4.69897", "predicted_g-4.69897"),
         c("mppc","mppi","beta.post","peaky.pred"))
    x[,maxcor:=pmax(corr_5v5,corr_10v10)]
    use <- which(x$maxcor>0.75)
    if("score" %in% names(x)) {
        raw <- x[use,.(baitID,preyID,N,residual,mppc,beta.post,score)]
        setnames(raw,"score","chicago")
    } else {
        (load(file=file.path(d,"../derived", paste0(nm,"-chicago.RData"))))
        raw <- x[use,.(baitID,preyID,N,residual,mppc,beta.post)]
        raw <- merge(raw,chic,by=c("baitID","preyID"))
    }
    maxchic <- raw[,max(chicago),by="baitID"]
    maxchic <- maxchic[V1>=5,]
    raw <- raw[baitID %in% maxchic$baitID,]
    message("writing to ",of)
    fwrite(raw,file=of)
    ## message("zipping")
    ## system(paste("gzip",of))
}

(load(file.path(d,"b2gene.RData")))
head(int.genes)
fwrite(int.genes,file=file.path(d,"../derived", paste0("bait2gene.csv")))

d2 <- "/mrc-bsu/scratch/cew54/peaky/summary/derived"
(load(file.path(d2,"hind-position.RData"))); hind[,hindID:=as.integer(hindID)]
head(hind)
fwrite(hind,file=file.path(d,"../derived", paste0("hind-positions.csv")))

system(paste0("cd ",d,"/../derived && tar zcvf supp-data.tgz bait2gene.csv hind-positions.csv suppdata-*.csv"))
