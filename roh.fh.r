#
system("plink --bfile ARIC.commonROH_lite --linear --recessive --pheno disc.resid.pheno --out ARIC.fh.rec")
system("plink --bfile ARIC.roh --linear --recessive --pheno disc.resid.pheno --out ARIC.roh.rec")

assoc <- read.table("ARIC.fh.rec.assoc.linear",header=TRUE)
nsub <- 8749

h2.dist <- vector()
for(i in 1:1000){
temp <- rnorm(nsub,mean=mean(assoc$BETA,na.rm=TRUE),sd=sd(assoc$BETA,na.rm=TRUE))
test <- t.test(temp,mu=0)
h2.dist[i] <- test$statistic
print(i)
}

t2 <- abs(mean(h2.dist))^2
h2 <- t2/(t2+nsub)
# Total model:

a.ind <- read.table("ARIC.commonROH_lite_snp65.hom.indiv",header=TRUE)
phen <- read.table("disc.resid.pheno",header=FALSE)
names(phen) <- c("FID","IID","PHE")
a <- merge(a.ind,phen,by="IID")
het <- read.table("ARIC.commonROH_lite.het",header=TRUE)
a <- merge(a,het,by="IID")
a$set <- rep("a",nrow(a))
a <- a[order(a$KB),]
a$FROH <- a$KB/2160000
a <- a[duplicated(a$IID)==FALSE,]


g.ind <- read.table("GNHS.commonROH_lite_snp65.hom.indiv",header=TRUE)
g <- merge(g.ind,phen,by="IID")
g.het <- read.table("GNHS.commonROH_lite.het",header=TRUE)
g <- merge(g,g.het,by="IID")
g$set <- rep("g",nrow(g))
g <- g[order(g$KB),]
g <- g[which(log(g$KB)>8),] # remove defective subjects
g$FROH <- g$KB/2160000
g <- g[duplicated(g$IID)==FALSE,]

wg.ind <- read.table("WHI.GRU.commonROH_lite_snp65.hom.indiv",header=TRUE)
wg <- merge(wg.ind,phen,by="IID")
wg.het <- read.table("WHI.GRU.commonROH_lite.het",header=TRUE)
wg <- merge(wg,wg.het,by="IID")
wg$set <- rep("wg",nrow(wg))
wg <- wg[order(wg$KB),]
wg$FROH <- wg$KB/2160000

wn.ind <- read.table("WHI.NPU.commonROH_lite_snp65.hom.indiv",header=TRUE)
wn <- merge(wn.ind,phen,by="IID")
wn.het <- read.table("WHI.NPU.commonROH_lite.het",header=TRUE)
wn <- merge(wn,wn.het,by="IID")
wn$set <- rep("wn",nrow(wn))
wn <- wn[order(wn$KB),]
wn$FROH <- wn$KB/2160000

w <- rbind(wn,wg)
w <- w[duplicated(w$IID)==FALSE,]

a.t <- rbind(a,g)
b.t <- rbind(wg,wn)
t.t <- rbind(a.t,b.t)
t.t$FROH <- t.t$KB/2160000

t.t <- t.t[duplicated(t.t$IID)==FALSE,]


summary(lm(t.t$PHE.y~t.t$FROH))
summary(lm(PHE.y~F,data=t.t))

summary(lm(a$PHE.y~a$FROH))
summary(lm(a$PHE.y~a$F))

summary(lm(w$PHE.y~w$FROH))
summary(lm(w$PHE.y~w$F))

summary(lm(g$PHE.y~g$FROH))
summary(lm(g$PHE.y~g$F))

install.packages("zoo")
library('zoo')

o.index <- rev(order(g$FROH))
plot(t.t$F[o.index],ylim=c(-.5,.5),xlab="sorted by FROH",ylab="%",main="Total")
points(t.t$FROH[o.index],col="red")


phen <- rollmean(t.t$PHE.y[o.index],500)
points(phe,col="blue")



mod4 <- lm(scale(proh)~scale(as.numeric(phen$SCORE)))
d4 <- summary(mod4)
plot(scale(proh),scale(phen$SCORE),xlab="proh_avg",ylab="SCORE",main=paste("est.= ",round(d4$coefficients[2,1],4)," t= ",round(d4$coefficients[2,3],2)," p= ",round(d4$coefficients[2,4],3),sep=""))
abline(mod4)
