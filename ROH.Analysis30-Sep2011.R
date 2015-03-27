
# START HERE FOR IQ
################################################
#23) ROH mapping approach to finding if ROHs in particular regions are predictive of case-control status
#1MB bins
setwd("/STATGEN/home/howrigan/ROH_IQ")
# Read in data with runs and covariate information:
ncs <- count.fields("ROH.matrix.data")
runs.data <- read.table("ROH.matrix.data",header=TRUE,fill=TRUE,colClasses=c(rep('factor',4),rep('numeric',ncs[1]-4)))# runs start at col 5
runs.only <- as.matrix(runs.data[,5:ncol(runs.data)]) # generate matrix without covariates; 4134 columns
tot.runs <- apply(runs.only,2,sum,na.rm=TRUE)
case.runs <- apply(runs.only[runs.data$PHE==2,],2,sum,na.rm=TRUE)
control.runs <- apply(runs.only[runs.data$PHE==1,],2,sum,na.rm=TRUE)

#Read in key
key <- read.table("ROH.matrix.key",header=FALSE)

# Read in map data of regression results on actual data:
ncs2 <- count.fields("ROH.map.data")
map.data <- read.table("ROH.map.data",header=FALSE,colClasses=c('character',rep('numeric',ncs2[1]-1))) # read in un-permuted results
map.data2 <- as.matrix(map.data[,2:ncol(map.data)])
dimnames(map.data2) <- list(map.data[,1],NULL)
position <- 1:ncol(map.data2)
chr <- as.numeric(key[1,3:ncol(key)])
rel.position <- as.numeric(key[2,3:ncol(key)])
#chr[chr==0] <- rel.position[chr==0] <- NA
nlog10pval <- -1*log10(map.data2[4,])
map.data3 <- rbind(map.data2,nlog10pval,position,rel.position,chr,tot.runs,case.runs,control.runs)
#map.data3 <- map.data3[, -ncol(map.data3)]
map.data4 <- as.data.frame(t(map.data3))
chr.ends <- c(max(which(chr==1)),max(which(chr==2)),max(which(chr==3)),max(which(chr==4)),max(which(chr==5)),max(which(chr==6)),max(which(chr==7)),max(which(chr==8)),max(which(chr==9)),max(which(chr==10)),max(which(chr==11)),max(which(chr==12)),max(which(chr==13)),max(which(chr==14)),max(which(chr==15)),max(which(chr==16)),max(which(chr==17)),max(which(chr==18)),max(which(chr==19)),max(which(chr==20)),max(which(chr==21)),max(which(chr==22)))


#look at data
hist(map.data4$Beta1_P_val,breaks=200)  #overly many small p-values
(small.pval <- map.data3[,order(map.data3['Beta1_P_val',])][,1:100])  #look at the 100 smallest p-values
summary(small.pval['Beta1',] > 0) #of the 100 smallest p-values, 41 of the associated slopes are positive (ROHs are risk factors for SZ)
summary(map.data3['Beta1',] > 0)  #of all the slopes, 1386 are positive and 1303 are negative or 49.5%

# UNSURE WHAT NEXT LINE MEANS (comment by Matt S.)
pbinom(q=95,size=100,prob=3512/(3512+1718),lower.tail=FALSE,log.p=FALSE) #probability of 95/100 occurring by chance given 67% base-rate: 1.24e-12

#Read in permuted data to get genome-wide correction
permute.data.orig <- read.table("ROH.map.permute",header=FALSE,colClasses=c('character',rep('numeric',ncs2[1]-1))) # read in results of permutation
permute.data <- as.matrix(permute.data.orig[,2:ncol(permute.data.orig)])
dimnames(permute.data) <- list(permute.data.orig[,1],NULL)

permuted.beta1 <- permute.data[seq(2,nrow(permute.data),by=8),]
permuted.pval <- permute.data[seq(4,nrow(permute.data),by=8),]
storage.mode(permuted.beta1) <- "numeric"
storage.mode(permuted.pval) <- "numeric"

#permute.data2 <- rbind(permute.data2,position,chr)
#permute.data2 <- permute.data2[,- ncol(permute.data2)]

#Find lowest beta1 p-values for each permutation (row)
lowest.pvals <- rep(1.1,1000)
for (i in 1:1000){lowest.pvals[i] <- min(permuted.pval[i,],na.rm=TRUE)}

#get 95% cut-off for these values
(pval.cutoff <- lowest.pvals[order(lowest.pvals)][50]) 
(pval.cutoff2 <- lowest.pvals[order(lowest.pvals)][100]) 

#useful function
mid <- function(x1,x2){abs(x2-x1)/2 + min(x1,x2)}
mid2 <- function(x1,x2){(x2-x1)/2 + x1}
mid.vec <- function(x){
  z <- vector(length=length(x)-1)
 for (i in 1:length(z)){z[i] <- mid(x[i+1],x[i])}
 return(z)
}


na.exclude(map.data4[-1*log10(map.data4$Beta1_P_val) > 3,])


#PLOT
pdf("Figure5-ROH.mapping.pdf",height=7,width=11)
#x11(height=7,width=11)  #for formatting trials

op <- par(fig=c(0,1,.3,1),mar=c(3,4,2,.5),xpd=FALSE)

#manhattan plot
l.cols <- c('red','blue')[(map.data3['Beta1',]<0)+1]
plot(-1*log10(map.data3['Beta1_P_val',]),type='p',col=l.cols,ylab="- log10 p-value of case-control difference in # ROHs",xlab="Chromosome & Mb position within chromosomes",axes=FALSE,ylim=c(0,4))

axis(1,1:22,at=mid.vec(c(0,chr.ends)))
axis(2)
#abline(h=-1*log10(.05/2895),lty=2)  #bonferoni correction; 1 location is sig
abline(h=-1*log10(pval.cutoff)) #genome-wide correction; 0 locations are sig
abline(h=-1*log10(pval.cutoff2),lty=2) #genome-wide correction; 11 locations are sig
abline(v=chr.ends+.5)
legend(x=1950,y=4,legend=c('ROHs more common in cases','ROHs more common in controls','Genomewide significance threshold','Suggestive threshold'),pch=1,col=c('red','blue','white','white'),bg='white')
lines(x=c(1970,2090),y=c(3.3,3.3),lty=1)
lines(x=c(1970,2090),y=c(3.1,3.1),lty=2)



#plot of # of ROHs per region
par(new=TRUE)
par(fig=c(0,1,0,.45),mar=c(5,4,5,.5))
plot((tot.runs/21847)*100,type='l',ylim=c(0,1.6),axes=FALSE,ylab="% having ROH",xlab="Chromosome & Mb position within chromosomes",col='green')
axis(1,1:22,at=mid.vec(c(0,chr.ends)))
axis(2)
abline(v=chr.ends+.5)

max.pct <- round(max((tot.runs/21847)*100),2)
max.loc <- which(tot.runs==max(tot.runs))
par(xpd=TRUE)
text(max.loc+145,1.8,max.pct)
arrows(max.loc+80,1.7,max.loc,1.61,length=.05)
par(op)

dev.off()
################################################










################################################
#24) ROH mapping approach to finding if ROHs in particular regions are predictive of case-control status
#500 kb bins - NO COVARIATES

setwd("~/HD3/SNP.Homozygosity/PGC.SZ/mapping/Logistic_test4")

# Read in data with runs and covariate information: (this is same as above, so no need to read in)
#ncs <- count.fields("ROH.matrix.data")
#runs.data <- read.table("ROH.matrix.data",header=TRUE,fill=TRUE,colClasses=c(rep('factor',6),rep('numeric',ncs[1]-6)))# runs start at col 27
#runs.only <- as.matrix(runs.data[,27:ncol(runs.data)]) # generate matrix without covariates; 5742 columns
#tot.runs <- apply(runs.only,2,sum,na.rm=TRUE)
#case.runs <- apply(runs.only[runs.data$PHE==2,],2,sum,na.rm=TRUE)
#control.runs <- apply(runs.only[runs.data$PHE==1,],2,sum,na.rm=TRUE)

#Read in key (same as above)
#key <- read.table("ROH.matrix.key",header=FALSE)

# Read in map data of regression results on actual data:
ncs2NC <- count.fields("ROH.map.noCOV.data")
map.dataNC <- read.table("ROH.map.noCOV.data",header=FALSE,colClasses=c('character',rep('numeric',ncs2[1]-1))) # read in un-permuted results
map.data2NC <- as.matrix(map.dataNC[,2:ncol(map.dataNC)])
dimnames(map.data2NC) <- list(map.dataNC[,1],NULL)
positionNC <- 1:ncol(map.data2NC)
chrNC <- as.numeric(key[1,3:ncol(key)])
rel.positionNC <- as.numeric(key[2,3:ncol(key)])
#chr[chr==0] <- rel.position[chr==0] <- NA
nlog10pvalNC <- -1*log10(map.data2NC[4,])
map.data3NC <- rbind(map.data2NC,nlog10pvalNC,positionNC,rel.positionNC,chrNC,tot.runs,case.runs,control.runs)
#map.data3 <- map.data3[, -ncol(map.data3)]
map.data4NC <- as.data.frame(t(map.data3NC))
chr.ends <- c(498,982,1375,1757,2117,2458,2774,3067,3346,3615,3884,4147,4375,4586,4784,4962,5119,5272,5398,5524,5618,5716)


#look at data
hist(map.data4NC$Beta1_P_val,breaks=200)  #overly many small p-values
(small.pvalNC <- map.data3NC[,order(map.data3NC['Beta1_P_val',])][,1:100])  #look at the 100 smallest p-values
summary(small.pvalNC['Beta1',] > 0) #of the 100 smallest p-values, 91 of the associated slopes are positive (ROHs are risk factors for SZ)
summary(map.data3NC['Beta1',] > 0)  #of all the slopes, 3352 are positive and 1878 are negative, or 64%
pbinom(q=91,size=100,prob=3352/(3352+1878),lower.tail=FALSE,log.p=FALSE) #probability of 95/100 occurring by chance given 67% base-rate: 1.24e-12

#Read in permuted data to get genome-wide correction
permute.data.origNC <- read.table("ROH.map.noCOV.permute",header=FALSE,colClasses=c('character',rep('numeric',ncs2[1]-1))) # read in results of permutation
permute.dataNC <- as.matrix(permute.data.origNC[,2:ncol(permute.data.origNC)])
dimnames(permute.dataNC) <- list(permute.data.origNC[,1],NULL)

permuted.beta1NC <- permute.dataNC[seq(2,nrow(permute.dataNC),by=8),]
permuted.pvalNC <- permute.dataNC[seq(4,nrow(permute.dataNC),by=8),]
storage.mode(permuted.beta1NC) <- "numeric"
storage.mode(permuted.pvalNC) <- "numeric"

#permute.data2 <- rbind(permute.data2,position,chr)
#permute.data2 <- permute.data2[,- ncol(permute.data2)]

#Find lowest beta1 p-values for each permutation (row)
lowest.pvalsNC <- rep(1.1,1000)
for (i in 1:1000){lowest.pvalsNC[i] <- min(permuted.pvalNC[i,],na.rm=TRUE)}

#get 95% cut-off for these values
(pval.cutoffNC <- lowest.pvalsNC[order(lowest.pvalsNC)][50]) 
(pval.cutoff2NC <- lowest.pvalsNC[order(lowest.pvalsNC)][100]) 


#PLOT
layout(matrix(c(1,1,2,2),nrow=2,byrow=TRUE),heights=c(2,.75))

#manhattan plot
l.cols <- c('red','blue')[(map.data3NC['Beta1',]<0)+1]
plot(-1*log10(map.data3NC['Beta1_P_val',]),type='p',col=l.cols,ylab="- log10 p-value of case-control difference in # ROHs",xlab="Chromosome & Mb position within chromosomes",axes=FALSE,ylim=c(0,4))

axis(1,1:22,at=mid.vec(c(0,chr.ends)))
axis(2)
#abline(h=-1*log10(.05/2895),lty=2)  #bonferoni correction; 1 location is sig
abline(h=-1*log10(pval.cutoffNC)) #genome-wide correction; 0 locations are sig
abline(h=-1*log10(pval.cutoff2NC),lty=2) #genome-wide correction; 11 locations are sig
abline(v=chr.ends+.5)
legend(x=1500,y=4,legend=c('ROHs more common in cases','ROHs more common in controls','Bonferoni threshold','Permutation threshold'),pch=1,col=c('red','blue','white','white'),bg='white')



#plot of # of ROHs per region
plot((tot.runs/21847)*100,type='l',ylim=c(0,1.6),axes=FALSE,ylab="% having ROH",xlab="Chromosome & Mb position within chromosomes",col='green')
axis(1,1:22,at=mid.vec(c(0,chr.ends)))
axis(2)
abline(v=chr.ends+.5)

max.pct <- round(max((tot.runs/21847)*100),2)
max.loc <- which(tot.runs==max(tot.runs))
op <- par(xpd=TRUE)
text(max.loc+65,1.85,max.pct)
arrows(max.loc+55,1.75,max.loc,1.61,length=.05)
par(op)

#the below makes the graph look too messy and is left out:
#op <- par(new=TRUE)
#plot(control.runs,type='l',col='blue',ylim=c(0,400),axes=FALSE)
#par(new=TRUE)
#plot(case.runs,type='l',col='red',ylim=c(0,400),axes=FALSE)
#par(op)
################################################









################################################
#25) Look at specific regions from above


#Look at two most significant regions
sig.data <- map.data4[map.data4$Beta1_P_val < pval.cutoff2,]
(sig.data2 <- na.omit(sig.data))
#Rerun the above Froh analyses after removign these two regions by manually rerunning from above!


#look at some spikes where lots of ROHs are
#Chr2
map.data4[map.data4$chr==2 & map.data4$tot.runs>150,]  #lactose gene is 136.55-136.59; the last two hits are at 136 & 137 Mb
#Chr 6
map.data4[map.data4$chr==6,][1:75,]


#why are there so few ROHs under the LAC gene? Gibson (2006) found 27% for 500kb around LAC. This is probably bc we prune for LD.
HOM.lac <- HOM[HOM$CHR==2 & HOM$POS2 < 138e6 & HOM$POS1 > 135e6,]
#see how dense the SNPs are at various locations
summary(HOM.lac$DENSITY);summary(HOM$DENSITY) #32 kb/SNP in ROHs this region vs. 19.5 kb/SNP across genome
summary(HOM.lac);summary(HOM)
setwd("~/HD3/SNP.Homozygosity/PGC.SZ/DATA/AUG2011.Imputed")  #data from imputed SNPs
bim <- read.table("Fin.Imputed1a7.bim",header=FALSE)
bim2 <- bim[bim$V1==2 & bim$V4 > 136e6 & bim$V4 < 137e6,]
2.77e9/nrow(bim)       #1 snp per 29.2 kb across entire genome
1e6/nrow(bim2)         #1 snp per 91 kb in the LAC region; over 3 x less dense - this is why so few ROHs here relative to expectation of 27%
#you could also redo the PLINK --hom analysis using only --homozyg-kb 500 argument to show that ~27% of our sample is indeed homozygous here

bim6 <- bim[bim$V1==6 & bim$V4 > 26e6 & bim$V4 < 34e6,]
8e6/nrow(bim6)       #1 snp per 13.8 kb across HLA
2.77e9/nrow(bim)       #1 snp per 29.2 kb across entire genome


#Look at MHC region
hla <- HOM_I[HOM_I$CHR==6 & HOM_I$POS1 < 30.5e6 & HOM_I$POS2 > 27e6,]



################################################









################################################
#26) Look at overlap between CNVs and ROHs

cnv.gain <- read.table("~/HD3/SNP.Homozygosity/PGC.SZ/CNV.MGS.Calls/output_gain_ea_broad_final_cnvsonly_colorado.txt",header=TRUE)
cnv.nongain <- read.table("~/HD3/SNP.Homozygosity/PGC.SZ/CNV.MGS.Calls/output_nongain_broad_final_cnvsonly_colorado.txt",header=TRUE)
cnvs <- rbind(cnv.gain,cnv.nongain)
dels <- cnvs[cnvs$CN==1,] 
dels$length <- dels$END - dels$START
dels$IID2 <- paste("NG-",dels$ID,sep='')
summary(dels$length)  #median length of del = 10 kb;
summary(HOM_I$KB)       #median length of ROH = 2084 kb

HOM.mgs <- HOM_I[HOM_I$dataset=='mgs2',]

COMB <- merge(HOM.mgs,dels,by.x=c('IID','CHR'),by.y=c('IID2','CHR'),all=TRUE)
COMB2 <- na.omit(COMB)
COMB2$del.midpt <- (COMB2$END-COMB2$START)/2 + COMB2$START
COMB2$roh.midpt <- (COMB2$POS2-COMB2$POS1)/2 + COMB2$POS1

COMB3 <- COMB2[(COMB2$del.midpt > COMB2$POS1 & COMB2$del.midpt < COMB2$POS2) | (COMB2$roh.midpt > COMB2$START & COMB2$roh.midpt < COMB2$END),]
COMB3$LENGTH.del <- COMB3$LENGTH/1000
COMB3$POS1.roh <- COMB3$POS1/1000
COMB3$POS2.roh <- COMB3$POS2/1000
COMB3$START.cnv <- COMB3$START/1000
COMB3$END.cnv <- COMB3$END/1000
COMB3 <- COMB3[order(COMB3$IID),]
COMB3[,c('IID','CHR','POS1.roh','POS2.roh','START.cnv','END.cnv','KB','LENGTH.del','del.midpt','roh.midpt')]

COMB3.shortdels <- COMB3[COMB3$LENGTH.del < COMB3$KB,]
COMB3.longdels <- COMB3[COMB3$LENGTH.del > COMB3$KB,]

COMB3.shortdels[1:125,c('IID','CHR','POS1.roh','POS2.roh','START.cnv','END.cnv','KB','LENGTH.del','del.midpt','roh.midpt','LENGTH.del.comb')]

ind1 <- which(diff(COMB3.shortdels$KB)==0)
ind2 <- ind1+1
COMB3.shortdels$LENGTH.del.comb <- COMB3.shortdels$LENGTH.del
COMB3.shortdels$LENGTH.del.comb[ind1] <- COMB3.shortdels$LENGTH.del.comb[ind1]+COMB3.shortdels$LENGTH.del.comb[ind2]

COMB4.shortdels <- COMB3.shortdels[- (ind1[which((diff(ind1) !=1))]+1),]

ind1 <- which(diff(COMB4.shortdels$KB)==0)
ind2 <- ind1+1
COMB4.shortdels$LENGTH.del.comb <- COMB4.shortdels$LENGTH.del
COMB4.shortdels$LENGTH.del.comb[ind1] <- COMB4.shortdels$LENGTH.del.comb[ind1]+COMB4.shortdels$LENGTH.del.comb[ind2]


COMB5.shortdels <- COMB4.shortdels[- (ind1[which((diff(ind1) !=1))]+1),]

ind1 <- which(diff(COMB5.shortdels$KB)==0)
ind2 <- ind1+1
COMB5.shortdels$LENGTH.del.comb <- COMB5.shortdels$LENGTH.del
COMB5.shortdels$LENGTH.del.comb[ind1] <- COMB5.shortdels$LENGTH.del.comb[ind1]+COMB5.shortdels$LENGTH.del.comb[ind2]
COMB5.shortdels$len.roh.remain <- COMB5.shortdels$KB - COMB5.shortdels$LENGTH.del.comb

COMB5.shortdels[,c('IID','CHR','POS1.roh','POS2.roh','START.cnv','END.cnv','KB','LENGTH.del','del.midpt','roh.midpt','LENGTH.del.comb')]

removed.rohs <- COMB5.shortdels[COMB5.shortdels$len.roh.remain < 500,]  #7 removed here
COMB3.longdels #3 removed here
nrow(HOM.mgs)
#So with rule to remove any ROHs that are < 500kb after removing any overlapping deletion portions, we remove 10 out of a total of 6480
#this is the same rule used in McQuillan et al, 2008
10/6480  #0.15% of ROHs removed, this compares to <0.3% in McQuillan et al, 2008

sum(COMB5.shortdels$LENGTH.del)/sum(COMB5.shortdels$KB)
################################################











################################################
#27) Look at distribution of ROH lengths and Froh


#PLOT
pdf("Figure1-ROH.Froh.Distns.pdf",height=5,width=10)
#x11(height=5,width=10)  #for formatting trials

op <- par(fig=c(0,.5,0,1))
hist(HOM_I$KB/1000,breaks=100,main="ROH Lengths",xlab="Mb of ROHs",col='lightblue')
par(new=TRUE)
par(fig=c(.5,1,0,1))
hist(MAIN2$tot2.1.imp[MAIN2$tot2.1.imp < .0625],breaks=100,main="Distribution of Froh < .0625",xlab="Froh",col='red')
par(op)
dev.off()















