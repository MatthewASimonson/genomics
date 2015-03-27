



################################################
#1) Import functions needed for post-processing
#source("/home/mmkeller/HD3/SNP.Homozygosity/ROH.ANALYSIS.FUNCTIONS5.R")
source("~/HD3/SNP.Homozygosity/ROH.ANALYSIS.FUNCTIONS6.R")

require(lme4)

require(sciplot)
################################################





################################################
#2) Grab the RAW PGC data -
# *_a12 is using the light pruning (VIF 10), 60 SNP, 0 het thresholds
# *_a9 is using moderate pruning (VIF 2), 40 SNP, 0 het thresholds
#the script for making the ROHs is at lisa.sara:/home/matthew/ROH/SHELL.FILES/*_A.sh
#the ROH data is located at statgen2:~/HD3/SNP.Homozygosity/PGC.SZ/DATA/AUG2011/*_a12.*

DATS <- c('ab','bon','bulg','carwtc','cat2','dk','dub','edi','mgs2','muc','port','sw1','sw2','top3','ucla','ucl','zhh')

SNP.CUTOFF <- 65  #must be > 65 hmz SNPs in a row to call a ROH - this was the highest power in Dan's paper
KB.CUTOFF <- 500  #all ROHs must be longer than this 
MISS.CUTOFF <- .10 #ROHs cannot have > 10% missingness
SHORT.CUTOFF <- 2200 #Long ROHs are > 1.67Mb (the median)
RARE.CUTOFF <- 1     #Rare ROHs must be unique (median & 1st & 3rd quartile = 1)
OUTLIER.CUTOFF <- 10000 #for removing ROHs longer than 16.6Mb (~expected length of ROH from 2nd cous inbreeding)
HET.CUTOFF <- .05

#These regions come from the paper by Levinson and colleagues
REGIONS.TO.REMOVE <- matrix(c(11,133.65e6,133.69e6,   #GLB1L3/2 del from MGS
                            3, 197.2e6,198.83e6,    #3q29 del from MGS
                            3, 165.61e6, 165.66e6,  #3q26.1 del from MGS
                            1, 144.6e6, 146.3e6,    #1q21.1 del
                            15, 28.7e6, 30.3e6,     #15q13.3 del
                            22, 17.1e6, 20.2e6,     #22q11.2 del
                            16, 29.5e6, 30.1e6,     #16p11.2 del
                            2, 50e6,    51.1e6),     #NRXN1 del
                            ncol=3,byrow=T)

# Run this if you do NOT want to remove regions; otherwise, comment out
#REGIONS.TO.REMOVE <- matrix(c(1,1, 100,2,1, 100),ncol=3,byrow=T)

for (i in 1:length(DATS)){
setwd("~/HD3/SNP.Homozygosity/PGC.SZ/DATA/AUG2011")  #data for "low error" recommendations
ind <- read.table(paste(DATS[i],"_a12.hom.indiv",sep=''),header=TRUE)
hom <-  read.table(paste(DATS[i],"_a12.hom",sep=''),header=TRUE)
over <- read.table(paste(DATS[i],"_a12.hom.overlap",sep=''),header=TRUE)
het <- read.table(paste("../March2011/SZ.PGC/",DATS[i],"8.het",sep=''),header=TRUE)
mds <- read.table(paste("../March2011/SZ.PGC/",DATS[i],"6.mds",sep=''),header=TRUE)
miss <- read.table(paste("../March2011/SZ.PGC/",DATS[i],"2.imiss",sep=''),header=TRUE)

#change the FIDs so that they match up with imputed FIDs
if (DATS[i] %in% c('ab','bulg','dub','edi','port','sw1','sw2','ucl')){
ind$FID <- gsub("_NOXLS","",ind$FID)
hom$FID <- gsub("_NOXLS","",hom$FID)
over$FID <- gsub("_NOXLS","",over$FID)
het$FID <- gsub("_NOXLS","",het$FID)
mds$FID <- gsub("_NOXLS","",mds$FID)
miss$FID <- gsub("_NOXLS","",miss$FID)
}

if (DATS[i]=='carwtc'){
ind$FID <- gsub("*1_case_scz_car_eur_A500K","",ind$FID,fixed=TRUE)
hom$FID <- gsub("*1_case_scz_car_eur_A500K","",hom$FID,fixed=TRUE)
over$FID <- gsub("*1_case_scz_car_eur_A500K","",over$FID,fixed=TRUE)
het$FID <- gsub("*1_case_scz_car_eur_A500K","",het$FID,fixed=TRUE)
mds$FID <- gsub("*1_case_scz_car_eur_A500K","",mds$FID,fixed=TRUE)
miss$FID <- gsub("*1_case_scz_car_eur_A500K","",miss$FID,fixed=TRUE)
ind$FID <- gsub("*0_control_bip_wtc_eur_NOXLS_A500k","",ind$FID,fixed=TRUE)
hom$FID <- gsub("*0_control_bip_wtc_eur_NOXLS_A500k","",hom$FID,fixed=TRUE)
over$FID <- gsub("*0_control_bip_wtc_eur_NOXLS_A500k","",over$FID,fixed=TRUE)
het$FID <- gsub("*0_control_bip_wtc_eur_NOXLS_A500k","",het$FID,fixed=TRUE)
mds$FID <- gsub("*0_control_bip_wtc_eur_NOXLS_A500k","",mds$FID,fixed=TRUE)
miss$FID <- gsub("*0_control_bip_wtc_eur_NOXLS_A500k","",miss$FID,fixed=TRUE)
}

#Assign these files to separate files
assign(paste(DATS[i],".ind",sep=""), ind)
assign(paste(DATS[i],".hom",sep=""), hom)
assign(paste(DATS[i],".overlap",sep=""), over)
assign(paste(DATS[i],".het",sep=""), het)
assign(paste(DATS[i],".mds",sep=""), mds)
assign(paste(DATS[i],".miss",sep=""), miss)


#get final ROH calls that meet cutoffs above; note: NA's in FID of cmm are not a problem; these are ROHs that do not overlap with anyone else
cmm <- find.commonality(over,hom,kb.cutoff=KB.CUTOFF,phet.cutoff=HET.CUTOFF,snp.cutoff=SNP.CUTOFF,miss.cutoff=MISS.CUTOFF)
assign(paste(DATS[i],".comm",sep=""), remove.regions2(cmm,REGIONS.TO.REMOVE,KB.CUTOFF))
assign("comm",eval(parse(text=paste(DATS[i],".comm",sep=""))))
#make person level data from ROHs
assign("roh", make.roh.data2(comm,ind,total.genome.length=2.77e9,rare.com=RARE.CUTOFF,short=SHORT.CUTOFF,suffix=1))

#create above except where no long ROHs
hom.temp <- hom[hom$KB < OUTLIER.CUTOFF,]
over.temp <- over[over$KB < OUTLIER.CUTOFF,]
comm.temp <- find.commonality(over.temp,hom.temp,kb.cutoff=KB.CUTOFF,phet.cutoff=HET.CUTOFF,snp.cutoff=SNP.CUTOFF,miss.cutoff=MISS.CUTOFF)
comm.temp <- remove.regions(comm.temp,REGIONS.TO.REMOVE,KB.CUTOFF)
#make person level data from ROHs
assign("roh2", make.roh.data2(comm.temp,ind,total.genome.length=2.77e9,rare.com=RARE.CUTOFF,short=SHORT.CUTOFF,suffix=1))

#merge minor allele load, MDS, and missingness info
roh <- merge(roh,mds,by=c("FID","IID"),all=TRUE) #b
roh <- merge(roh,miss,by=c("FID","IID"),all=TRUE) #c
roh <- merge(roh,het,by=c("FID","IID"),all=TRUE) #c

roh2 <- merge(roh2,mds,by=c("FID","IID"),all=TRUE) #b
roh2 <- merge(roh2,miss,by=c("FID","IID"),all=TRUE) #c
roh2 <- merge(roh2,het,by=c("FID","IID"),all=TRUE) #c

#control for MDS WITHIN each dataset info:
roh$tot2.1a <- lm(roh$tot2.1 ~ roh$C1 + roh$C2 + roh$C3 + roh$C4 + roh$C5)$residuals
roh2$tot2.1a <- lm(roh2$tot2.1 ~ roh2$C1 + roh2$C2 + roh2$C3 + roh2$C4 + roh2$C5)$residuals

assign(paste(DATS[i],".roh",sep=""), roh)
assign(paste(DATS[i],".roh2",sep=""), roh2)

}

setwd("~/HD3/SNP.Homozygosity/PGC.SZ/Results.April2011")
################################################






################################################
#3) Put all the ROH-level datasets together

#this is for the post-processed ROH level datasets
COMM <- rbind(ab.comm,bon.comm,bulg.comm,carwtc.comm,cat2.comm,dk.comm,dub.comm,edi.comm,mgs2.comm,muc.comm,port.comm,sw1.comm,sw2.comm,top3.comm,ucla.comm,ucl.comm,zhh.comm)

#this is for the pre-processed ROH level datasets
HOM <- rbind(ab.hom,bon.hom,bulg.hom,carwtc.hom,cat2.hom,dk.hom,dub.hom,edi.hom,mgs2.hom,muc.hom,port.hom,sw1.hom,sw2.hom,top3.hom,ucla.hom,ucl.hom,zhh.hom)
HOM$PMISS <- round(1 - HOM$PHOM - HOM$PHET,3)
HOM$NMISS <- HOM$PMISS*HOM$NSNP
fff <- strsplit(as.character(HOM$FID),"_",fixed=TRUE)
HOM$case <- as.factor(unlist(lapply(fff,function(x) x[1])))
HOM$dataset <- as.factor(unlist(lapply(fff,function(x) x[3])))
HOM$dataset.case <- paste(HOM$dataset,HOM$case,sep=".")
HOM$platform <- NA
HOM$platform[HOM$dataset %in% c('ab','carwtc','cat2','port','sw1','ucl','zhh')] <- 'A5'
HOM$platform[HOM$dataset %in% c('bulg','dub','edi','mgs2','sw2','top3')] <- 'A6'
HOM$platform[HOM$dataset %in% c('bon','ucla')] <- 'I550'
HOM$platform[HOM$dataset %in% c('dk')] <- 'I650'
HOM$platform[HOM$dataset %in% c('muc')] <- 'I317'
HOM$platform <- as.factor(HOM$platform)
HOM2 <- HOM[HOM$PMISS <= MISS.CUTOFF & HOM$KB >= KB.CUTOFF & HOM$NSNP >= SNP.CUTOFF,]

#ROHs lost for each reason:
HOM.lost.kb <- HOM[HOM$KB >= KB.CUTOFF,]
HOM.lost.snp <- HOM.lost.kb[HOM.lost.kb$NSNP >= SNP.CUTOFF,]
HOM.lost.miss <- HOM.lost.snp[HOM.lost.snp$PMISS < MISS.CUTOFF,]
nrow(HOM)- nrow(HOM.lost.kb); (nrow(HOM)- nrow(HOM.lost.kb))/nrow(HOM)  #number lost due to KB CUTOFF
nrow(HOM.lost.kb) - nrow(HOM.lost.snp); (nrow(HOM.lost.kb) - nrow(HOM.lost.snp))/nrow(HOM)  #number lost due to SNP CUTOF
nrow(HOM.lost.snp) - nrow(HOM.lost.miss); (nrow(HOM.lost.snp) - nrow(HOM.lost.miss))/nrow(HOM) #number lost due to MISS CUTOFF
################################################






################################################
#4) Put all the person-level datasets together
DAT <- rbind(ab.roh,bon.roh,bulg.roh,carwtc.roh,cat2.roh,dk.roh,dub.roh,edi.roh,mgs2.roh,muc.roh,port.roh,sw1.roh,sw2.roh,top3.roh,ucla.roh,ucl.roh,zhh.roh)
fff <- strsplit(as.character(DAT$FID),"_",fixed=TRUE)
DAT$case <- as.factor(unlist(lapply(fff,function(x) x[1])))
DAT$dataset <- as.factor(unlist(lapply(fff,function(x) x[3])))
DAT$dataset.case <- paste(DAT$dataset,DAT$case,sep=".")

#DAT2 has the outlier runs (>16.7Mb) removed
DAT2 <- rbind(ab.roh2,bon.roh2,bulg.roh2,carwtc.roh2,cat2.roh2,dk.roh2,dub.roh2,edi.roh2,mgs2.roh2,muc.roh2,port.roh2,sw1.roh2,sw2.roh2,top3.roh2,ucla.roh2,ucl.roh2,zhh.roh2)
DAT2 <- DAT2[,c(1,3,5,10,13,16,19,20)]

#both of these together = TOT.med; the ".out" suffix variables are ones with outlier runs removed
TOT.med <- merge(DAT,DAT2,by="UID",all=TRUE,suffix=c('','.out'))
TOT.med$SZ <- (TOT.med$case=='case')*1

#make the measures of %ROH & Fh
TOT.med$pct <- TOT.med$tot2.1*100
TOT.med$pct.rare <- (TOT.med$proh.short.r2.1 + TOT.med$proh.big.r2.1)*100
TOT.med$pct.comm <- (TOT.med$proh.short.c2.1 + TOT.med$proh.big.c2.1)*100
TOT.med$pct.out <- TOT.med$tot2.1.out*100
TOT.med$pct.rare.out <- (TOT.med$proh.short.r2.1.out + TOT.med$proh.big.r2.1.out)*100
TOT.med$pct.comm.out <- (TOT.med$proh.short.c2.1.out + TOT.med$proh.big.c2.1.out)*100
TOT.med$Fh <- TOT.med$F*100

#Add the platform information
#NOTE: all subjects in port dataset have 4 grandparents from the Azores or Madeira islands; there is no case-control confounding of this
TOT.med$platform <- 'A500'
TOT.med$platform[TOT.med$dataset %in% c('ab','port','sw1','ucl')] <- 'A5'
TOT.med$platform[TOT.med$dataset %in% c('bulg','dub','edi','mgs2','sw2','top3')] <- 'A6'
#TOT.med$platform[TOT.med$dataset %in% c('bon','ucla')] <- 'I550'
TOT.med$platform[TOT.med$dataset %in% c('dk')] <- 'I650'
TOT.med$platform[TOT.med$dataset %in% c('muc')] <- 'I317'

#Grab the first 20 PCs & merge
#NOTE: PC1 - PC20 are PCs ACROSS the datasets (done using common SNPs) and can be used ACROSS datasets
pcs <- read.table("~/HD3/SNP.Homozygosity/PGC.SZ/DATA/roh.TEST2/SCZ17b.mds_cov_cp",header=TRUE)
names(pcs)[4:23] <- paste("PC",1:20,sep="")
pcs$FID <- gsub("*1_case_scz_car_eur_A500K","",pcs$FID,fixed=TRUE)
pcs$FID <- gsub("*0_control_bip_wtc_eur_NOXLS_A500k","",pcs$FID,fixed=TRUE)
pcs$FID <- gsub("_NOXLS","",pcs$FID)
pcs$UID <- paste(pcs$FID,pcs$IID,sep='')
TOT.med2 <- merge(TOT.med,pcs,by='UID',all=TRUE)

#look at distn's of Fh & ROH & F_MISS
quantile(TOT.med2$Fh,probs=seq(0,1,.01))
quantile(TOT.med2$pct,probs=seq(0,1,.01))
quantile(TOT.med2$F_MISS,probs=seq(0,1,.01))

#Make new dataset removing those missing in the PC dimensions; this gives us the subjects used in previous PGC publications
#MAIN <- TOT.med2[! is.na(TOT.med2$st16),]
MAIN <- TOT.med2

#st16 variances between outlier & no outlier data
sd(MAIN$pct[MAIN$pct<6.25]);sd(MAIN$pct) # 71% of the SD; should increase SE(beta) by ~1.4
sd(MAIN$pct[MAIN$pct<3.125]);sd(MAIN$pct) # 1/2 of the SD; should increase SE(beta) by ~2
################################################

 



################################################
#5) Look at descriptive stats across datasets
DATS <- c('ab','bon','bulg','carwtc','cat2','dk','dub','edi','mgs2','muc','port','sw1','sw2','top3','ucla','ucl','zhh')
clnms <- c('n','nSNP','mean.%','sd.%','max.%','max.snp','mean.snp','min.len','max.len','mean.len','min.F','max.F','min.Fz','max.Fz','sd.F','r(%,Fh)','r(%,miss)','r(SZ,miss)','r(%,SZ)','r(Fh,miss)','r(%,Fh)')
RESULTS <- matrix(nrow=length(DATS),ncol=length(clnms),dimnames=list(DATS,clnms))

kk <- 0
for (REP in DATS){
kk <- kk+1
mypct <- MAIN$pct[MAIN$dataset==REP]
max.snp <- max(HOM2$NSNP[HOM2$dataset==REP])
mean.snp <- mean(HOM2$NSNP[HOM2$dataset==REP])
min.len <- min(HOM2$KB[HOM2$dataset==REP])
max.len <- max(HOM2$KB[HOM2$dataset==REP])
mean.len <- mean(HOM2$KB[HOM2$dataset==REP])
MyF <- MAIN$Fh[MAIN$dataset==REP]
MyMiss <- MAIN$F_MISS[MAIN$dataset==REP]
sz <- MAIN$SZ[MAIN$dataset==REP]
min.F <- min(MyF)
max.F <- max(MyF)
min.Fz <- min(scale(MyF))
max.Fz <- max(scale(MyF))
sd.F <- sd(MyF)
nsnp <- as.numeric(strsplit(system(paste("wc -l ~/HD3/SNP.Homozygosity/PGC.SZ/DATA/March2011/SZ.PGC/",REP,"4.bim",sep=''),intern=TRUE),' ')[[1]][1])
nnn <- nrow(MAIN[MAIN$dataset==REP,])
RESULTS[kk,] <- c(nnn,nsnp,round(mean(mypct),2),round(sd(mypct),2),round(max(mypct),1),round(max.snp,1),round(mean.snp,1),round(min.len,1),round(max.len,1),round(mean.len,1),round(min.F,2),round(max.F,2),round(min.Fz,2),round(max.Fz,2),round(sd.F,2),round(cor(mypct,MyF),2),round(cor(mypct,MyMiss,use='pairwise.complete.obs'),2),round(cor(sz,MyF),2),round(cor(sz,mypct),2),round(cor(MyF,MyMiss),2),round(cor(MyF,mypct),2))
}

RESULTS <- as.data.frame(RESULTS)
RESULTS$platform <- 'A500'
RESULTS$platform[rownames(RESULTS) %in% c('ab','port','sw1','ucl')] <- 'A5'
RESULTS$platform[rownames(RESULTS) %in% c('bulg','dub','edi','mgs2','sw2','top3')] <- 'A6'
RESULTS$platform[rownames(RESULTS) %in% c('bon','ucla')] <- 'I550'
RESULTS$platform[rownames(RESULTS) %in% c('dk')] <- 'I650'
RESULTS$platform[rownames(RESULTS) %in% c('muc')] <- 'I317'
RESULTS
################################################











## @@@@@@@@@@@@   ******   @@@@@@@@@@@@   ******   @@@@@@@@@@@@   ******   @@@@@@@@@@@@










################################################
#6) Add in the imputed data
#Red.Imputed1a,b,c,d = high QC (.99 r2, .90 r2 individual datasets)
# *4 data is VIF 2, *7 is VIF 10;   it makes no difference which one is used, both are sig
# *9 is VIF 2, 35 SNPs, 0 het
# *10 is VIF 10, 60 SNPs, 0 het

DATS <- c('Red.Imputed4a','Red.Imputed4b','Red.Imputed4c','Red.Imputed4d')  

SNP.CUTOFF <- 65  #must be > 65 hmz SNPs in a row to call a ROH - this was the highest power in Dan's paper
KB.CUTOFF <- 500  #all ROHs must be longer than this 
MISS.CUTOFF <- .10 #ROHs cannot have > 10% missingness (basically ignored)
SHORT.CUTOFF <- 2000 #Long ROHs are > 2 Mb (the median)
RARE.CUTOFF <- 6     #Rare ROHs must be shared < 5 times in sample of ~8000 ROHs (8000 ROHs typically observed in sample of 5570) 
OUTLIER.CUTOFF <- 10000 #for suppl. analysis removing ROHs longer than 10Mb (~expected length of ROH from inbreeding from 5 generations or more back)
HET.CUTOFF <- .05    #ignored

#These regions come from the paper by Levinson and colleagues
REGIONS.TO.REMOVE <- matrix(c(11,133.65e6,133.69e6,   #GLB1L3/2 del from MGS
                            3, 197.2e6,198.83e6,    #3q29 del from MGS
                            3, 165.61e6, 165.66e6,  #3q26.1 del from MGS
                            1, 144.6e6, 146.3e6,    #1q21.1 del
                            15, 28.7e6, 30.3e6,     #15q13.3 del
                            22, 17.1e6, 20.2e6,     #22q11.2 del
                            16, 29.5e6, 30.1e6,     #16p11.2 del
                            2, 50e6,    51.1e6),     #NRXN1 del
                            ncol=3,byrow=T)

# Run this if you do NOT want to remove regions; otherwise, comment out
#REGIONS.TO.REMOVE <- matrix(c(1,1, 100,2,1, 100),ncol=3,byrow=T)

for (i in 1:length(DATS)){
setwd("~/HD3/SNP.Homozygosity/PGC.SZ/DATA/AUG2011.Imputed")  #data from imputed SNPs
assign(paste(DATS[i],".ind",sep=""), read.table(paste(DATS[i],"10.hom.indiv",sep=''),header=TRUE))
assign(paste(DATS[i],".hom",sep=""), read.table(paste(DATS[i],"10.hom",sep=''),header=TRUE))
assign(paste(DATS[i],".overlap",sep=""), read.table(paste(DATS[i],"10.hom.overlap",sep=''),header=TRUE))
assign(paste(DATS[i],".het",sep=""), read.table(paste(DATS[i],"8.het",sep=''),header=TRUE))
assign(paste(DATS[i],".miss",sep=""), read.table(paste(DATS[i],"2.imiss",sep=''),header=TRUE))

assign("over",eval(parse(text=paste(DATS[i],".overlap",sep=""))))
assign("hom",eval(parse(text=paste(DATS[i],".hom",sep=""))))
assign("miss",eval(parse(text=paste(DATS[i],".miss",sep=""))))
assign("ind",eval(parse(text=paste(DATS[i],".ind",sep=""))))
assign("het",eval(parse(text=paste(DATS[i],".het",sep=""))))

#get final ROH calls that meet cutoffs above; note: NA's in FID of cmm are not a problem; these are ROHs that do not overlap with anyone else
cmm <- find.commonality(over,hom,kb.cutoff=KB.CUTOFF,phet.cutoff=HET.CUTOFF,snp.cutoff=SNP.CUTOFF,miss.cutoff=MISS.CUTOFF)
assign(paste(DATS[i],".comm",sep=""), remove.regions2(cmm,REGIONS.TO.REMOVE,KB.CUTOFF))
assign("comm",eval(parse(text=paste(DATS[i],".comm",sep=""))))
#make person level data from ROHs
assign("roh", make.roh.data2(comm,ind,total.genome.length=2.77e9,rare.com=RARE.CUTOFF,short=SHORT.CUTOFF,suffix=1))

#create above except where no long ROHs
hom.temp <- hom[hom$KB < OUTLIER.CUTOFF,]
over.temp <- over[over$KB < OUTLIER.CUTOFF,]
comm.temp <- find.commonality(over.temp,hom.temp,kb.cutoff=KB.CUTOFF,phet.cutoff=HET.CUTOFF,snp.cutoff=SNP.CUTOFF,miss.cutoff=MISS.CUTOFF)
comm.temp <- remove.regions(comm.temp,REGIONS.TO.REMOVE,KB.CUTOFF)
#make person level data from ROHs
assign("roh2", make.roh.data2(comm.temp,ind,total.genome.length=2.77e9,rare.com=RARE.CUTOFF,short=SHORT.CUTOFF,suffix=1))

#merge missingness and Fh information
roh <- merge(roh,miss,by=c("FID","IID"),all=TRUE) #b
roh <- merge(roh,het,by=c("FID","IID"),all=TRUE) #c
roh2 <- merge(roh2,miss,by=c("FID","IID"),all=TRUE) #b
roh2 <- merge(roh2,het,by=c("FID","IID"),all=TRUE) #c

assign(paste(DATS[i],".roh",sep=""), roh)
assign(paste(DATS[i],".roh2",sep=""), roh2)

}
################################################







################################################
#7) Put all the ROH-level datasets together

#this is for the post-processed ROH level datasets
COMM_I <- rbind(Red.Imputed1a.comm,Red.Imputed1b.comm,Red.Imputed1c.comm,Red.Imputed1d.comm)   #CHANGE
#COMM_I <- rbind(Full.Imputed.a.comm,Full.Imputed.b.comm,Full.Imputed.c.comm,Full.Imputed.d.comm)   #CHANGE

#this is for the pre-processed ROH level datasets
HOM_I <- rbind(Red.Imputed1a.hom,Red.Imputed1b.hom,Red.Imputed1c.hom,Red.Imputed1d.hom)  #CHANGE
#HOM_I <- rbind(Full.Imputed.a.hom,Full.Imputed.b.hom,Full.Imputed.c.hom,Full.Imputed.d.hom)  #CHANGE
HOM_I$PMISS <- round(1 - HOM_I$PHOM - HOM_I$PHET,3)
HOM_I$NMISS <- HOM_I$PMISS*HOM_I$NSNP
fff <- strsplit(as.character(HOM_I$FID),"_",fixed=TRUE)
HOM_I$case <- as.factor(unlist(lapply(fff,function(x) x[1])))
HOM_I$dataset <- as.factor(unlist(lapply(fff,function(x) x[3])))
HOM_I$dataset.case <- paste(HOM_I$dataset,HOM_I$case,sep=".")
HOM_I$platform <- NA
HOM_I$platform[HOM_I$dataset %in% c('ab','carwtc','cat2','port','sw1','ucl','zhh')] <- 'A5'
HOM_I$platform[HOM_I$dataset %in% c('bulg','dub','edi','mgs2','sw2','top3')] <- 'A6'
HOM_I$platform[HOM_I$dataset %in% c('bon','ucla')] <- 'I550'
HOM_I$platform[HOM_I$dataset %in% c('dk')] <- 'I650'
HOM_I$platform[HOM_I$dataset %in% c('muc')] <- 'I317'
HOM_I$platform <- as.factor(HOM_I$platform)
HOM_I2 <- HOM_I[HOM_I$KB >= KB.CUTOFF & HOM_I$NSNP >= SNP.CUTOFF,]

#ROHs lost for each reason:
HOM_I.lost.kb <- HOM_I[HOM_I$KB >= KB.CUTOFF,]
HOM_I.lost.snp <- HOM_I.lost.kb[HOM_I.lost.kb$NSNP >= SNP.CUTOFF,]
nrow(HOM_I)- nrow(HOM_I.lost.kb); (nrow(HOM_I)- nrow(HOM_I.lost.kb))/nrow(HOM_I)  #number lost due to KB CUTOFF
nrow(HOM_I.lost.kb) - nrow(HOM_I.lost.snp); (nrow(HOM_I.lost.kb) - nrow(HOM_I.lost.snp))/nrow(HOM_I)  #number lost due to SNP CUTOF
################################################







################################################
#8) Put all the person-level datasets together
DAT_I <- rbind(Red.Imputed1a.roh,Red.Imputed1b.roh,Red.Imputed1c.roh,Red.Imputed1d.roh)  #CHANGE
#DAT_I <- rbind(Full.Imputed.a.roh,Full.Imputed.b.roh,Full.Imputed.c.roh,Full.Imputed.d.roh)  #CHANGE

fff <- strsplit(as.character(DAT_I$FID),"_",fixed=TRUE)
DAT_I$case <- as.factor(unlist(lapply(fff,function(x) x[1])))
DAT_I$dataset <- as.factor(unlist(lapply(fff,function(x) x[3])))
DAT_I$dataset.case <- paste(DAT_I$dataset,DAT_I$case,sep=".")

#DAT_I2 has the outlier runs (>16.7Mb) removed
DAT_I2 <- rbind(Red.Imputed1a.roh2,Red.Imputed1b.roh2,Red.Imputed1c.roh2,Red.Imputed1d.roh2)
DAT_I2 <- DAT_I2[,c(1,3,5,10,13,16,19,20)]

#both of these together = TOT.med_I; the ".out" suffix variables are ones with outlier runs removed
TOT.med_I <- merge(DAT_I,DAT_I2,by="UID",all=TRUE,suffix=c('','.out'))
TOT.med_I$SZ <- (TOT.med_I$case=='case')*1

#make the measures of %ROH & Fh
TOT.med_I$pct <- TOT.med_I$tot2.1*100
TOT.med_I$pct.rare <- (TOT.med_I$proh.short.r2.1 + TOT.med_I$proh.big.r2.1)*100
TOT.med_I$pct.comm <- (TOT.med_I$proh.short.c2.1 + TOT.med_I$proh.big.c2.1)*100
TOT.med_I$pct.out <- TOT.med_I$tot2.1.out*100
TOT.med_I$pct.rare.out <- (TOT.med_I$proh.short.r2.1.out + TOT.med_I$proh.big.r2.1.out)*100
TOT.med_I$pct.comm.out <- (TOT.med_I$proh.short.c2.1.out + TOT.med_I$proh.big.c2.1.out)*100
TOT.med_I$Fh <- TOT.med_I$F*100

#look at distn's of Fh & ROH & F_MISS
quantile(TOT.med_I$Fh,probs=seq(0,1,.01))
quantile(TOT.med_I$pct,probs=seq(0,1,.01))
quantile(TOT.med_I$F_MISS,probs=seq(0,1,.01))

#Make new dataset removing those missing in the PC dimensions; this gives us the subjects used in previous PGC publications
MAIN_I <- TOT.med_I

#st16 variances between outlier & no outlier data
sd(MAIN_I$pct[MAIN_I$pct<6.25]);sd(MAIN_I$pct) # 71% of the SD; should increase SE(beta) by ~1.4
sd(MAIN_I$pct[MAIN_I$pct<3.125]);sd(MAIN_I$pct) # 1/2 of the SD; should increase SE(beta) by ~2
################################################







################################################
#9) Look at descriptive stats across datasets
DATS <- c('ab','bon','bulg','carwtc','cat2','dk','dub','edi','mgs2','muc','port','sw1','sw2','top3','ucla','ucl','zhh')
clnms <- c('n','nSNP','mean.%','sd.%','max.%','max.snp','mean.snp','min.len','max.len','mean.len','min.F','max.F','min.Fz','max.Fz','sd.F','r(%,Fh)','r(%,miss)','r(SZ,miss)','r(%,SZ)','r(Fh,miss)','r(%,Fh)')
RESULTS_I <- matrix(nrow=length(DATS),ncol=length(clnms),dimnames=list(DATS,clnms))

kk <- 0
for (REP in DATS){
kk <- kk+1
mypct <- MAIN_I$pct[MAIN_I$dataset==REP]
max.snp <- max(HOM_I2$NSNP[HOM_I2$dataset==REP])
mean.snp <- mean(HOM_I2$NSNP[HOM_I2$dataset==REP])
min.len <- min(HOM_I2$KB[HOM_I2$dataset==REP])
max.len <- max(HOM_I2$KB[HOM_I2$dataset==REP])
mean.len <- mean(HOM_I2$KB[HOM_I2$dataset==REP])
MyF <- MAIN_I$Fh[MAIN_I$dataset==REP]
MyMiss <- MAIN_I$F_MISS[MAIN_I$dataset==REP]
sz <- MAIN$SZ[MAIN$dataset==REP]
min.F <- min(MyF)
max.F <- max(MyF)
min.Fz <- min(scale(MyF))
max.Fz <- max(scale(MyF))
sd.F <- sd(MyF)
nsnp <- as.numeric(strsplit(system(paste("wc -l ~/HD3/SNP.Homozygosity/PGC.SZ/DATA/March2011/SZ.PGC/",REP,"4.bim",sep=''),intern=TRUE),' ')[[1]][1])
nnn <- nrow(MAIN_I[MAIN_I$dataset==REP,])
RESULTS_I[kk,] <- c(nnn,nsnp,round(mean(mypct),2),round(sd(mypct),2),round(max(mypct),1),round(max.snp,1),round(mean.snp,1),round(min.len,1),round(max.len,1),round(mean.len,1),round(min.F,2),round(max.F,2),round(min.Fz,2),round(max.Fz,2),round(sd.F,2),round(cor(mypct,MyF),2),round(cor(mypct,MyMiss,use='pairwise.complete.obs'),2),round(cor(sz,MyF),2),round(cor(sz,mypct),2),round(cor(MyF,MyMiss),2),round(cor(MyF,mypct),2))
}

RESULTS_I <- as.data.frame(RESULTS_I)
RESULTS_I$platform <- 'A500'
RESULTS_I$platform[rownames(RESULTS_I) %in% c('ab','port','sw1','ucl')] <- 'A5'
RESULTS_I$platform[rownames(RESULTS_I) %in% c('bulg','dub','edi','mgs2','sw2','top3')] <- 'A6'
RESULTS_I$platform[rownames(RESULTS_I) %in% c('bon','ucla')] <- 'I550'
RESULTS_I$platform[rownames(RESULTS_I) %in% c('dk')] <- 'I650'
RESULTS_I$platform[rownames(RESULTS_I) %in% c('muc')] <- 'I317'
RESULTS_I
################################################








################################################
#10) Merge the Imputed & raw datasets together

#rename MAIN_I
names(MAIN_I) <- paste(names(MAIN_I),".imp",sep='')

#Merge
MAIN2 <- merge(MAIN,MAIN_I,by.x="UID",by.y="UID.imp",all=TRUE)
################################################








################################################
#11) Preliminary results - IMPUTED DATA

#Overall result
summary(md1 <- glmer(SZ ~ tot2.1.imp + (1|dataset),data=MAIN2,family='binomial'))

#Result controlling for 20 PCs and 2 data quality metrics (Fh and missingness in the raw data); this is the base model going forward
summary(md2 <- glmer(SZ ~ tot2.1.imp+F_MISS+Fh.imp+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+(1|dataset),data=MAIN2,family='binomial'))

#Some alternative regression models just to ensure that the above results are stable across different ways of controlling for confounds
  #using Fh from original (raw) data
summary(glmer(SZ ~ tot2.1.imp + F_MISS+Fh+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+ (1|dataset),data=MAIN2,family='binomial'))
  #not controlling for snp quality metrics
summary(glmer(SZ ~ tot2.1.imp+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+ (1|dataset),data=MAIN2,family='binomial'))
  #not controlling for Fh
summary(glmer(SZ ~ tot2.1.imp+F_MISS+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+ (1|dataset),data=MAIN2,family='binomial'))
  #not controlling for missingness
summary(glmer(SZ ~ tot2.1.imp+Fh.imp+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+ (1|dataset),data=MAIN2,family='binomial'))
   #controlling for 10 PCs instead of 20
summary(glmer(SZ ~ tot2.1.imp+F_MISS+Fh.imp+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+(1|dataset),data=MAIN2,family='binomial'))
   #all significant covariate-by-dataset interactions; this model cannot be fit using glmer (takes hours); fixed effect model should be similar
summary(glm(SZ ~ tot2.1.imp +F_MISS*dataset+Fh.imp+PC1*dataset+PC2*dataset+PC3*dataset+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=MAIN2,family='binomial'))
   #to see which interactions are significant (justification for the above analysis)
summary(red <- glm(SZ ~ tot2.1.imp +dataset + F_MISS+Fh.imp+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20,data=MAIN2,family='binomial'))
summary(full1 <- glm(SZ ~ tot2.1.imp +dataset + F_MISS*dataset+Fh.imp+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20,data=MAIN2,family='binomial'))
summary(full2 <- glm(SZ ~ tot2.1.imp +dataset + F_MISS+Fh.imp*dataset+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20,data=MAIN2,family='binomial'))
summary(full3 <- glm(SZ ~ tot2.1.imp +dataset + F_MISS+Fh.imp+PC1*dataset+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20,data=MAIN2,family='binomial'))
summary(full4 <- glm(SZ ~ tot2.1.imp +dataset + F_MISS+Fh.imp+PC1+PC2*dataset+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20,data=MAIN2,family='binomial'))
summary(full5 <- glm(SZ ~ tot2.1.imp +dataset + F_MISS+Fh.imp+PC1+PC2+PC3*dataset+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20,data=MAIN2,family='binomial'))
summary(full6 <- glm(SZ ~ tot2.1.imp +dataset + F_MISS+Fh.imp+PC1+PC2+PC3+PC4+PC5+PC6*dataset+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20,data=MAIN2,family='binomial'))
summary(full7 <- glm(SZ ~ tot2.1.imp +dataset + F_MISS+Fh.imp+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9*dataset+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20,data=MAIN2,family='binomial'))
summary(full8 <- glm(SZ ~ tot2.1.imp +dataset + F_MISS+Fh.imp+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14*dataset+PC15+PC16+PC17+PC18+PC19+PC20,data=MAIN2,family='binomial'))
anova(red,full1,test='Chisq') #miss by dataset is sig
anova(red,full2,test='Chisq') #Fh by dataset is not sig
anova(red,full3,test='Chisq') #PC1 by dataset is sig
anova(red,full4,test='Chisq') #PC2 by dataset is sig
anova(red,full5,test='Chisq') #PC3 by dataset is sig
anova(red,full6,test='Chisq') #PC6 by dataset is not sig
anova(red,full7,test='Chisq') #PC9 by dataset is not sig
anova(red,full8,test='Chisq') #PC14 by dataset is not sig

#removing outliers on Froh to see if results hold; they appear to, though become non-sig at low thresholds due to little var in pred
summary(glm(SZ ~ tot2.1.imp+F_MISS+Fh.imp+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+dataset,data=MAIN2[MAIN2$tot2.1.imp < .125+mean(MAIN2$tot2.1.imp),],family='binomial'))  #avuncular-neice/nephew or half-sib level
summary(glm(SZ ~ tot2.1.imp+F_MISS+Fh.imp+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+dataset,data=MAIN2[MAIN2$tot2.1.imp < .0625+mean(MAIN2$tot2.1.imp),],family='binomial'))  #cousin-cousin level
summary(glm(SZ ~ tot2.1.imp+F_MISS+Fh.imp+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+dataset,data=MAIN2[MAIN2$tot2.1.imp < .03125+mean(MAIN2$tot2.1.imp),],family='binomial'))   # half cousin level
summary(glm(SZ ~ tot2.1.imp+F_MISS+Fh.imp+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+dataset,data=MAIN2[MAIN2$tot2.1.imp < .005+mean(MAIN2$tot2.1.imp),],family='binomial'))  # ancient inbreeding level

#removing long (> 10Mb) ROHs (used fixed effects model bc the rand effects won't converge)
summary(glm(SZ ~ tot2.1.out.imp+F_MISS+Fh.imp+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+dataset,data=MAIN2,family='binomial'))

#turn the problem around so that we're predicting Froh rather than SZ - just to make sure the glmer algorithm is working appropriately (it does, giving ~ identical answers)
summary(lmer(tot2.1.imp ~ SZ +(1|dataset),data=MAIN2))
summary(lmer(tot2.1.imp ~ SZ +F_MISS+Fh.imp+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+(1|dataset),data=MAIN2))

#Do common or rare ROHs drive this effect? They both appear to
summary(md2.rare <- glm(SZ ~ pct.rare.imp+F_MISS+Fh.imp+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+dataset,data=MAIN2,family='binomial'))
summary(md2.comm <- glm(SZ ~ pct.comm.imp+F_MISS+Fh.imp+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+dataset,data=MAIN2,family='binomial'))

#Do short or long ROHs drive this effect? They both appear to, although short even more than long, suggesting that it is ancient inbreeding driving the effect
summary(md2.rare <- glm(SZ ~ I(proh.short.c2.1.out.imp + proh.short.r2.1.out.imp)+F_MISS+Fh.imp+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+dataset,data=MAIN2,family='binomial'))
summary(md2.comm <- glm(SZ ~ I(proh.big.c2.1.out.imp + proh.big.r2.1.out.imp)+F_MISS+Fh.imp+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+dataset,data=MAIN2,family='binomial'))
###############################################








################################################
#12) Preliminary results - RAW DATA

#Overall result
summary(md1 <- glmer(SZ ~ tot2.1 + (1|dataset),data=MAIN2,family='binomial'))

#Result controlling for 20 PCs and 2 data quality metrics (Fh and missingness in the raw data); this is the base model going forward
summary(md2 <- glmer(SZ ~ tot2.1+F_MISS+Fh+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+(1|dataset),data=MAIN2,family='binomial'))

#Some alternative regression models just to ensure that the above results are stable across different ways of controlling for confounds
  #using Fh from original (raw) data
summary(glmer(SZ ~ tot2.1 + F_MISS+Fh+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+ (1|dataset),data=MAIN2,family='binomial'))
  #not controlling for snp quality metrics
summary(glmer(SZ ~ tot2.1+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+ (1|dataset),data=MAIN2,family='binomial'))
  #not controlling for Fh
summary(glmer(SZ ~ tot2.1+F_MISS+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+ (1|dataset),data=MAIN2,family='binomial'))
  #not controlling for missingness
summary(glmer(SZ ~ tot2.1+Fh+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+ (1|dataset),data=MAIN2,family='binomial'))
   #controlling for 10 PCs instead of 20
summary(glmer(SZ ~ tot2.1+F_MISS+Fh+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+(1|dataset),data=MAIN2,family='binomial'))
   #all significant covariate-by-dataset interactions; this model cannot be fit using glmer (takes hours); fixed effect model should be similar
summary(glm(SZ ~ tot2.1 +F_MISS*dataset+Fh+PC1*dataset+PC2*dataset+PC3*dataset+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=MAIN2,family='binomial'))
   #to see which interactions are significant (justification for the above analysis)
summary(red <- glm(SZ ~ tot2.1 +dataset + F_MISS+Fh+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20,data=MAIN2,family='binomial'))
summary(full1 <- glm(SZ ~ tot2.1 +dataset + F_MISS*dataset+Fh+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20,data=MAIN2,family='binomial'))
summary(full2 <- glm(SZ ~ tot2.1 +dataset + F_MISS+Fh*dataset+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20,data=MAIN2,family='binomial'))
summary(full3 <- glm(SZ ~ tot2.1 +dataset + F_MISS+Fh+PC1*dataset+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20,data=MAIN2,family='binomial'))
summary(full4 <- glm(SZ ~ tot2.1 +dataset + F_MISS+Fh+PC1+PC2*dataset+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20,data=MAIN2,family='binomial'))
summary(full5 <- glm(SZ ~ tot2.1 +dataset + F_MISS+Fh+PC1+PC2+PC3*dataset+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20,data=MAIN2,family='binomial'))
summary(full6 <- glm(SZ ~ tot2.1 +dataset + F_MISS+Fh+PC1+PC2+PC3+PC4+PC5+PC6*dataset+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20,data=MAIN2,family='binomial'))
summary(full7 <- glm(SZ ~ tot2.1 +dataset + F_MISS+Fh+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9*dataset+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20,data=MAIN2,family='binomial'))
summary(full8 <- glm(SZ ~ tot2.1 +dataset + F_MISS+Fh+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14*dataset+PC15+PC16+PC17+PC18+PC19+PC20,data=MAIN2,family='binomial'))
anova(red,full1,test='Chisq') #miss by dataset is sig
anova(red,full2,test='Chisq') #Fh by dataset is not sig
anova(red,full3,test='Chisq') #PC1 by dataset is sig
anova(red,full4,test='Chisq') #PC2 by dataset is sig
anova(red,full5,test='Chisq') #PC3 by dataset is sig
anova(red,full6,test='Chisq') #PC6 by dataset is not sig
anova(red,full7,test='Chisq') #PC9 by dataset is not sig
anova(red,full8,test='Chisq') #PC14 by dataset is not sig

#removing outliers on Froh to see if results hold; they appear to, though become non-sig at low thresholds due to little var in pred
summary(glm(SZ ~ tot2.1+F_MISS+Fh.imp+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+dataset,data=MAIN2[MAIN2$tot2.1.imp < .125+mean(MAIN2$tot2.1),],family='binomial'))  #avuncular-neice/nephew or half-sib level
summary(glm(SZ ~ tot2.1+F_MISS+Fh.imp+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+dataset,data=MAIN2[MAIN2$tot2.1.imp < .0625+mean(MAIN2$tot2.1),],family='binomial'))  #cousin-cousin level
summary(glm(SZ ~ tot2.1+F_MISS+Fh.imp+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+dataset,data=MAIN2[MAIN2$tot2.1.imp < .03125+mean(MAIN2$tot2.1),],family='binomial'))   # half cousin level
summary(glm(SZ ~ tot2.1+F_MISS+Fh.imp+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+dataset,data=MAIN2[MAIN2$tot2.1.imp < .005+mean(MAIN2$tot2.1),],family='binomial'))  # ancient inbreeding level

#removing long (> 10Mb) ROHs (used fixed effects model bc the rand effects won't converge)
summary(glm(SZ ~ tot2.1.out+F_MISS+Fh+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+dataset,data=MAIN2,family='binomial'))

#turn the problem around so that we're predicting Froh rather than SZ - just to make sure the glmer algorithm is working appropriately (it does, giving ~ identical answers)
summary(lmer(tot2.1 ~ SZ +(1|dataset),data=MAIN2))
summary(lmer(tot2.1 ~ SZ +F_MISS+Fh+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+(1|dataset),data=MAIN2))

#Do common or rare ROHs drive this effect? They both appear to
summary(md2.rare <- glm(SZ ~ pct.rare+F_MISS+Fh+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+dataset,data=MAIN2,family='binomial'))
summary(md2.comm <- glm(SZ ~ pct.comm+F_MISS+Fh+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+dataset,data=MAIN2,family='binomial'))

#Do short or long ROHs drive this effect? In the raw data, long ROHs appear to drive the effect
summary(md2.short <- glm(SZ ~ I(proh.short.c2.1.out + proh.short.r2.1.out)+F_MISS+Fh+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+dataset,data=MAIN2,family='binomial'))
summary(md2.long <- glm(SZ ~ I(proh.big.c2.1.out + proh.big.r2.1.out)+F_MISS+Fh+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+dataset,data=MAIN2,family='binomial'))
###############################################









################################################
#13) RESULTS SEPARATED BY DATASET - all data points, IMPUTED data
#Look at results across datasets, controlling for F & F_MISS
MAIN2$pct.res[! is.na(MAIN2$PC1)] <- lm(tot2.1.imp ~ F_MISS+Fh.imp+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20,data=MAIN2)$residuals

#create matrix to put results in
GLM.RESULTS <- as.data.frame(matrix(NA,nrow=length(DATS),ncol=5,dimnames=list(DATS,c('Est','StErr','tval','df','platform'))))
GLM.RESULTS[,5] <-  c("A5","I550","A6","A500","A500","I650","A6","A6","A6","I317","A5","A5","A6","A6","I550","A5","A500")
kk <- 0
for (DD in DATS){
  kk <- kk+1
x <- summary(glm(SZ~pct.res,data=MAIN2[MAIN2$dataset==DD,],family='binomial'))
GLM.RESULTS[kk,1:3] <- x$coefficients[2,1:3]
GLM.RESULTS[kk,4] <- x$df[2]
}

#put overall result in 
mod.reg <- summary(glmer(SZ ~ tot2.1.imp+F_MISS+Fh.imp+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+(1|dataset),data=MAIN2,family='binomial'))
GLM.RESULTS2 <- rbind.data.frame(GLM.RESULTS,c(mod.reg@coefs[2,1:3],21856,NA))
GLM.RESULTS2[nrow(GLM.RESULTS2),ncol(GLM.RESULTS2)] <- 'xxx'
rownames(GLM.RESULTS2)[c(4,nrow(GLM.RESULTS2))] <- c("wtc","Total")

t.crits <- qt(.975,GLM.RESULTS2$df)
GLM.RESULTS2$low.bound <- GLM.RESULTS2$Est - t.crits*GLM.RESULTS2$StErr
GLM.RESULTS2$up.bound <- GLM.RESULTS2$Est + t.crits*GLM.RESULTS2$StErr

#include correct SE from bootstrap
#GLM.RESULTS2["Total",c(6,7)] <- c(.047697,.203393)

#PLOT results
#png("Froh.Results1.png",height=800,width=1000)
op <- par(mfrow=c(1,1),xpd=TRUE)

#plot the betas & SE(betas)
mycols <- c('red','blue','orange','green','purple','yellow','black')[as.numeric(as.factor(GLM.RESULTS2$platform))]
plot(GLM.RESULTS2$Est,ylim=c(-120,200),axes=FALSE,xlab='PGC Dataset',ylab='Slope of %ROH predicting SZ',col=mycols,main='Slopes and 95% CIs of SZ~%ROH across 17 PGC datasets, imputed data',pch=20)
axis(2)
axis(1,at=1:nrow(GLM.RESULTS2),labels=rownames(GLM.RESULTS2))

for (jj in 1:nrow(GLM.RESULTS2)){lines(x=c(jj,jj),y=c(GLM.RESULTS2$low.bound[jj],GLM.RESULTS2$up.bound[jj]),col=mycols[jj],lwd=2)}
legend(1,200,c('Affy5','Affy500','Affy6','Illum317','Illum550','Illum650'),col= c('red','blue','orange','green','purple','yellow','black'),lty=1,pch=20,title='SNP Platform')
par(xpd=FALSE)
abline(h=0,lty=1,lwd=.5)
abline(h=GLM.RESULTS2["Total","Est"],lwd=.5,lty=2)
text(1:nrow(GLM.RESULTS2),y=-125,GLM.RESULTS2$df+2)
text(0.5,y=-125,'n:')
#dev.off()

################################################








################################################
#14) RESULTS SEPARATED BY DATASET - no cous-cous outliers, IMPUTED data
#Look at results across datasets, controlling for F & F_MISS,
GLM.RESULTS.F <- as.data.frame(matrix(NA,nrow=length(DATS),ncol=5,dimnames=list(DATS,c('Est','StErr','tval','df','platform'))))
GLM.RESULTS.F[,5] <-  c("A5","I550","A6","A500","A500","I650","A6","A6","A6","I317","A5","A5","A6","A6","I550","A5","A500")
kk <- 0
for (DD in DATS){
  kk <- kk+1
x <- summary(glm(SZ~pct.res,data=MAIN2[MAIN2$dataset==DD & MAIN2$tot2.1.imp < .0625+mean(MAIN2$tot2.1),],family='binomial'))
GLM.RESULTS.F[kk,1:3] <- x$coefficients[2,1:3]
GLM.RESULTS.F[kk,4] <- x$df[2]
}

#put overall result in 
mod.reg2 <- summary(glmer(SZ ~ tot2.1.imp+F_MISS+Fh.imp+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+(1|dataset),data=MAIN2[MAIN2$tot2.1.imp < .0625+mean(MAIN2$tot2.1),],family='binomial'))
GLM.RESULTS.F2 <- rbind.data.frame(GLM.RESULTS.F,c(mod.reg2@coefs[2,1:3],21854,NA))
GLM.RESULTS.F2[nrow(GLM.RESULTS.F2),ncol(GLM.RESULTS.F2)] <- 'xxx'
rownames(GLM.RESULTS.F2)[c(4,nrow(GLM.RESULTS.F2))] <- c("wtc","Total")

t.crits <- qt(.975,GLM.RESULTS.F2$df)
GLM.RESULTS.F2$low.bound <- GLM.RESULTS.F2$Est - t.crits*GLM.RESULTS.F2$StErr
GLM.RESULTS.F2$up.bound <- GLM.RESULTS.F2$Est + t.crits*GLM.RESULTS.F2$StErr

#include correct SE from bootstrap
#GLM.RESULTS.F2["Total",c(6,7)] <- c(.047697,.203393)

#PLOT results
#png("Froh.Results1.png",height=800,width=1000)
op <- par(mfrow=c(1,1),xpd=TRUE)

#plot the betas & SE(betas)
mycols <- c('red','blue','orange','green','purple','yellow','black')[as.numeric(as.factor(GLM.RESULTS.F2$platform))]
plot(GLM.RESULTS.F2$Est,ylim=c(-120,200),axes=FALSE,xlab='PGC Dataset',ylab='Slope of %ROH predicting SZ',col=mycols,main='Slopes and 95% CIs of SZ~%ROH across 17 PGC datasets, imputed data',pch=20)
axis(2)
axis(1,at=1:nrow(GLM.RESULTS.F2),labels=rownames(GLM.RESULTS.F2))

for (jj in 1:nrow(GLM.RESULTS.F2)){lines(x=c(jj,jj),y=c(GLM.RESULTS.F2$low.bound[jj],GLM.RESULTS.F2$up.bound[jj]),col=mycols[jj],lwd=2)}
legend(1,200,c('Affy5','Affy500','Affy6','Illum317','Illum550','Illum650'),col= c('red','blue','orange','green','purple','yellow','black'),lty=1,pch=20,title='SNP Platform')
par(xpd=FALSE)
abline(h=0,lty=1,lwd=.5)
abline(h=GLM.RESULTS.F2["Total","Est"],lwd=.5,lty=2)
text(1:nrow(GLM.RESULTS.F2),y=-125,GLM.RESULTS.F2$df+2)
text(0.5,y=-125,'n:')
#dev.off()

################################################








################################################
#15) RESULTS SEPARATED BY DATASET - all data points, RAW data
#Look at results across datasets, controlling for F & F_MISS
MAIN2$pct.res[! is.na(MAIN2$PC1)] <- lm(tot2.1 ~ F_MISS+Fh+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20,data=MAIN2)$residuals

#create matrix to put results in
GLM.RESULTS.R <- as.data.frame(matrix(NA,nrow=length(DATS),ncol=5,dimnames=list(DATS,c('Est','StErr','tval','df','platform'))))
GLM.RESULTS.R[,5] <-  c("A5","I550","A6","A500","A500","I650","A6","A6","A6","I317","A5","A5","A6","A6","I550","A5","A500")
kk <- 0
for (DD in DATS){
  kk <- kk+1
x <- summary(glm(SZ~pct.res,data=MAIN2[MAIN2$dataset==DD,],family='binomial'))
GLM.RESULTS.R[kk,1:3] <- x$coefficients[2,1:3]
GLM.RESULTS.R[kk,4] <- x$df[2]
}

#put overall result in 
mod.reg <- summary(glmer(SZ ~ tot2.1+F_MISS+Fh+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+(1|dataset),data=MAIN2,family='binomial'))
GLM.RESULTS.R2 <- rbind.data.frame(GLM.RESULTS.R,c(mod.reg@coefs[2,1:3],21854,NA))
GLM.RESULTS.R2[nrow(GLM.RESULTS.R2),ncol(GLM.RESULTS.R2)] <- 'xxx'
rownames(GLM.RESULTS.R2)[c(4,nrow(GLM.RESULTS.R2))] <- c("wtc","Total")

t.crits <- qt(.975,GLM.RESULTS.R2$df)
GLM.RESULTS.R2$low.bound <- GLM.RESULTS.R2$Est - t.crits*GLM.RESULTS.R2$StErr
GLM.RESULTS.R2$up.bound <- GLM.RESULTS.R2$Est + t.crits*GLM.RESULTS.R2$StErr

#include correct SE from bootstrap
#GLM.RESULTS.R2["Total",c(6,7)] <- c(.047697,.203393)


#PLOT results
#png("Froh.Results1.png",height=800,width=1000)
op <- par(mfrow=c(1,1),xpd=TRUE)

#plot the betas & SE(betas)
mycols <- c('red','blue','orange','green','purple','yellow','black')[as.numeric(as.factor(GLM.RESULTS.R2$platform))]
plot(GLM.RESULTS.R2$Est,ylim=c(-120,200),axes=FALSE,xlab='PGC Dataset',ylab='Slope of %ROH predicting SZ',col=mycols,main='Slopes and 95% CIs of SZ~%ROH across 17 PGC datasets, imputed data',pch=20)
axis(2)
axis(1,at=1:nrow(GLM.RESULTS.R2),labels=rownames(GLM.RESULTS.R2))

for (jj in 1:nrow(GLM.RESULTS.R2)){lines(x=c(jj,jj),y=c(GLM.RESULTS.R2$low.bound[jj],GLM.RESULTS.R2$up.bound[jj]),col=mycols[jj],lwd=2)}
legend(1,200,c('Affy5','Affy500','Affy6','Illum317','Illum550','Illum650'),col= c('red','blue','orange','green','purple','yellow','black'),lty=1,pch=20,title='SNP Platform')
par(xpd=FALSE)
abline(h=0,lty=1,lwd=.5)
abline(h=GLM.RESULTS.R2["Total","Est"],lwd=.5,lty=2)
text(1:nrow(GLM.RESULTS.R2),y=-125,GLM.RESULTS.R2$df+2)
text(0.5,y=-125,'n:')
#dev.off()

################################################








################################################
#16) EFFECT SIZE
#get Nagelkerke's R-squared
require("fmsb")
summary(reg.mod.full <-glm(SZ ~ tot2.1.imp+F_MISS+Fh.imp+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+dataset,data=MAIN2,family='binomial'))
summary(reg.mod.red <- glm(SZ ~ F_MISS+Fh.imp+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+dataset,data=MAIN2,family='binomial'))
sqrt(NagelkerkeR2(reg.mod.full)$R2 - NagelkerkeR2(reg.mod.red)$R2)

#predict cousin-cousin effect
summary(mod.mn <- glmer(SZ ~ tot2.1.imp+F_MISS+Fh.imp+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+(1|dataset),data=MAIN2,family='binomial'))
odds.null <- exp(fixef(mod.mn) %*% c(1,mean(MAIN2$tot2.1.imp,na.rm=T),mean(MAIN2$F_MISS,na.rm=T),mean(MAIN2$Fh.imp,na.rm=T),rep(0,20)))  #odds of having SZ if mean %ROH & mean on all covariates = .9073
odds.cous <- exp(fixef(mod.mn) %*% c(1,mean(MAIN2$tot2.1.imp,na.rm=T)+.0625,mean(MAIN2$F_MISS,na.rm=T),mean(MAIN2$Fh.imp,na.rm=T),rep(0,20)))  #odds of having SZ if mean %ROH & mean on all covariates = 2.403
odds.cous/odds.null  #2.65 increase in odds
################################################









################################################
#17) Look at just the extremely long ROHs - are these over-represented among cases?
#Answer: Yes
(case.percent.obs <- mean(MAIN2$SZ))


#For 10Mb ROHs
LENGTH <- 10000
(pct.cases <- sum(HOM_I[HOM_I$KB>LENGTH,'case']=='case')/sum(HOM_I$KB>LENGTH))
(obs.cases <- sum(HOM_I[HOM_I$KB>LENGTH,'case']=='case'))
(exp.cases <- case.percent.obs*sum(HOM_I$KB>LENGTH))
pbinom(q=obs.cases,size=sum(HOM_I$KB>LENGTH),prob=case.percent.obs,lower.tail=FALSE,log.p=FALSE)
x <- HOM_I[HOM_I$KB>LENGTH,c(2,4,7:10,14,16,17)]
#x[order(x$KB),]


#For 16.66Mb ROHs (expected size of 2nd cous inbreeding)
LENGTH <- 16667
(pct.cases <- sum(HOM_I[HOM_I$KB>LENGTH,'case']=='case')/sum(HOM_I$KB>LENGTH))
(obs.cases <- sum(HOM_I[HOM_I$KB>LENGTH,'case']=='case'))
(exp.cases <- case.percent.obs*sum(HOM_I$KB>LENGTH))
(tot.num <- sum(HOM_I$KB>LENGTH) )
pbinom(q=obs.cases,size=tot.num,prob=case.percent.obs,lower.tail=FALSE,log.p=FALSE)
x <- HOM_I[HOM_I$KB>LENGTH,c(1,2,4,7:10,14,16,17)]
#x[order(x$KB),]

#correct above p-value by just looking at unique individuals rather than by ROH
(pct.cases <- sum(z[z$KB>LENGTH,'case']=='case')/sum(z$KB>LENGTH))
(obs.cases <- sum(z[z$KB>LENGTH,'case']=='case'))
(exp.cases <- case.percent.obs*sum(z$KB>LENGTH))
(tot.num <- sum(z$KB>LENGTH) )
pbinom(q=obs.cases,size=tot.num,prob=case.percent.obs,lower.tail=FALSE,log.p=FALSE)
################################################







################################################
#18) Look at those with %ROH > 5%

#Look at cases/controls who have %ROH > 5%
inb <- MAIN2[MAIN2$tot2.1.imp>.05,]
(pct.cases <- sum(inb$SZ==1)/nrow(inb));case.percent.obs
(obs.cases <- sum(inb$SZ==1))
(exp.cases <- case.percent.obs*nrow(inb))
(tot.num <- nrow(inb))
pbinom(q=obs.cases,size=tot.num,prob=case.percent.obs,lower.tail=FALSE,log.p=FALSE)
barplot(c(case.percent.obs,pct.cases),col=c('blue','red'),names.arg=c('% cases overall','% cases with Froh > 5%'),ylim=c(0,1))





#Look at distribution of ROHs among the most inbred individuals
HOM_I$UID <- paste(HOM_I$FID,HOM_I$IID,sep='')
INBRED <- HOM_I[HOM_I$UID %in% MAIN2$UID[MAIN2$tot2.1.imp>.05],]
INBRED <- INBRED[order(INBRED$UID),]

#make the proper lengths of ROHs:
INBRED$POS1x <- as.numeric(paste(INBRED$CHR*1000,INBRED$POS1,sep=''))
INBRED$POS2x <- as.numeric(paste(INBRED$CHR*1000,INBRED$POS2,sep=''))
pos1 <- INBRED$POS1x[2:nrow(INBRED)]
pos2 <- INBRED$POS2x[1:(nrow(INBRED)-1)]
uid1 <- INBRED$UID[2:nrow(INBRED)]
uid2 <- INBRED$UID[1:(nrow(INBRED)-1)]
cont <- (pos1-pos2) < 1e6  & (pos1-pos2)>0 & uid1==uid2    #if two 'separate' ROHs within a person are within 1 Mb of each other, they are called 'the same'
INBRED$cont1 <- c(0,cont)*1
INBRED$cont2 <- c(0,0,cont[-length(cont)])*1
INBRED$cont3 <- (INBRED$cont1+INBRED$cont2==2)*1
INBRED$POS1n <- round(INBRED$POS1/1e6,2)
INBRED$POS2n <- round(INBRED$POS2/1e6,2)
INBRED[,c(1,4,26,27,23,24,25)]   #look



write.table(INBRED,file='inbred',quote=FALSE,row.names=FALSE,col.names=TRUE,sep=',') #write this out and do it by hand in excel (can't think how to program it)

#grab the modified file that has the ROH lenghts of all people with Froh > 5
x <- read.table("~/Documents/Academics/Presentations/PGC_SZ/inbred2.txt",header=TRUE)
x$ROH.length <- (x$POS2 - x$POS1)/1e6


hist(x$ROH.length[x$ROH.length<50],breaks=20,probability=TRUE,main='',xlab='',col='lightblue1')
op <- par(new=TRUE)
plot(seq(0,50,length.out=100),fun(seq(0,50,length.out=100),lmbd=1/25),type='l',xlab='',ylab='',axes=FALSE,col='red',lwd=3,ylim=c(0,.14))
par(op)



################################################










































################################################
#11) Is Fh a measure of stratification & technical artifacts not picked up by C1:C20? Doesn't appear to be
F.mod1 <- lm(Fh~dataset,data=TOT.med3)
F.mod2 <- lm(Fh~C1+C2+C3+C4+C5+C6+C8+C9+C10+C11+C12+C13+C14+C15+C16+C18+C19+C20,data=TOT.med3)
F.mod3 <- lm(Fh~dataset+C1+C2+C3+C4+C5+C6+C8+C9+C10+C11+C12+C13+C14+C15+C16+C18+C19+C20,data=TOT.med3)
anova(F.mod2,F.mod3)
anova(F.mod1)
#The effect of dataset on Fh is stronger after controlling for C1:C20

summary(F.mod1 <- lm(Fh~platform,data=TOT.med3))
F.mod2 <- lm(Fh~C1+C2+C3+C4+C5+C6+C8+C9+C10+C11+C12+C13+C14+C15+C16+C18+C19+C20,data=TOT.med3)
summary(F.mod3 <- lm(Fh~platform+C1+C2+C3+C4+C5+C6+C8+C9+C10+C11+C12+C13+C14+C15+C16+C18+C19+C20,data=TOT.med3))
anova(F.mod2,F.mod3)
anova(F.mod1)

#What about %ROH?
summary(Pct.mod1 <- lm(pct~dataset,data=TOT.med3))
Pct.mod2 <- lm(pct~C1+C2+C3+C4+C5+C6+C8+C9+C10+C11+C12+C13+C14+C15+C16+C18+C19+C20,data=TOT.med3)
summary(Pct.mod3 <- lm(pct~dataset+C1+C2+C3+C4+C5+C6+C8+C9+C10+C11+C12+C13+C14+C15+C16+C18+C19+C20,data=TOT.med3))
anova(Pct.mod2,Pct.mod3)
anova(Pct.mod1)
#The effect of dataset on pct is weaker after controlling for C1:C20

summary(Pct.mod1 <- lm(pct~platform,data=TOT.med3))
Pct.mod2 <- lm(pct~C1+C2+C3+C4+C5+C6+C8+C9+C10+C11+C12+C13+C14+C15+C16+C18+C19+C20,data=TOT.med3)
summary(Pct.mod3 <- lm(pct~platform+C1+C2+C3+C4+C5+C6+C8+C9+C10+C11+C12+C13+C14+C15+C16+C18+C19+C20,data=TOT.med3))
anova(Pct.mod2,Pct.mod3)
anova(Pct.mod1)


id <- TOT.med3$IID.x[seq(1,nrow(TOT.med3),length.out=400)]
fid <- TOT.med3$FID.x[seq(1,nrow(TOT.med3),length.out=400)]
dats <- TOT.med3$dataset[seq(1,nrow(TOT.med3),length.out=400)]
cs <- TOT.med3$case[seq(1,nrow(TOT.med3),length.out=400)]
cbind.data.frame(id,fid,dats,cs)
################################################








################################################
#12) PLATE EFFECTS
#Look at plate effects for those 7000 subjects that we have plate information on
TOT.med3.plates <- TOT.med3[TOT.med3$dataset %in% c('ab','bulg','dub','edi','port','sw1','sw2','ucl','mgs2'),]
TOT.med3.plates <- TOT.med3.plates[! (TOT.med3.plates$dataset=='ucl' & TOT.med3.plates$case=='control'),]
fff <- strsplit(as.character(TOT.med3.plates$IID.x),"_",fixed=TRUE)
TOT.med3.plates$plate <- as.factor(unlist(lapply(fff,function(x) x[1])))

summary(as.factor(TOT.med3.plates$plate),maxsum=10000)
length(summary(as.factor(TOT.med3.plates$plate),maxsum=10000)) #87 plates


#Pull in plate info on MGS:
mgs.plates <- read.table("~/HD3/SNP.Homozygosity/PGC.SZ/MGS.plates/MGS.plates2",header=TRUE)
gain.plates <- read.table("~/HD3/SNP.Homozygosity/PGC.SZ/MGS.plates/GAIN.plates2",header=TRUE)
colnames(gain.plates) <- colnames(mgs.plates)
mgs.plt <- rbind(mgs.plates[,c(1,3)],gain.plates[,c(1,3)])
TOT.med3.plates <- merge(TOT.med3.plates,mgs.plt,by.x='IID.x',by.y=1,all.x=TRUE,all.y=FALSE)
TOT.med3.plates$plate <- as.character(TOT.med3.plates$plate)
TOT.med3.plates$Plate <- as.character(TOT.med3.plates$Plate)
TOT.med3.plates[TOT.med3.plates$dataset=='mgs2','plate'] <- TOT.med3.plates[TOT.med3.plates$dataset=='mgs2','Plate']
TOT.med3.plates$plate <- as.factor(TOT.med3.plates$plate)

TOT.med3[TOT.med3$dataset=='mgs2','IID.x'][seq(1,5000,length.out=200)]
mgs.plt[,1][seq(1,5000,length.out=200)]


#look at case/control mix on each plate; pretty good mix on each plate
TOT.med3.plates$plate.case <- paste(TOT.med3.plates$dataset,TOT.med3.plates$plate,TOT.med3.plates$case,sep=".")
summary(as.factor(TOT.med3.plates$plate.case),maxsum=10000)

#pct is more affected by plate than Fh. This is probably because type 2 errors are highly sensitive to quality
summary(lm(Fh~plate,data=TOT.med3.plates))
summary(lm(pct~plate,data=TOT.med3.plates))
summary(lm(pct~plate+Fh,data=TOT.med3.plates))

summary(glm(SZ~pct+C1+C2+C3+C4+C5+C6+C8+C9+C14+C15+C20,data=TOT.med3.plates))
summary(glm(SZ~pct+C1+C2+C3+C4+C5+C6+C8+C9+C14+C15+C20+plate,data=TOT.med3.plates))
summary(glm(SZ~pct+Fh+C1+C2+C3+C4+C5+C6+C8+C9+C14+C15+C20+plate,data=TOT.med3.plates))

summary(mmm0 <- lm(pct~1,data=TOT.med3.plates))
summary(mmm1 <- lm(pct~Fh,data=TOT.med3.plates))
summary(mmm2 <- lm(pct~C1+C2+C3+C4+C5+C6+C8+C9+C14+C15+C20 + dataset + plate,data=TOT.med3.plates))
summary(mmm3 <- lm(pct~C1+C2+C3+C4+C5+C6+C8+C9+C14+C15+C20 + dataset + plate + Fh,data=TOT.med3.plates))

.3982^2-.3573^2 ; .4183^2 - .3746^2 #Fh explains abt as much var in %ROH whether control for covariates or not

#add plates to TOT.med3
TOT.med3.plates2 <- TOT.med3.plates[,c('UID','plate')]
TOT.med4 <- merge(TOT.med3,TOT.med3.plates2,by='UID',all.x=TRUE,all.y=FALSE)
TOT.med4$plate <- as.character(TOT.med4$plate)
TOT.med4$plate[is.na(TOT.med4$plate)] <- 'unknown'
TOT.med4$plate <- as.factor(TOT.med4$plate)


#Look at models with and without outliers, including plates
summary(glm(SZ ~ pct + C1+C2+C3+C4+C5+C6+C8+C9+C14+C15+C20+ dataset,data=TOT.med4,family='binomial'))
summary(glm(SZ ~ pct + C1+C2+C3+C4+C5+C6+C8+C9+C14+C15+C20+ dataset+plate,data=TOT.med4,family='binomial'))

#Fh
summary(reg.mod.full.F <- glm(SZ ~ pct +Fh+ C1+C2+C3+C4+C5+C6+C8+C9+C14+C15+C20+ dataset+plate,data=TOT.med4,family='binomial'))

#Separated by rare vs. common and short vs. long
summary(reg.mod.full.F.rare <- glm(SZ ~ pct.long.rare +Fh+ C1+C2+C3+C4+C5+C6+C8+C9+C14+C15+C20+ dataset+plate,data=TOT.med4,family='binomial'))
summary(reg.mod.full.F.rare <- glm(SZ ~ pct.short.rare +Fh+ C1+C2+C3+C4+C5+C6+C8+C9+C14+C15+C20+ dataset+plate,data=TOT.med4,family='binomial'))
summary(reg.mod.full.F.comm <- glm(SZ ~ pct.long.comm + Fh+C1+C2+C3+C4+C5+C6+C8+C9+C14+C15+C20+ dataset+plate,data=TOT.med4,family='binomial'))
summary(reg.mod.full.F.comm <- glm(SZ ~ pct.short.comm + Fh+C1+C2+C3+C4+C5+C6+C8+C9+C14+C15+C20+ dataset+plate,data=TOT.med4,family='binomial'))



#LEARNED:
#1) So Fh is capturing something more than plates, platforms, and dataset differences; its mediating effect on %ROH is above and beyond these three potentially technical artifacts
#2) Plates have a big effect on %ROH, and controlling for them removes a lot of noise, making the effect
# significant in even these datasets, where the effect is quite small
################################################






################################################
#13) Look at just the extremely long ROHs - are these over-represented among cases?
#Answer: Yes, they appear to be, but the p-vals here can't be trusted bc non-indep & not controlling for dataset
(case.percent <- mean(TOT.med3$SZ))

#For 1Mb ROHs
LENGTH <- 1000
(pct.cases <- sum(HOM[HOM$KB>LENGTH,'case']=='case')/sum(HOM$KB>LENGTH))
(obs.cases <- sum(HOM[HOM$KB>LENGTH,'case']=='case'))
exp.cases <- case.percent*sum(HOM$KB>LENGTH)
pbinom(q=obs.cases,size=sum(HOM$KB>LENGTH),prob=case.percent,lower.tail=FALSE,log.p=FALSE)


#For 2.5Mb ROHs
LENGTH <- 2500
(pct.cases <- sum(HOM[HOM$KB>LENGTH,'case']=='case')/sum(HOM$KB>LENGTH))
(obs.cases <- sum(HOM[HOM$KB>LENGTH,'case']=='case'))
exp.cases <- case.percent*sum(HOM$KB>LENGTH)
pbinom(q=obs.cases,size=sum(HOM$KB>LENGTH),prob=case.percent,lower.tail=FALSE,log.p=FALSE)


#For 5Mb ROHs
LENGTH <- 5000
(pct.cases <- sum(HOM[HOM$KB>LENGTH,'case']=='case')/sum(HOM$KB>LENGTH))
(obs.cases <- sum(HOM[HOM$KB>LENGTH,'case']=='case'))
exp.cases <- case.percent*sum(HOM$KB>LENGTH)
pbinom(q=obs.cases,size=sum(HOM$KB>LENGTH),prob=case.percent,lower.tail=FALSE,log.p=FALSE)


#For 10Mb ROHs
LENGTH <- 10000
(pct.cases <- sum(HOM[HOM$KB>LENGTH,'case']=='case')/sum(HOM$KB>LENGTH))
(obs.cases <- sum(HOM[HOM$KB>LENGTH,'case']=='case'))
exp.cases <- case.percent*sum(HOM$KB>LENGTH)
pbinom(q=obs.cases,size=sum(HOM$KB>LENGTH),prob=case.percent,lower.tail=FALSE,log.p=FALSE)
x <- HOM[HOM$KB>LENGTH,c(2,4,7:10,14,16,17)]
x[order(x$KB),]


#For 16.66Mb ROHs (expected size of 2nd cous inbreeding)
LENGTH <- 16667
(pct.cases <- sum(HOM[HOM$KB>LENGTH,'case']=='case')/sum(HOM$KB>LENGTH))
(obs.cases <- sum(HOM[HOM$KB>LENGTH,'case']=='case'))
exp.cases <- case.percent*sum(HOM$KB>LENGTH)
pbinom(q=obs.cases,size=sum(HOM$KB>LENGTH),prob=case.percent,lower.tail=FALSE,log.p=FALSE)
x <- HOM[HOM$KB>LENGTH,c(2,4,7:10,14,16,17)]
x[order(x$KB),]
z <- x[! duplicated(x$UID),]
pbinom(q=115,size=nrow(z),prob=case.percent,lower.tail=FALSE,log.p=FALSE)

barplot(c(.43,.667),col=c('blue','red'))





#For 20Mb ROHs
LENGTH <- 20000
(pct.cases <- sum(HOM[HOM$KB>LENGTH,'case']=='case')/sum(HOM$KB>LENGTH))
(obs.cases <- sum(HOM[HOM$KB>LENGTH,'case']=='case'))
exp.cases <- case.percent*sum(HOM$KB>LENGTH)
pbinom(q=obs.cases,size=sum(HOM$KB>LENGTH),prob=case.percent,lower.tail=FALSE,log.p=FALSE)
x <- HOM[HOM$KB>LENGTH,c(1,2,4,7:10,14,16,17)]
x$UID <- paste(x$FID,x$IID,sep='.')
length(unique(x$UID))

x <- HOM[HOM$KB>LENGTH,c(1,2,4,7:10,14,16,17)]
x$UID <- paste(x[x$KB>LENGTH,'FID'],x[x$KB>LENGTH,'IID'],sep='')
z <- x[! duplicated(x$UID),]

pbinom(q=sum(z$case=='case'),size=nrow(z),prob=case.percent,lower.tail=FALSE,log.p=FALSE)



#grab the modified file that has the ROH lenghts of all people with Froh > 5
x <- read.table("~/Documents/Academics/Presentations/PGC_SZ/inbred2.txt",header=TRUE)
x$ROH.length <- (x$POS2 - x$POS1)/1e6


hist(x$ROH.length[x$ROH.length<50],breaks=20,probability=TRUE,main='',xlab='',col='lightblue1')
op <- par(new=TRUE)
plot(seq(0,50,length.out=100),fun(seq(0,50,length.out=100),lmbd=1/25),type='l',xlab='',ylab='',axes=FALSE,col='red',lwd=3,ylim=c(0,.14))
par(op)



x <- fun(seq(0,50,length.out=1000),lmbd=1/4)

mnx <- sum(x*seq(0,50,length.out=1000))



plot(seq(0,5,length.out=1000),fun(seq(0,5,length.out=1000),lmbd=.5),type='l',xlab='',ylab='',axes=T,col='orange',lwd=3)
par(new=TRUE)
plot(seq(0,5,length.out=1000),fun(seq(0,5,length.out=1000),lmbd=1),type='l',xlab='',ylab='',axes=T,col='purple',lwd=3)
par(new=TRUE)
plot(seq(0,5,length.out=1000),fun(seq(0,5,length.out=1000),lmbd=1.5),type='l',xlab='',ylab='',axes=T,col='blue',lwd=3)



fun(seq(0,5,length.out=10),lmbd=.25)
dexp(seq(0,5,length.out=10),rate=.25)

x <- rexp(1000,rate=.25)


#For 40Mb ROHs
LENGTH <- 40000
(pct.cases <- sum(HOM[HOM$KB>LENGTH,'case']=='case')/sum(HOM$KB>LENGTH))
(obs.cases <- sum(HOM[HOM$KB>LENGTH,'case']=='case'))
exp.cases <- case.percent*sum(HOM$KB>LENGTH)
pbinom(q=obs.cases,size=sum(HOM$KB>LENGTH),prob=case.percent,lower.tail=FALSE,log.p=FALSE)
x <- HOM2[HOM2$KB>LENGTH,c(2,4,7:10,14,16,17)]


summary(glm(SZ~NSEG.1+dataset+C1+C2+C3+C4+C5+C6+C8+C9+C14+C15+C20,family='binomial',data=TOT.med3))
summary(glm(SZ~NSEG.1+C1+C2+C3+C4+C5+C6+C8+C9+C14+C15+C20,family='binomial',data=TOT.med3))



sum(TOT.med3$pct > 6.)


#Look at distribution of ROHs among the most inbred individuals
HOM$UID <- paste(HOM$FID,HOM$IID,sep='')
INBRED <- HOM[HOM$UID %in% TOT.med3$UID[TOT.med3$pct > 5],]

#make the proper lengths of ROHs:
INBRED$POS1x <- as.numeric(paste(INBRED$CHR,INBRED$POS1,sep=''))
INBRED$POS2x <- as.numeric(paste(INBRED$CHR,INBRED$POS2,sep=''))
pos1 <- INBRED$POS1x[2:nrow(INBRED)]
pos2 <- INBRED$POS2x[1:(nrow(INBRED)-1)]
cont <- abs(pos1-pos2) < 500000
INBRED$cont <- c(cont,0)*1
INBRED$cont2 <- c(0,cont)*1
INBRED$cont3 <- (INBRED$cont+INBRED$cont2>0)*1

write.table(INBRED,file='inbred',quote=FALSE,row.names=FALSE,col.names=TRUE,sep=',')

#after modifying the file by hand
inb2 <- read.table('inbred.txt',


#case-control p-value for the most inbred individuals
pbinom(q=14,size=21,prob=case.percent,lower.tail=FALSE,log.p=FALSE)
################################################









################################################
#14) Get bootstrap CI
require(foreach)
require(doMC)
registerDoMC(cores=10)


num.iter <- 1000
PERM.REG <- matrix(NA,nrow=num.iter,ncol=4,dimnames=list(NULL,c('Est','StErr','z.val','p.val')))
PERM.REG.F <- matrix(NA,nrow=num.iter,ncol=4,dimnames=list(NULL,c('Est','StErr','z.val','p.val')))
PERM.NO <- matrix(NA,nrow=num.iter,ncol=4,dimnames=list(NULL,c('Est','StErr','z.val','p.val')))
PERM.NO.F <- matrix(NA,nrow=num.iter,ncol=4,dimnames=list(NULL,c('Est','StErr','z.val','p.val')))



TOT.boot <- TOT.med3[order(TOT.med3$dataset),c('SZ','pct','Fh','C1','C2','C3','C4','C5','C6','C8','C9','C14','C15','C20','dataset','dataset.case')]
st <- which(! duplicated(TOT.boot$dataset.case))
en <- which(! duplicated(TOT.boot$dataset.case,fromLast=TRUE))

PERM.REG <- foreach (T = 1:num.iter,.combine='rbind') %dopar% {
ind <- vector()
for (i in 1: length(en)) {ind <- c(ind,sample(st[i]:en[i],size=en[i]-st[i]+1,replace=TRUE))}
boot <- TOT.boot[ind,]
summary(reg.mod.full <- glmer(SZ ~ pct + C1+C2+C3+C4+C5+C6+C8+C9+C14+C15+C20+ (1|dataset),data=boot,family='binomial'))@coefs[2,]
}


PERM.REG.F <- foreach (T = 1:num.iter,.combine='rbind') %dopar% {
ind <- vector()
for (i in 1: length(en)) {ind <- c(ind,sample(st[i]:en[i],size=en[i]-st[i]+1,replace=TRUE))}
boot <- TOT.boot[ind,]
 summary(reg.mod.fh <- glmer(SZ ~ pct +Fh+ C1+C2+C3+C4+C5+C6+C8+C9+C14+C15+C20+ (1|dataset),data=boot,family='binomial'))@coefs[2,]
}



#Bootstrap for no outlier data
TOT.boot.no <- TOT.med3.no[order(TOT.med3.no$dataset),c('SZ','pct','Fh','C1','C2','C3','C4','C5','C6','C8','C9','C14','C15','C20','dataset','dataset.case')]
st <- which(! duplicated(TOT.boot.no$dataset.case))
en <- which(! duplicated(TOT.boot.no$dataset.case,fromLast=TRUE))

PERM.NO <- foreach (T = 1:num.iter,.combine='rbind') %dopar% {
ind <- vector()
for (i in 1: length(en)) {ind <- c(ind,sample(st[i]:en[i],size=en[i]-st[i]+1,replace=TRUE))}
boot.no <- TOT.boot.no[ind,]
summary(reg.mod.full.no <- glmer(SZ ~ pct + C1+C2+C3+C4+C5+C6+C8+C9+C14+C15+C20+ (1|dataset),data=boot.no,family='binomial'))@coefs[2,]
}

PERM.NO.F <- foreach (T = 1:num.iter,.combine='rbind') %dopar% {
ind <- vector()
for (i in 1: length(en)) {ind <- c(ind,sample(st[i]:en[i],size=en[i]-st[i]+1,replace=TRUE))}
boot.no <- TOT.boot.no[ind,]
summary(reg.mod.fh.no <- glmer(SZ ~ pct +Fh+ C1+C2+C3+C4+C5+C6+C8+C9+C14+C15+C20+ (1|dataset),data=boot.no,family='binomial'))@coefs[2,]
}



summary(glmer(SZ ~ pct + C1+C2+C3+C4+C5+C6+C8+C9+C14+C15+C20+ (1|dataset),data=TOT.med3,family='binomial'))
mean(PERM.REG[,1])
sd(PERM.REG[,1])
(sum(PERM.REG[,1]<0)/1000)*2  #bootstrap p-value = .002
hist(PERM.REG[,1],100)
preg <- sort(PERM.REG[,1])
c(preg[25],preg[975])  #95% CI

mean(PERM.REG.F[,1])
sd(PERM.REG.F[,1])
(sum(PERM.REG.F[,1]<0)/1000)*2  #bootstrap p-value < .001
hist(PERM.REG.F[,1],100)
pregf <- sort(PERM.REG.F[,1])
c(pregf[25],pregf[975])  #95% CI

summary(glmer(SZ ~ pct + C1+C2+C3+C4+C5+C6+C8+C9+C14+C15+C20+ (1|dataset),data=TOT.med3.no,family='binomial'))
mean(PERM.NO[,1])
sd(PERM.NO[,1])
(sum(PERM.NO[,1]<0)/1000)*2  #bootstrap p-value = .022
hist(PERM.NO[,1],100)
pregn <- sort(PERM.REG[,1])
c(pregn[25],pregn[975])  #95% CI

mean(PERM.NO.F[,1])
sd(PERM.NO.F[,1])
(sum(PERM.NO.F[,1]<0)/1000)*2  #bootstrap p-value < .001
hist(PERM.NO.F[,1],100)
################################################










#look at duplicated individuals
#this was sent to Dan H for him to get error rates in his project, but might be useful for me later:

IND <- rbind(ab.ind,bon.ind,bulg.ind,carwtc.ind,cat2.ind,dk.ind,dub.ind,edi.ind,mgs2.ind,muc.ind,port.ind,sw1.ind,sw2.ind,top3.ind,ucla.ind,ucl.ind,zhh.ind)

IND <- IND[,1:3]
DDD <- DAT[,2:3]
DDD$analyzed.set <- 'yes'
IND2 <- merge(IND,DDD,by=c(1,2),all=TRUE)
fff <- strsplit(as.character(IND2$FID),"_",fixed=TRUE)
IND2$dataset <- as.factor(unlist(lapply(fff,function(x) x[3])))

IND3 <- IND2[IND2$dataset=='mgs2',]
save(IND3,file='mgs.dupd')

save(TOT.med3,file='~/HD3/SNP.Homozygosity/PGC.SZ/tot.RData')
