# MERGE phenotype for ROH analysis and residualize all covariates:
#
# model: HEIGHT ~ AGE, SEX, PCAS, BATCH


# Phenotype and covariate data for ARIC, MESA,and CARDIA:
setwd("/home/simonsom/ROH_pathway/TOTAL")
load("phe.model.Rdata")
fam.dat <- read.table("MERGE.clean.FINAL.fam",header=FALSE)
FID.IID <- fam.dat[,1:2]
names(FID.IID) <- c('FID','IID')
# remove NA's from total data
total <- total[(!is.na(total$IID) & !is.na(total$BMI) & !is.na(total$SEX) & !is.na(total$AGE) & !is.na(total$BATCH) & !is.na(total$C1)),]
# remove any duplicates:
dup.index <- duplicated(total$IID)
full.data <- total[which(dup.index==FALSE),]

# read in ARIC height:

phen.total <- read.csv("ARIC_CARE_pheno.GRU.csv",header=TRUE)
aric.full <- merge(full.data,phen.total,by="IID") # full ARIC phenotype data

aric.data <- cbind.data.frame(aric.full[,1:26],aric.full$anta01)
names(aric.data) <- c(names(aric.full)[1:26],'HEIGHT')

aric.ht.pheno <- 
# read in MESA height:

mesa.pheno <- read.table("MESA_pheno.txt",header=TRUE, fill=TRUE)
names(mesa.pheno) <- c('dbGAP','IID','HEIGHT')
mesa.full <- merge(full.data,mesa.pheno,by="IID") # full MESA phenotype data
mesa.data <- cbind.data.frame(mesa.full[,1:26],mesa.full$HEIGHT)
names(mesa.data) <- c(names(mesa.full)[1:26],'HEIGHT')

# read in CARDIA height:

cardia.pheno <- read.table("CARDIA_phe_ht.txt", header=TRUE,fill=TRUE)
names(cardia.pheno) <- c('dbGAP','IID','a','b','HEIGHT')
cardia.full <- merge(full.data,cardia.pheno,by="IID") # full MESA phenotype data
cardia.data <- cbind.data.frame(cardia.full[,1:26],cardia.full$HEIGHT)
names(cardia.data) <- c(names(cardia.full)[1:26],'HEIGHT')

# read in WHI height
setwd("/STATGEN/home/simonsom/ROH_pathway/TOTAL/WHI")
annot <- read.csv("Sample_annotation.csv",header=TRUE)
names(annot) <- c("SOURCE_SUBJID",names(annot)[2:length(names(annot))])
phe1 <- read.csv("WHI_AGE.csv",header=TRUE)
phe2 <- read.csv("WHI_height.csv",header=TRUE)
dup.index <- duplicated(phe2$dbGaP.SubjID)
phe2.u <- phe2[!dup.index,]

phe.tot <- merge(phe1,phe2,by="dbGaP.SubjID")
dup.tindex <- duplicated(phe.tot$dbGaP.SubjID)
phe.tot.u <- phe.tot[!dup.tindex,]

gen.ids <- read.table("WHI.clean.final.fam",header=FALSE)
names(gen.ids) <- c("SOURCE_SUBJID","IID")
id.key <- read.csv("WHI_ID_KEY.csv",header=TRUE)

id.gen.key <- merge(gen.ids,id.key,by="SOURCE_SUBJID") # ignore warnings, this works

tot.phe.id <- merge(id.gen.key,phe.tot.u,by="dbGaP.SubjID")
tot.annot.phe <- merge(tot.phe.id,annot,by="SOURCE_SUBJID")

MDS <- read.table("whi_principal_components_all_subjects.csv",header=TRUE,sep=",")
names(MDS) <- c('SOURCE_SUBJID','IID', 'C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14','C15','C16','C17','C18','C19','C20')

tot.phe <- merge(tot.annot.phe,MDS,by="SOURCE_SUBJID")

whi.data <- cbind.data.frame(tot.phe$IID.y,tot.phe$SUBJID,tot.phe$sex,tot.phe$AGE,tot.phe$BSP.PlateName,tot.phe[,39:58],rep(4,nrow(tot.phe)),tot.phe$HEIGHTX)
names(whi.data) <- c(names(cardia.full)[1:26],'HEIGHT')

# read in MGS height:

setwd("/home/simonsom/ROH_pathway/TOTAL/MGS/mgs_plink_raw")
mgs.ht <- read.table("mgs.ht",header=TRUE) # height
names(mgs.ht) <- c('db','IID','HEIGHT')
mgs.controls <- read.table("mgs.controls",header=FALSE)
names(mgs.controls) <- c('FID','IID')
mgs.plate <- read.table("plate.covar" ,header=FALSE) # plate
names(mgs.plate) <- c('FID','IID','BATCH')
mgs.phe <- read.csv("phs000021.v2.pht000068.v1.p1.c1.EA_GAIN_phenotype_controls_main.GRU.csv",header=TRUE) 
mgs.age <- cbind.data.frame(mgs.phe$individual_id,mgs.phe$individual_id,mgs.phe$sampling_age) # age
names(mgs.age) <- c('FID','IID','AGE')
mgs.sex <- read.table("sex.covar",header=FALSE) # sex
names(mgs.sex) <- c('FID','IID','SEX')
mgs.eig <- read.table("MGS.eigenvec",header=FALSE)
names(mgs.eig) <- c('FID','IID', 'C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14','C15','C16','C17','C18','C19','C20')

mgs.data.1 <- merge(mgs.ht,mgs.controls,by='IID')[,1:3]
mgs.data.2 <- merge(mgs.data.1,mgs.plate,by='IID')[,c(1,2,3,5)]
mgs.data.3 <- merge(mgs.data.2,mgs.age,by='IID')[,c(1,2,3,4,6)]
mgs.data.4 <- merge(mgs.data.3,mgs.sex,by='IID')[,c(1,2,3,4,5,7)]
mgs.data.5 <- merge(mgs.data.4,mgs.eig,by="IID")
mgs.data <- cbind.data.frame(mgs.data.5,rep(5,nrow(mgs.data.5)))
names(mgs.data) <- c(names(mgs.data[1:27]),'set')

# read in GENEVA height:

setwd("/home/simonsom/ROH_pathway/TOTAL/GENEVA_diab/GEI")
plate.gei.1 <- read.csv("Sample_annotation_consent_0.csv",header=TRUE)
plate.gei.2 <- read.csv("Sample_annotation_consent_1.csv",header=TRUE)
pca.gei <- read.csv("Principal_components_eur_subset.csv",header=TRUE)

setwd("/home/simonsom/ROH_pathway/TOTAL/GENEVA_diab/NHS")
plate.nhs <- read.csv("Sample_annotation.csv",header=TRUE)
pca.nhs <- read.csv("Principal_components_white_subset.csv",header=TRUE)

setwd("/home/simonsom/ROH_pathway/TOTAL/GENEVA_diab/TOTAL")
id.key <- read.csv("ID.key.csv",header=TRUE)
ht.g <-read.table("ht.phe",header=TRUE)
names(ht.g) <- c('FID','sample.num','HT','AGE')
# merge plate and pca from NHS and gei together:

plate.gei <- rbind(plate.gei.1,plate.gei.2)
plate.pca <- merge(plate.gei,pca.gei,by='sample.num')
ppca <- plate.pca[,c(1,2,6,8,22:41)]

plate.pca2 <- merge(plate.nhs,pca.nhs,by='sample.num')
ppca2 <- plate.pca2[,c(1,2,6,8,26:45)]
names(ppca2) <- names(ppca)

ppca.tot <- rbind(ppca,ppca2)

geneva.dat <- merge(ppca.tot,ht.g,by='sample.num')
geneva.dat <- cbind.data.frame(geneva.dat,rep(6,nrow(geneva.dat)))
names(geneva.dat) <- c('IID','GID','SEX','PLATE','C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14','C15','C16','C17','C18','C19','C20','FID','HEIGHT','AGE','set')


# Merge training set and replication set data:

# replication set:

m.g.dat <- cbind.data.frame(geneva.dat$FID,geneva.dat$IID,geneva.dat$SEX,geneva.dat$AGE,geneva.dat$PLATE,geneva.dat$set,geneva.dat$HEIGHT,geneva.dat$C1,geneva.dat$C2,geneva.dat$C3,geneva.dat$C4,geneva.dat$C5,geneva.dat$C6,geneva.dat$C7,geneva.dat$C8,geneva.dat$C9,geneva.dat$C10,geneva.dat$C11,geneva.dat$C12,geneva.dat$C13,geneva.dat$C14,geneva.dat$C15,geneva.dat$C16,geneva.dat$C17,geneva.dat$C18,geneva.dat$C19,geneva.dat$C20)
names(m.g.dat) <- c('FID','IID','SEX','AGE','PLATE','set','HEIGHT','C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14','C15','C16','C17','C18','C19','C20')

mgs.g.dat <- cbind.data.frame(mgs.data$FID,mgs.data$IID,mgs.data$SEX,mgs.data$AGE,mgs.data$BATCH,mgs.data$set,mgs.data$HEIGHT,mgs.data$C1,mgs.data$C2,mgs.data$C3,mgs.data$C4,mgs.data$C5,mgs.data$C6,mgs.data$C7,mgs.data$C8,mgs.data$C9,mgs.data$C10,mgs.data$C11,mgs.data$C12,mgs.data$C13,mgs.data$C14,mgs.data$C15,mgs.data$C16,mgs.data$C17,mgs.data$C18,mgs.data$C19,mgs.data$C20)
names(mgs.g.dat) <- c('FID','IID','SEX','AGE','PLATE','set','HEIGHT','C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14','C15','C16','C17','C18','C19','C20')

rep.dat <- rbind(m.g.dat,mgs.g.dat)
rep.dat$set[which(rep.dat$set==5)] <- -1/2
rep.dat$set[which(rep.dat$set==6)] <- 1/2

# height replication data model:

rep.dat <- rep.dat[(!is.na(rep.dat$IID) & !is.na(rep.dat$HEIGHT) & !is.na(rep.dat$SEX) & !is.na(rep.dat$AGE) & !is.na(rep.dat$PLATE) & !is.na(rep.dat$C1)),]
rep.mod <- lm(as.numeric(HEIGHT) ~ as.numeric(AGE) + as.factor(SEX) + as.numeric(C1) + as.numeric(C2) + as.numeric(C3) + as.numeric(C4) + as.numeric(C5) + as.numeric(C6) + as.numeric(C7) + as.numeric(C8) + as.numeric(C9) + as.numeric(C10) + as.numeric(C11) + as.numeric(C12) + as.numeric(C13) + as.numeric(C14) + as.numeric(C15) + as.numeric(C16) + as.numeric(C17) + as.numeric(C18) + as.numeric(C19) + as.numeric(C20) + set + as.factor(PLATE),data=rep.dat)

rep.t <- resid(rep.mod)
rep.phe.total <- cbind.data.frame(rep.dat$FID,rep.dat$IID,rep.t)
names(rep.phe.total) <- c('FID','IID','HT.res')

rep.fam <- read.table("rep.fam",header=FALSE)
names(rep.fam) <- c('FID','IID')
rep.tot <- merge(rep.phe.total,rep.fam,by='IID')

rep.phe.total <- rep.tot[,c(4,1,3)]
names(rep.phe.total) <- c('FID','IID','HT.res')
setwd("/home/simonsom/ROH_pathway/TOTAL/TOTAL_mapping")
write.table(rep.phe.total,file="rep.set.phe",quote=FALSE,row.names=FALSE,col.names=TRUE)

# training set:
aric.data <- merge(aric.data,FID.IID,by='IID')
m.a.dat <- cbind.data.frame(aric.data$FID.x,aric.data$IID,aric.data$SEX,aric.data$AGE,aric.data$BATCH,aric.data$set,aric.data$HEIGHT*.393701,aric.data$C1,aric.data$C2,aric.data$C3,aric.data$C4,aric.data$C5,aric.data$C6,aric.data$C7,aric.data$C8,aric.data$C9,aric.data$C10,aric.data$C11,aric.data$C12,aric.data$C13,aric.data$C14,aric.data$C15,aric.data$C16,aric.data$C17,aric.data$C18,aric.data$C19,aric.data$C20)
names(m.a.dat) <- c('FID','IID','SEX','AGE','PLATE','set','HEIGHT','C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14','C15','C16','C17','C18','C19','C20')
m.a.dat <- m.a.dat[(!is.na(m.a.dat$IID) & !is.na(m.a.dat$HEIGHT) & !is.na(m.a.dat$SEX) & !is.na(m.a.dat$AGE) & !is.na(m.a.dat$PLATE) & !is.na(m.a.dat$C1)),]

mesa.data <- merge(mesa.data,FID.IID,by='IID')
m.mesa.dat <- cbind.data.frame(mesa.data$FID,mesa.data$IID,mesa.data$SEX,mesa.data$AGE,mesa.data$BATCH,mesa.data$set,mesa.data$HEIGHT,mesa.data$C1,mesa.data$C2,mesa.data$C3,mesa.data$C4,mesa.data$C5,mesa.data$C6,mesa.data$C7,mesa.data$C8,mesa.data$C9,mesa.data$C10,mesa.data$C11,mesa.data$C12,mesa.data$C13,mesa.data$C14,mesa.data$C15,mesa.data$C16,mesa.data$C17,mesa.data$C18,mesa.data$C19,mesa.data$C20)
names(m.mesa.dat) <- c('FID','IID','SEX','AGE','PLATE','set','HEIGHT','C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14','C15','C16','C17','C18','C19','C20')
m.mesa.dat <- m.mesa.dat[(!is.na(m.mesa.dat$IID) & !is.na(m.mesa.dat$HEIGHT) & !is.na(m.mesa.dat$SEX) & !is.na(m.mesa.dat$AGE) & !is.na(m.mesa.dat$PLATE) & !is.na(m.mesa.dat$C1)),]

cardia.data <- merge(cardia.data,FID.IID,by='IID')
m.c.dat <-  cbind.data.frame(cardia.data$FID,cardia.data$IID,cardia.data$SEX,cardia.data$AGE,cardia.data$BATCH,cardia.data$set,cardia.data$HEIGHT,cardia.data$C1,cardia.data$C2,cardia.data$C3,cardia.data$C4,cardia.data$C5,cardia.data$C6,cardia.data$C7,cardia.data$C8,cardia.data$C9,cardia.data$C10,cardia.data$C11,cardia.data$C12,cardia.data$C13,cardia.data$C14,cardia.data$C15,cardia.data$C16,cardia.data$C17,cardia.data$C18,cardia.data$C19,cardia.data$C20)
names(m.c.dat) <- c('FID','IID','SEX','AGE','PLATE','set','HEIGHT','C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14','C15','C16','C17','C18','C19','C20')
m.c.dat <- m.c.dat[(!is.na(m.c.dat$IID) & !is.na(m.c.dat$HEIGHT) & !is.na(m.c.dat$SEX) & !is.na(m.c.dat$AGE) & !is.na(m.c.dat$PLATE) & !is.na(m.c.dat$C1)),]


setwd("/home/simonsom/ROH_pathway/TOTAL/TOTAL_mapping")
disc.fam <- read.table("disc.fam",header=FALSE)
names(disc.fam) <- c('FID','IID')
whi.data<- merge(disc.fam,whi.data,by='IID')
m.w.dat <-  cbind.data.frame(whi.data$FID,whi.data$IID,whi.data$SEX,whi.data$AGE,whi.data$BATCH,whi.data$set,whi.data$HEIGHT*.393701,whi.data$C1,whi.data$C2,whi.data$C3,whi.data$C4,whi.data$C5,whi.data$C6,whi.data$C7,whi.data$C8,whi.data$C9,whi.data$C10,whi.data$C11,whi.data$C12,whi.data$C13,whi.data$C14,whi.data$C15,whi.data$C16,whi.data$C17,whi.data$C18,whi.data$C19,whi.data$C20)
names(m.w.dat) <- c('FID','IID','SEX','AGE','PLATE','set','HEIGHT','C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14','C15','C16','C17','C18','C19','C20')
m.w.dat <- m.w.dat[(!is.na(m.w.dat$IID) & !is.na(m.w.dat$HEIGHT) & !is.na(m.w.dat$SEX) & !is.na(m.w.dat$AGE) & !is.na(m.w.dat$PLATE) & !is.na(m.w.dat$C1)),]

disc.dat1 <- rbind(m.a.dat,m.mesa.dat)
disc.dat2 <- rbind(disc.dat1,m.c.dat)
disc.dat <- as.data.frame(rbind(as.matrix(m.w.dat),as.matrix(disc.dat2)))
disc.dat <- as.data.frame(as.matrix(disc.dat[(!is.na(disc.dat$IID) & !is.na(disc.dat$HEIGHT) & !is.na(disc.dat$SEX) & !is.na(disc.dat$AGE) & !is.na(disc.dat$PLATE) & !is.na(disc.dat$C1) & !is.na(disc.dat$set)),]))

disc.mod <- lm(as.numeric(as.character(HEIGHT)) ~ as.numeric(AGE) + as.factor(SEX) + as.numeric(C1) + as.numeric(C2) + as.numeric(C3) + as.numeric(C4) + as.numeric(C5) + as.numeric(C6) + as.numeric(C7) + as.numeric(C8) + as.numeric(C9) + as.numeric(C10) + as.numeric(C11) + as.numeric(C12) + as.numeric(C13) + as.numeric(C14) + as.numeric(C15) + as.numeric(C16) + as.numeric(C17) + as.numeric(C18) + as.numeric(C19) + as.numeric(C20) + as.factor(set) + as.factor(PLATE),data=disc.dat)

disc.t <- resid(disc.mod)
disc.phe.total <- cbind.data.frame(disc.dat$FID,disc.dat$IID,disc.t)
names(disc.phe.total) <- c('FID','IID','V2')

setwd("/home/simonsom/ROH_pathway/TOTAL/TOTAL_mapping")
write.table(disc.phe.total,file="disc.set.phe",quote=FALSE,row.names=FALSE,col.names=TRUE)





