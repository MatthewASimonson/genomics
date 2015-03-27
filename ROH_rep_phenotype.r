###########################################################################
# Phenotype and covariate data: object "full.data" contains all variables #
###########################################################################

load("phe.model.Rdata")
fam.dat <- read.table("MERGE.clean.FINAL.fam",header=FALSE)
FID.IID <- fam.dat[,1:2]
names(FID.IID) <- c('FID','IID')
# remove NA's from total data
total <- total[(!is.na(total$IID) & !is.na(total$BMI) & !is.na(total$SEX) & !is.na(total$AGE) & !is.na(total$BATCH) & !is.na(total$C1)),]
# remove any duplicates:
dup.index <- duplicated(total$IID)
full.data <- total[which(dup.index==FALSE),]


# Read in total phenotypes for ARIC:

phen.total <- read.csv("ARIC_CARE_pheno.GRU.csv",header=TRUE)
phen.full <- merge(phen.total,full.data,by="IID") # full ARIC phenotype data with all covariates

# get residualized phenotype for all variables:

resid.phen.full <- list()

for(i in 225:364){
variable <- phen.full[,i]
var.na.index <- is.na(variable)
full.mod <-  lm(as.numeric(variable[var.na.index==FALSE]) ~ as.numeric(AGE) + as.factor(SEX) + as.numeric(C1) + as.numeric(C2) + as.numeric(C3) + as.numeric(C4) + as.numeric(C5) +as.numeric(C6)+ as.numeric(C7)+ as.numeric(C8)+ as.numeric(C9)+ as.numeric(C10)+ as.numeric(C11)+ as.numeric(C12)+ as.numeric(C13)+ as.numeric(C14)+ as.numeric(C15)+ as.numeric(C16) +as.numeric(C17)+ as.numeric(C18)+as.numeric(C19) +as.numeric(C20)+as.factor(BATCH) + as.factor(centerid),data=phen.full[var.na.index==FALSE,]) # ALSO CONTROL FOR WHICH ARIC FIELD CENTER SUBJECT WAS FROM
resids <- resid(full.mod)
inds <- phen.full[var.na.index==FALSE,1]
resid.pheno <- cbind.data.frame(inds,resids)
names(resid.pheno) <- c('IID','PHE')
resid.phen.full[[i-6]] <- resid.pheno # phenotypes start at 7th column and exlude female only
print(i)
}
names(resid.phen.full) <- names(phen.full)[7:364]

save.image("ROH.Rdata")
load("ROH.Rdata")


# Examine if any phenotypes are significant with ROH burden:

#setwd("/STATGEN/home/simonsom/ROH_pathway/TOTAL/")

# levels of LD-pruning & SNP density:

hom.i.light <- read.table("MERGE.clean.FINALMERGE_ROH_lite_snp65.hom.indiv",header=TRUE)

#hom.i.35mod <- read.table("MERGE.clean.FINALMERGE_ROH_mod_snp35.hom.indiv",header=TRUE)

#hom.i.45mod <- read.table("MERGE.clean.FINALMERGE_ROH_mod_snp45.hom.indiv",header=TRUE)

#hom.i.50mod <- read.table("MERGE.clean.FINALMERGE_ROH_mod_snp50.hom.indiv",header=TRUE)

####### loop through each phenotype:

# examine light pvals
# This is what was found to be optimal ROH detection parameters in previous papers

light.pvals <- vector()
light.betas <- vector()
for( i in 1:length(resid.phen.full)){
  if(is.null(resid.phen.full[[i]])==FALSE){
temp.dat <- merge(resid.phen.full[[i]],hom.i.light,by="IID")
modl <- lm(PHE.x~KB,data=temp.dat)
light.pvals[i] <- summary(modl)$coefficients[2,4] # p-value
light.betas[i] <- summary(modl)$coefficients[2,1]
}else{
  light.pvals[i] <- NA
  light.betas[i] <- NA
}
print(i)
}

light.dat <- cbind.data.frame(light.betas,round(light.pvals,10))
row.names(light.dat) <- names(phen.full)[7:364]
o.index <- order(light.pvals)

light.frame <- light.dat[o.index,]

# examine 35mod pvals

mod35.pvals <- vector()
mod35.betas <- vector()
for( i in 1:length(resid.phen.full)){
  if(is.null(resid.phen.full[[i]])==FALSE){
temp.dat <- merge(resid.phen.full[[i]],hom.i.35mod,by="IID")
modl <- lm(PHE.x~KB,data=temp.dat)
mod35.pvals[i] <- summary(modl)$coefficients[2,4] # p-value
mod35.betas[i] <- summary(modl)$coefficients[2,1]
}else{
  mod35.pvals[i] <- NA
  mod35.betas[i] <- NA
}
print(i)
}

mod35.dat <- cbind.data.frame(mod35.betas,round(mod35.pvals,10))
row.names(mod35.dat) <- phe.labels
o.index <- order(mod35.pvals)

mod35.frame <- mod35.dat[o.index,]

# examine 45mod pvals

mod45.pvals <- vector()
mod45.betas <- vector()
for( i in 1:length(resid.phen.full)){
  if(is.null(resid.phen.full[[i]])==FALSE){
temp.dat <- merge(resid.phen.full[[i]],hom.i.45mod,by="IID")
modl <- lm(PHE.x~KB,data=temp.dat)
mod45.pvals[i] <- summary(modl)$coefficients[2,4] # p-value
mod45.betas[i] <- summary(modl)$coefficients[2,1]
}else{
  mod45.pvals[i] <- NA
  mod45.betas[i] <- NA
}
print(i)
}

mod45.dat <- cbind.data.frame(mod45.betas,round(mod45.pvals,10))
row.names(mod45.dat) <- phe.labels
o.index <- order(mod45.pvals)

mod45.frame <- mod45.dat[o.index,]

# examine 50mod pvals

mod50.pvals <- vector()
mod50.betas <- vector()
for( i in 1:length(resid.phen.full)){
  if(is.null(resid.phen.full[[i]])==FALSE){
temp.dat <- merge(resid.phen.full[[i]],hom.i.50mod,by="IID")
modl <- lm(PHE.x~KB,data=temp.dat)
mod50.pvals[i] <- summary(modl)$coefficients[2,4] # p-value
mod50.betas[i] <- summary(modl)$coefficients[2,1]
}else{
  mod50.pvals[i] <- NA
  mod50.betas[i] <- NA
}
print(i)
}

mod50.dat <- cbind.data.frame(mod50.betas,round(mod50.pvals,10))
row.names(mod50.dat) <- phe.labels
o.index <- order(mod50.pvals)

mod50.frame <- mod50.dat[o.index,]

save.image("ROH.Rdata")

#####################################
# Examine just height ROH in ARIC:  #
#####################################

height.mod <-  lm(as.numeric(phen.full$anta01) ~ as.numeric(AGE) + as.factor(SEX) + as.numeric(C1) + as.numeric(C2) + as.numeric(C3) + as.numeric(C4) + as.numeric(C5) +as.numeric(C6)+ as.numeric(C7)+ as.numeric(C8)+ as.numeric(C9)+ as.numeric(C10)+ as.numeric(C11)+ as.numeric(C12)+ as.numeric(C13)+ as.numeric(C14)+ as.numeric(C15)+ as.numeric(C16) +as.numeric(C17)+ as.numeric(C18)+as.numeric(C19) +as.numeric(C20)+as.factor(BATCH) + as.factor(centerid),data=phen.full) #

ht.resid <- resid(height.mod)
inds <- phen.full[,1]
resid.pheno <- cbind.data.frame(inds,ht.resid)
names(resid.pheno) <- c('IID','HT')
hom.i.light$logKB <- log10(hom.i.light$KB)
ht.dat <- merge(resid.pheno,hom.i.light,by="IID")
ht.dat <- ht.dat[which(ht.dat$KB!=0),]
ht.mod <- lm(HT~KB,data=ht.dat)

#####################
# WHI phenotype:    #
#####################

# read in phenotype and covariate data:

fam.data <- read.table("WHI.clean.final.fam",header=FALSE)
names(fam.data) <- c('ID1','ShareID')

whi.phe1 <- read.csv("WHI.phe1.GRU.csv",header=TRUE)
whi.phe2 <- read.csv("WHI.phe2.GRU.csv",header=TRUE)
pcas <- read.csv("whi_principal_components_all_subjects.csv",header=TRUE)
annot <- read.csv("Sample_annotation.csv",header=TRUE)
ID.key <- read.csv("ID_key.csv",header=TRUE)

phe.id1 <- whi.phe1[,1:2]
names(phe.id1) <- c('V1','V2')

phe.id2 <- whi.phe2[,1:2]
names(phe.id2) <- c('V1','V2')

phe.id.tot <- rbind(phe.id1,phe.id2)
names(phe.id.tot) <- c('ID2','ShareID')
test <- merge(fam.data,phe.id.tot,by="ShareID")

###########################
# CARDIA                  #
###########################
load("phe.model.Rdata")
fam.dat <- read.table("MERGE.clean.FINAL.fam",header=FALSE)
FID.IID <- fam.dat[,1:2]
names(FID.IID) <- c('FID','IID')
# remove NA's from total data
total <- total[(!is.na(total$IID) & !is.na(total$BMI) & !is.na(total$SEX) & !is.na(total$AGE) & !is.na(total$BATCH) & !is.na(total$C1)),]
# remove any duplicates:
dup.index <- duplicated(total$IID)
full.data <- total[which(dup.index==FALSE),]

# Read in total phenotypes for CARDIA:

phen.card <- read.csv("CARDIA.pheno.csv",header=TRUE)
card.ht <- cbind.data.frame(phen.card[,2],phen.card$A20HGT)
names(card.ht) <- c('IID','HT')
phen.full <- merge(card.ht,full.data,by="IID") 

#######################################
# Examine just height ROH in CARDIA:  #
#######################################
hom.i.light <- read.table("MERGE.clean.FINALMERGE_ROH_lite_snp65.hom.indiv",header=TRUE)

height.mod <-  lm(HT ~ as.numeric(AGE) + as.factor(SEX) + as.numeric(C1) + as.numeric(C2) + as.numeric(C3) + as.numeric(C4) + as.numeric(C5) +as.numeric(C6)+ as.numeric(C7)+ as.numeric(C8)+ as.numeric(C9)+ as.numeric(C10)+ as.numeric(C11)+ as.numeric(C12)+ as.numeric(C13)+ as.numeric(C14)+ as.numeric(C15)+ as.numeric(C16) +as.numeric(C17)+ as.numeric(C18)+as.numeric(C19) +as.numeric(C20)+as.factor(BATCH),data=phen.full) #

ht.resid <- resid(height.mod)
inds <- phen.full[,1]
resid.pheno <- cbind.data.frame(inds,ht.resid)
names(resid.pheno) <- c('IID','HT')
ht.dat <- merge(resid.pheno,hom.i.light,by="IID")
ht.mod <- lm(HT~KB,data=ht.dat)
