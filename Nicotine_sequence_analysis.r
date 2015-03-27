# Created by Matthew A. Simonson 3/13/2012:

###############################################################################################################
# Examine if number of variants is predictive of substance abuse phenotype; seperately examine CHRNA4, CHRNB2 #
###############################################################################################################

#########################
# Step 1: Organize Data #
#########################

library('foreign') # load foreign library so SPSS data can be read

# Read in phenotype data:

phe.data <- read.spss("PhenotypeConnection.sav",to.data.frame=TRUE) # some warnings occur, but data is fine

# control is 1, case is 2 (check this)

# Read in polymorphism data:

gen1.data <- read.table("CHRNA4map.txt",header=TRUE,fill=TRUE)

gen2.data <- read.table("CHRNB2map.txt",header=TRUE,fill=TRUE)

# Read in genotype data:

gt1.data <- read.table("CHRNA4genotype.txt",header=TRUE,fill=TRUE)

gt2.data <- read.table("CHRNB2genotype.txt",header=TRUE,fill=TRUE)

# Make map format files:
# Format: 1 row per variant, cols= chrom, variant, intergenic, pos

map1 <- as.data.frame(cbind(as.character(gen1.data$chrom),as.character(gen1.data$Variant),as.character(gen1.data$intergenic),as.character(gen1.data$pos))) # CHRNA4 map
names(map1) <- c('chrom','Variant','intergenic','pos') # assign column names

map2 <- as.data.frame(cbind(as.character(gen2.data$chrom),as.character(gen2.data$Variant),as.character(gen2.data$intergenic),as.character(gen2.data$pos))) # CHRNB2 map
names(map2) <- c('chrom','Variant','intergenic','pos')# assign column names


# Make ped format files:
# Format: 1 row per subject; 3 columns, cols= Seq_ID, control, VariantMAcount ... N variants ; (columns are give the name of the variant)

t.geno1 <- as.data.frame(t(gt1.data)) # transpose genotypes 1

t.geno2 <- as.data.frame(t(gt2.data)) # transpose genotypes 2

ped1 <-  as.data.frame(cbind(as.character(phe.data$Seq_ID),as.character(phe.data$control),t.geno1))
names(ped1) <- c('Seq_ID','control',as.character(map1$Variant)) 

ped2 <-  as.data.frame(cbind(as.character(phe.data$Seq_ID),as.character(phe.data$control),t.geno2))
names(ped2) <- c('Seq_ID','control',as.character(map2$Variant))


##############################################
# Using functions from AssotesteR library:   #
##############################################

install.packages("AssotesteR")
library("AssotesteR")

# Convert data into correct format:

CHRNA4.phe<- as.numeric(ped1[,2])
A4o.index <- which(CHRNA4.phe==2)
CHRNA4.phe[A4o.index] <- 0
CHRNA4.gen <- as.matrix(ped1[,(3:ncol(ped1))])

CHRNB2.phe<- as.numeric(ped2[,2])
B2o.index <- which(CHRNB2.phe==2)
CHRNB2.phe[B2o.index] <- 0
CHRNB2.gen <- as.matrix(ped2[,(3:ncol(ped2))])

##############
# Run tests: #
##############

# ASCORE
CHRNA4asc = ASCORE(CHRNA4.phe, CHRNA4.gen, perm=1000)

CHRNB2asc = ASCORE(CHRNB2.phe, CHRNB2.gen, perm=1000)

save(CHRNA4asc,file="2CHRNA4asc.Rdata")
save(CHRNB2asc,file="2CHRNB2asc.Rdata")

# ASSU
#CHRNA4asu = ASSU(CHRNA4.phe, CHRNA4.gen, perm=1000)

#CHRNB2asu = ASSU(CHRNB2.phe, CHRNB2.gen, perm=1000)

#save(CHRNA4asuo,file="2CHRNA4asuo.Rdata")
#save(CHRNB2asuo,file="2CHRNB2asuo.Rdata")

# ASSUW
#CHRNA4assuw = ASSUW(CHRNA4.phe, CHRNA4.gen, perm=1000)

#CHRNB2assuw = ASSUW(CHRNB2.phe, CHRNB2.gen, perm=1000)

#save(CHRNA4assuw,file="2CHRNA4assuw.Rdata")
#save(CHRNB2assuw,file="2CHRNB2assuw.Rdata")

# ASUM
CHRNA4asum = ASUM(CHRNA4.phe, CHRNA4.gen, perm=1000)

CHRNB2asum = ASUM(CHRNB2.phe, CHRNB2.gen, perm=1000)

save(CHRNA4asum,file="2CHRNA4asum.Rdata")
save(CHRNB2asum,file="2CHRNB2asum.Rdata")

# BST
CHRNA4bst = BST(CHRNA4.phe, CHRNA4.gen, perm=1000)

CHRNB2bst = BST(CHRNB2.phe, CHRNB2.gen, perm=1000)

save(CHRNA4bst,file="2CHRNA4bst.Rdata")
save(CHRNB2bst,file="2CHRNB2bst.Rdata")

# C-ALPHA
CHRNA4calp = CALPHA(CHRNA4.phe, CHRNA4.gen, perm=1000)

CHRNB2calp = CALPHA(CHRNB2.phe, CHRNB2.gen, perm=1000)

save(CHRNA4calp,file="2CHRNA4calp.Rdata")
save(CHRNB2calp,file="2CHRNB2calp.Rdata")

# CARV
CHRNA4carv = CARV(CHRNA4.phe, CHRNA4.gen, waf=TRUE, approach="variable", perm=1000)

CHRNB2carv = CARV(CHRNB2.phe, CHRNB2.gen, waf=TRUE, approach="variable", perm=1000)

save(CHRNA4carv,file="2CHRNA4carv.Rdata")
save(CHRNB2carv,file="2CHRNB2carv.Rdata")

# CAST
#CHRNA4cast = CAST(CHRNA4.phe, CHRNA4.gen, test = "fisher", perm=1000)

#CHRNB2cast = CAST(CHRNB2.phe, CHRNB2.gen, test = "fisher", perm=1000)

#save(CHRNA4cast,file="2CHRNA4cast.Rdata")
#save(CHRNB2cast,file="2CHRNB2cast.Rdata")

# CMAT
CHRNA4cmat = CMAT(CHRNA4.phe, CHRNA4.gen, perm=1000)

CHRNB2cmat = CMAT(CHRNB2.phe, CHRNB2.gen, perm=1000)

save(CHRNA4cmat,file="2CHRNA4cmat.Rdata")
save(CHRNB2cmat,file="2CHRNB2cmat.Rdata")

# CMC
CHRNA4cmc = CMC(CHRNA4.phe, CHRNA4.gen, maf=0.05, perm=1000)

CHRNB2cmc = CMC(CHRNB2.phe, CHRNB2.gen, maf=0.05, perm=1000)

save(CHRNA4cmc,file="2CHRNA4cmc.Rdata")
save(CHRNB2cmc,file="2CHRNB2cmc.Rdata")

# GDBR
CHRNA4gdbr = GDBR(CHRNA4.phe, CHRNA4.gen, distance = "IBS", perm=1000)

CHRNB2gdbr = GDBR(CHRNB2.phe, CHRNB2.gen, distance = "IBS", perm=1000)

save(CHRNA4gdbr,file="2CHRNA4gdbr.Rdata")
save(CHRNB2gdbr,file="2CHRNB2gdbr.Rdata")

# ORWSS
CHRNA4orwss = ORWSS(CHRNA4.phe, CHRNA4.gen, c.param=NULL, perm=1000)

CHRNB2orwss = ORWSS(CHRNB2.phe, CHRNB2.gen, c.param=NULL, perm=1000)

save(CHRNA4orwss,file="2CHRNA4orwss.Rdata")
save(CHRNB2orwss,file="2CHRNB2orwss.Rdata")

# RARECOVER
CHRNA4rcvr = RARECOVER(CHRNA4.phe, CHRNA4.gen, maf=0.05, perm=1000)

CHRNB2rcvr = RARECOVER(CHRNB2.phe, CHRNB2.gen, maf=0.05, perm=1000)

save(CHRNA4rcvr,file="2CHRNA4rcvr.Rdata")
save(CHRNB2rcvr,file="2CHRNB2rcvr.Rdata")

# RBT
CHRNA4rbt = RBT(CHRNA4.phe, CHRNA4.gen, perm=1000)

CHRNB2rbt = RBT(CHRNB2.phe, CHRNB2.gen, perm=1000)

save(CHRNA4rbt,file="2CHRNA4rbt.Rdata")
save(CHRNB2rbt,file="2CHRNB2rbt.Rdata")

# RVT1
CHRNA4rvt1 = RVT1(CHRNA4.phe, CHRNA4.gen, maf=0.01, perm=1000)

CHRNB2rvt1 = RVT1(CHRNB2.phe, CHRNB2.gen, maf=0.01, perm=1000)

save(CHRNA4rvt1,file="2CHRNA4rvt1.Rdata")
save(CHRNB2rvt1,file="2CHRNB2rvt1.Rdata")

# RVT2
CHRNA4rvt2 = RVT2(CHRNA4.phe, CHRNA4.gen, maf=0.01, perm=1000)

CHRNB2rvt2 = RVT2(CHRNB2.phe, CHRNB2.gen, maf=0.01, perm=1000)

save(CHRNA4rvt2,file="2CHRNA4rvt2.Rdata")
save(CHRNB2rvt2,file="2CHRNB2rvt2.Rdata")

# RWAS
CHRNA4rwas = RWAS(CHRNA4.phe, CHRNA4.gen, perm=1000)

CHRNB2rwas = RWAS(CHRNB2.phe, CHRNB2.gen, perm=1000)

save(CHRNA4rwas,file="2CHRNA4rwas.Rdata")
save(CHRNB2rwas,file="2CHRNB2rwas.Rdata")

# SCORE
CHRNA4scor = SCORE(CHRNA4.phe, CHRNA4.gen, perm=1000)

CHRNB2scor = SCORE(CHRNB2.phe, CHRNB2.gen, perm=1000)

save(CHRNA4scor,file="2CHRNA4scor.Rdata")
save(CHRNB2scor,file="2CHRNB2scor.Rdata")

# SEQSUM
CHRNA4seqsm = SEQSUM(CHRNA4.phe, CHRNA4.gen, perm=1000)

CHRNB2seqsm = SEQSUM(CHRNB2.phe, CHRNB2.gen, perm=1000)

save(CHRNA4seqsm,file="2CHRNA4seqsm.Rdata")
save(CHRNB2seqsm,file="2CHRNB2seqsm.Rdata")

# SKAT
CHRNA4sk = SKAT(CHRNA4.phe, CHRNA4.gen, perm=1000)

CHRNB2sk = SKAT(CHRNB2.phe, CHRNB2.gen, perm=1000)

save(CHRNA4sk,file="2CHRNA4sk.Rdata")
save(CHRNB2sk,file="2CHRNB2sk.Rdata")

# SSU
CHRNA4ssu = SSU(CHRNA4.phe, CHRNA4.gen, perm=1000)

CHRNB2ssu = SSU(CHRNB2.phe, CHRNB2.gen, perm=1000)

save(CHRNA4ssu,file="2CHRNA4ssu.Rdata")
save(CHRNB2ssu,file="2CHRNB2ssu.Rdata")

# SSUW
#CHRNA4ssuw = SSUW(CHRNA4.phe, CHRNA4.gen, perm=1000)

#CHRNB2ssuw = SSUW(CHRNB2.phe, CHRNB2.gen, perm=1000)

#save(CHRNA4ssuw,file="2CHRNA4ssuw.Rdata")
#save(CHRNB2ssuw,file="2CHRNB2ssuw.Rdata")

# SUM
#CHRNA4sum = SUM(CHRNA4.phe, CHRNA4.gen, perm=1000)

#CHRNB2sum = SUM(CHRNB2.phe, CHRNB2.gen, perm=1000)

#save(CHRNA4sum,file="2CHRNA4wst.Rdata")
#save(CHRNB2sum,file="2CHRNB2sum.Rdata")

# T-TEST
#CHRNA4t.t = TTEST(CHRNA4.phe, CHRNA4.gen)

#CHRNB2t.t = TTEST(CHRNB2.phe, CHRNB2.gen)

#save(CHRNA4t.t,file="2CHRNA4t.t.Rdata")
#save(CHRNB2t.t,file="2CHRNB2t.t.Rdata")

# UMINP
CHRNA4uminp = UMINP(CHRNA4.phe, CHRNA4.gen, perm=1000)

CHRNB2uminp = UMINP(CHRNB2.phe, CHRNB2.gen, perm=1000)

save(CHRNA4uminp,file="2CHRNA4uminp.Rdata")
save(CHRNB2uminp,file="2CHRNB2uminp.Rdata")

# VT
CHRNA4vt = VT(CHRNA4.phe, CHRNA4.gen, maf=0.01, perm=1000)

CHRNB2vt = VT(CHRNB2.phe, CHRNB2.gen, maf=0.01, perm=1000)

save(CHRNA4vt,file="2CHRNA4vt.Rdata")
save(CHRNB2vt,file="2CHRNB2vt.Rdata")

# WSS
CHRNA4wss = WSS(CHRNA4.phe, CHRNA4.gen, perm=1000)

CHRNB2wss = WSS(CHRNB2.phe, CHRNB2.gen, perm=1000)

save(CHRNA4wss,file="2CHRNA4wss.Rdata")
save(CHRNB2wss,file="2CHRNB2wss.Rdata")

# WST
CHRNA4wst = WST(CHRNA4.phe, CHRNA4.gen, perm=1000)

CHRNB2wst = WST(CHRNB2.phe, CHRNB2.gen, perm=1000)

save(CHRNA4wst,file="2CHRNA4wst.Rdata")
save(CHRNB2wst,file="2CHRNB2wst.Rdata")

# Load saved data:

system("ls *.Rdata > files")
files <- read.table("files",header=FALSE)

for(i in 1:nrow(files)){
  load(as.character(files$V1[i]))
}
