
setwd("")

dataset.bfile <- 'ARIC2'

##########################################
#****************************************#
##########################################
# QC for individuals:

######################################################################
# Step 1: Remove individuals who are missing more than 5% of SNP calls
######################################################################
# First step is generating missingness data with PLINK

system(paste("plink --bfile ",dataset.bfile," --missing",sep=""))

ind.miss <- read.table("plink.imiss",header=TRUE)
above.five <- which(ind.miss[,6]>.05)
high.missing <- ind.miss[above.five,]

# Create a text file with IID and FID for individuals with missingness above 5%

miss.dat <- cbind(as.character(high.missing$FID),as.character(high.missing$IID))

write.table(miss.dat,file="remove.miss.list.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

system(paste("plink --bfile ",dataset.bfile," --remove remove.miss.list.txt --make-bed --out 2.step.miss",sep=""))

#######################################################
# Step 2:  Generate Covariate File ('.covar' extension)
#######################################################

# read in covariate data:

annot <- read.table("Sample_annotation_consent_1.txt",header=TRUE)

covar.file <- cbind(as.character(annot$sample.num),as.character(annot$sex),as.character(annot$V1AGE01),as.character(annot$plate))

# write out covariate files:

write.table(covar.file,file="batch.covar",quote=FALSE,row.names=FALSE,col.names=FALSE)

######################################################################
# Step 3: Prune data outliers with respect to estimated heterozygosity
######################################################################

system("plink --bfile 2.step.miss --het --out 3.step.het")

hetchk <- read.table("3.step.het.het",header=TRUE)

hist(hetchk$F[hetchk$F<.5], breaks=100) # plot distribution of heterozygosity

Three.sd <- 3*(sd(hetchk$F)) # 3 sd's 

drop.het.index <- c(which(hetchk$F>Three.sd), which(hetchk$F<(-1*Three.sd))) # create index of individuals outside 3sd of mean het
# read in het file and prune individuals outside of 3 standard deviations from the mean

drop.het.list <- hetchk[drop.het.index,1:2]

write.table(drop.het.list,file="remove.het.list.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

system("plink --bfile 2.step.miss --remove remove.het.list.txt --make-bed --out 3.step.drop.het")

##################################################################################
# Step 4: Remove individuals with discrepencies between reported and genotypic sex
##################################################################################

system("plink --bfile 3.step.drop.het --check-sex --out 4.step.sex")

sexchk <- read.table("4.step.sex.sexcheck",header=TRUE)
hist(sexchk$F[sexchk$F<.5], breaks=100)

problem.sex <-  sexchk[which(sexchk$STATUS!='OK'),]

FID.IID <- cbind(as.character(problem.sex$FID),as.character(problem.sex$IID))

write.table(FID.IID,file="remove.sex.list.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

# use system command to remove discrepant individuals from data set and write new data files

system("plink --bfile 3.step.drop.het --remove remove.sex.list.txt --make-bed --out 4.step.drop.sex")

########################################################
# Step 5: Prune data for very closely related indviduals
########################################################

#Scan for any individuals with high PIHAT values (e.g. greater than 0.125 (less related than 2nd degree relatives (3rd mean =.125, 2nd mean=.25 ) # USE GCTA for this:

system ("nohup gcta --bfile 4.step.drop.sex --autosome --make-grm --out 5.step.a &")
# START HERE
system ("nohup gcta --grm 5.step.a --grm-cutoff 0.125 --make-grm --out 5.step &")

FID.IID <- read.table("5.step.grm.id",header=FALSE) # check the number of subjects being kept in data set

# use system command to remove discrepant individuals from data set and write new data files
#
system("plink --bfile 4.step.drop.sex --make-founders --keep 5.step.grm.id --make-bed --out 5.step.unrelated")
                    
#######################################################
# Step 6: Multidimensional scaling analysis with HapMap
#######################################################
# NEXT STEP:
# merge dataset with HapMap, THEN prune for LE: 
## convert SNP_A# to rs# so SNPs data can be merged:

SNPlist <- read.table("5.step.unrelated.bim",header=FALSE)# read in list of all SNPs from Obesity1 data

# Read in file with rs# (marker-info file)SNPs:

SNPmi <- read.table("GenomeWideSNP_6.na27.annot.csv",header=TRUE,sep=",")

markers <- SNPmi[,c(1,3)]

# write out SNP names file:

write.table(markers,file="updateSNP.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

# Use --update-name option in plink under data management to rename SNPs so they have rs#'s

system("plink --bfile 5.step.unrelated --update-map updateSNP.txt --update-name --make-bed --out 6.step.update")

# check the new SNP names:

SNPlist <- read.table("6.step.update.bim",header=FALSE)

##################################################################################
# Now merge Data and Hapmap founders:

###Part A: Select 60 unrelated individuals from CEU, YRI, JPT/CHB

# read in fam files from filtered HapMap data on all 3 groups, these are already only founders, thus all unrelated.

CEU.ind <- read.table("hapmap_CEU_r23a_filtered.fam",header=FALSE) 
YRI.ind <- read.table("hapmap_YRI_r23a_filtered.fam",header=FALSE) 
JPT.CHB.ind <- read.table("hapmap_JPT_CHB_r23a_filtered.fam",header=FALSE)
index <- sample(1:90,60)
JPT.CHB.60ind <- JPT.CHB.ind[index,]

founder.list <- rbind(as.matrix(CEU.ind[,1:2]),as.matrix(YRI.ind[,1:2]),as.matrix(JPT.CHB.60ind[,1:2])) # combined list of founder individuals from each group

write.table(founder.list,file="hapmap.founders.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

###Part B: Remove non-founders from hapmap data:
# only keep SNPs included in Obesity1

SNPlist <- read.table("6.step.update.bim",header=FALSE) # read in list of Obesity1 SNPs
write.table(SNPlist[,2],file="Obesity1snp.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

system("plink --bfile hapmap_r23a --keep hapmap.founders.txt --extract Obesity1snp.txt --make-bed --out hapmap_f")

# Also for CEU hapmap:

system("plink --bfile hapmap_CEU_r23a_filtered --keep hapmap.founders.txt --extract Obesity1snp.txt --make-bed --out hapmapCEU_f")

###Part C: Merge hapmap founders with Obesity1 data and resolve strand flip and positional issues:

# Flip strands based on affy reference file:

SNPmi <- read.table("GenomeWideSNP_6.na27.annot.csv",header=TRUE,sep=",")
order.SNP.index <- order(as.numeric(as.character(SNPmi$Chromosome)),as.numeric(as.character(SNPmi$Physical.Position))) # order SNPs by chromosome and physical position
SNPmi <- SNPmi[order.SNP.index,]

rs.strands <- SNPmi[,c(3,6)]
flip.index <- which(rs.strands$Strand=='-')
flip.snps <- rs.strands$dbSNP.RS.ID[flip.index]
write.table(flip.snps,file="flip.list",quote=FALSE,row.names=FALSE,col.names=FALSE)

system("plink --bfile 6.step.update --flip flip.list --recode --make-bed --out 6.step.flip")

# get MAF from each CEU HapMap data and Obesity1 and drop any outliers as a final check:

system("plink --bfile hapmapCEU_f --freq --out hapMAF")

system("plink --bfile 6.step.flip --freq --out datasetMAF")

# read if MAF data and remove duplicate SNPs if any exist:
# 
datasetMAF <- read.table("datasetMAF.frq",header=TRUE)
hapMAF <- read.table("hapMAF.frq",header=TRUE)

dataset.map <- read.table("6.step.flip.bim",header=FALSE)
dataset.snp_pos <- paste(dataset.map$V2,dataset.map$V4,sep=".") # unique identifier for rs# and location

hap.map <- read.table("hapmapCEU_f.bim",header=FALSE)
hap.snp_pos <- paste(hap.map$V2,hap.map$V4,sep=".") # unique identifier for rs# and location

datasetMAF2 <- cbind(datasetMAF,dataset.snp_pos)
hapMAF2 <- cbind(hapMAF,hap.snp_pos)

names(datasetMAF2) <- c('CHR','SNP','A1','A2','MAF','NCHROBS','snp_pos')
names(hapMAF2) <- c('CHR','SNP','A1','A2','MAF','NCHROBS','snp_pos')

merged.MAF <- merge(hapMAF2,datasetMAF2,"snp_pos",all.x=TRUE,all.y=TRUE) # merge MAF for data sets

# find indeces of duplicate SNPs:

y.drops <- which(is.na(merged.MAF$CHR.y))
x.drops <- which(is.na(merged.MAF$CHR.x))

dup.index <- c(x.drops,y.drops)

dup.list <- c(as.character(merged.MAF$SNP.x[y.drops]),as.character(merged.MAF$SNP.y[x.drops])) # list of SNPs that don't have compliments in other dataset

mergedMAF <- merged.MAF[c(-dup.index),] # NEW OBJECT: mergedMAF created DIFF THAN merged.maf

# Find SNPs where marker chromosome is uncertain:

chr.drop.index <- which(merged.MAF$CHR.x!=merged.MAF$CHR.y)
chr.drop.list <- unique(c(as.character(merged.MAF$SNP.y[chr.drop.index]),as.character(merged.MAF$SNP.x[chr.drop.index]))) # list of SNPs where marker chromosome is uncertain

# Examine Allele frequency differences between sets:

MAFdiff <- abs((mergedMAF$MAF.x)-(mergedMAF$MAF.y)) # abs value of MAF differences 

MAFwDIFF <- cbind(mergedMAF,MAFdiff)# merge data

drop.index <- which((MAFwDIFF$MAFdiff>(.5*sd(MAFwDIFF$MAFdiff)))) # create drop index

droplist <- c(as.character(MAFwDIFF$SNP.y[drop.index]),as.character(dup.list),as.character(chr.drop.list))
#
write.table(droplist,file="dropissues.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
#
system("plink --bfile hapmap_f --exclude dropissues.txt --make-bed --out hapmap_f2") # drop SNPs likely to have issues

system("plink --bfile 6.step.flip --exclude dropissues.txt --make-bed --out 6.step.flip2") # drop SNPs likely to have issues

# Now merge the data (or detect residual strand flip problems):
system("plink --bfile hapmap_f2 --bmerge 6.step.flip2.bed 6.step.flip2.bim 6.step.flip2.fam --make-bed --out Obesity1.HM.merge")
#
strand.prob <- read.table("Obesity1.HM.merge.missnp",header=FALSE)

system("plink --bfile hapmap_f2 --exclude Obesity1.HM.merge.missnp --make-bed --out hapmap_f3") # drop detected strand issues

system("plink --bfile 6.step.flip2 --exclude Obesity1.HM.merge.missnp --make-bed --out 6.step.flip3") # drop detected strand issues

# Now merge the data Final time:
system("plink --bfile hapmap_f3 --bmerge 6.step.flip3.bed 6.step.flip3.bim 6.step.flip3.fam --make-bed --out Obesity1.HM.merge")

# Get list of SNPs in linkage equilibrium and create independent SNP data set
#############################################################################

system("nohup plink --bfile Obesity1.HM.merge --indep 200 5 2 &") 
# 
system("nohup plink --bfile Obesity1.HM.merge --extract plink.prune.in --make-bed --out 6.step.LE &")

system("plink --bfile 6.step.LE --thin 0.5 --recode --out 7_step_shellfish") # make sure fewer than 40k SNPs to avoid crash (number of SNPs read into shellfish gets confused when too large for some reason)

# NOTE: MUST CONVERT FROM BINARY TO PED/MAP B4 RUNNING GTOOL 
system("gtool -P --ped 7_step_shellfish.ped --map 7_step_shellfish.map --og shellfish_MDS.gen --os shellfish_MDS.sample")

# Only SNPs in LE were used to generate the 20 principle dimensions to be used as covariates
# Generate IBS matrix of genetic distance to be used for covariates:
# NOTE: gtool format files expected as input
system("nohup shellfish.py --pca --numpcs 20 --maxprocs 17 --file shellfish_MDS --out 7_shell_mds &") ##

# TRANSPOSE SHELLFISH OUTPUT USING BASH SCRIPT 

system({"awk -F ' ' '{
for (f = 1; f <= NF; f++)
a[NR, f] = $f
}
NF > nf { nf = NF }
END {
for (f = 1; f <= nf; f++)
for (r = 1; r <= NR; r++)
printf a[r, f] (r==NR ? RS : FS) 
}' 7_shell_mds.evecs  > 7_shell_mds.tevecs"})

# Read in PCA data:

PCA <- read.table("7_shell_mds.tevecs",header=FALSE)
names(PCA) <- c('C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14','C15','C16','C17','C18','C19','C20')

# Plot MDS data:

sample.index <- 181:nrow(MDS)
JPT.index <- 1:60
YRI.index <- 121:180
CEU.index <- 61:120

# Plot MDS dimensions 1 & 2 to show genetic relatedness of sample subjects relative to Hapmap
plot(MDS$C1[c(sample.index,YRI.index,JPT.index)], MDS$C2[c(sample.index,YRI.index,JPT.index)],main="MDS Genetic Relatedness of Individuals",xlab="Dimension1",ylab="Dimension2")
points(MDS$C1[sample.index],MDS$C2[sample.index],col="blue") # sample points are blue
#points(MDS$C1[CEU.index],MDS$C2[CEU.index],col="green") # CEU Points are green
points(MDS$C1[YRI.index],MDS$C2[YRI.index],col="orange") # YRI points are orange
points(MDS$C1[JPT.index],MDS$C2[JPT.index],col="purple") # JPT/CHB points are red
legend(mean(MDS$C1[YRI.index]),mean(MDS$C2[sample.index]),c('Sample (European Ancestry)','HapMap YRI','HapMap JPT+CHB'),text.col=c('blue','orange','purple'))
#
center.cau.x <- mean(c(MDS$C1[sample.index]))
center.cau.y <- mean(c(MDS$C2[sample.index]))
center.jptc.x <- mean(MDS$C1[JPT.index])
center.jptc.y <- mean(MDS$C2[JPT.index])
center.afc.x <- mean(MDS$C1[YRI.index])
center.afc.y <- mean(MDS$C2[YRI.index])  
segments(center.cau.x,center.cau.y,center.jptc.x,center.jptc.y,col='red')
segments(center.cau.x,center.cau.y,center.afc.x,center.afc.y,col='red')
SA.dist <- sqrt((center.afc.x-center.cau.x)^2)+((center.afc.y-center.cau.y)^2)
SJC.dist <- sqrt((center.jptc.x-center.cau.x)^2)+((center.jptc.y-center.cau.y)^2)

# Now add 2 circles to represent the boundry for genetic dissimilarity; both have a radius of 10% line distance between clusters. Any points outside either boundry are excluded
symbols(center.cau.x,center.cau.y,circles=.1*SA.dist,inches=FALSE, add=TRUE)
symbols(center.cau.x,center.cau.y,circles=.1*SJC.dist,inches=FALSE,add=TRUE)


# NOTE: THIS NEXT SECTION USES THE "fields" LIBRARY:
# install.packages(fields)
library(fields)

# Use distance function to compute distance of each individual from specified (caucasian) centroid:                                        
center.cau <- as.matrix(cbind(center.cau.x,center.cau.y)) # caucasian center coordinates
IND.coord <- as.matrix(cbind(MDS$C1,MDS$C2))
IND.dat <- read.table("6.step.LE.fam",header=FALSE) # read and merge IID's
distance <- cbind((as.vector(rdist(center.cau,IND.coord))),IND.dat[,1:2]) # create a vector of the distance of each individual away from the centroid of the sample cluster merged with FID & IID
names(distance) <- c('distance','FID','IID')

# generate list of individuals who are inside the inner radius to keep:
keep.list <- distance[which(distance$distance<(.1*SJC.dist)),2:3]

write.table(keep.list,file="keepstrat.list.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

##################################################################################
# Step 7: Prune data for Stratification & make sure all SNPs are in LE
##################################################################################

system("plink --bfile 6.step.LE --keep keep.strat.list.txt --make-bed --out 7.step.LEstrat")
system("plink --bfile 5.step.unrelated --keep keep.strat.list.txt --make-bed --out 7.step.LDstrat") # remove individuals from full LD data set that don't pass stratification

##########################################
#****************************************#
##########################################
# QC for SNPs:
#################################
# Step 8: Drop SNPs with MAF <.01
#################################

system("plink --bfile 7.step.LDstrat --maf 0.01 --make-bed --out 8.step.noMAF05")
# 
#########################################
# Step 9: Drops SNPs with call rates <.05
#########################################

system("plink --bfile 8.step.noMAF05 --missing")
snp.miss <- read.table("plink.lmiss",header=TRUE)
bad.fmiss.index <- which(snp.miss$F_MISS > .05) # greater than 5% missingness removed
bad.snp.list <- snp.miss[bad.fmiss.index,2] # create a list a SNPs with call rates <=.95
write.table(bad.snp.list,file="remove.highmiss.list.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

system("plink --bfile 8.step.noMAF05 --exclude remove.highmiss.list.txt --make-bed --out 9.step.snpcall")

################################
# Step 10: Prune SNPs out of HWE
################################

system("plink --bfile 9.step.snpcall --hardy")

HWE.data <- read.table("plink.hwe",header=TRUE)
out.HWE <- which(HWE.data[,9]<.001)# 10^-3
prune.HWE.snps <- HWE.data[out.HWE,2]

write.table(prune.HWE.snps,file="remove.HWE.list.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
#
system("plink --bfile 9.step.snpcall --exclude remove.HWE.list.txt --make-bed --out 10.step.HWE")

###########################################################################################
# Step 11: non-random genotyping failure, as inferred by the flanking haplotypic background 
###########################################################################################

# drop SNPs with P < 10^-10
system("plink --bfile 10.step.HWE --test-mishap") # generate mishap file

mishap <- read.table("plink.missing.hap",header=TRUE)

fail.index <- which(mishap$P<(10^-10))

drop.hap<- as.data.frame(mishap$SNP[fail.index])

write.table(drop.hap,file="remove.hap.list.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

system("plink --bfile 10.step.HWE --exclude remove.hap.list.txt --make-bed --out 11.step.hap")

#####################################################################
########## ***** END OF ALL DATA CLEANING PROCEDUERS ***** ##########
#####################################################################
