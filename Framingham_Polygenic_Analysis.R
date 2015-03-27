
setwd("/STATGEN/home/shared/km/Framingham/Merged_data")

##########################################
#****************************************#
##########################################
# QC for SNPs:

############################################################################################
# Step 1: Make sure physical positions of SNPs is correct and make all individuals founders
############################################################################################

#Check position of SNPs in Illumina annotation file and update position of all SNPs

# read in annotated affymetrix data and compare with previous SNP data:
#part 1: Sort data by chromosome

# Affy 500k contains 2 annotation files:

affy.annot1 <- read.table("Mapping250K_Nsp.na30.annot.csv", header=TRUE,sep=",", colClasses=c(rep("character",5),rep("NULL",22)))

affy.annot2 <- read.table("Mapping250K_Sty.na31.annot.csv", header=TRUE,sep=",", colClasses=c(rep("character",5),rep("NULL",22)))

# Merge 2 affy 500k annotation files:

affy.annot <- rbind(affy.annot1,affy.annot2) # merge 2 files

sorted.affy.index <-order(affy.annot$Chromosome,affy.annot$Physical.Position)

sorted.affy <- affy.annot[sorted.affy.index,] 

# remove SNPs with no chromosome or physical marker info
keep.rows <- which(sorted.affy$Chromosome!='---')
sorted.affy <- sorted.affy[keep.rows,1:5]# remove after column 7
affy.snps <- substr(sorted.affy[,2],start=2,stop=10) ## USED SUBSTR HERE BECAUSE OF WEIRD FRAMINGHAM SNP FORMAT
# read in Map file from original data

map <- read.table("merge.NHLBI.bim", header=FALSE)
snps <- as.vector(substr(map[,2],start=2,stop=10))

# match snps in old data with affy list

snps.2.compare <- match(snps,affy.snps)
relevent.affy <- sorted.affy[snps.2.compare,] # only relevent snps
diff.position <- which(map[,4!=relevent.affy[,5],]) # indeces of snps with changed base postions
snp.diff <- abs(as.numeric(relevent.affy[diff.position,5])-map[diff.position,4]) # how large is the difference
large.diff.index <- which(snp.diff>100)
remove.snps <- as.data.frame(map[diff.position,2][large.diff.index])
write.table(remove.snps,file="diff.position.txt", quote=FALSE,sep="",row.names=FALSE,col.names=FALSE) # write list of snps to remove that have significantly different base positions

# Write new data file that excludes specified SNPs
system("plink --bfile merge.NHLBI --make-founders --exclude diff.position.txt --make-bed --out 2.step.annot")

##################################################################################
# Step 2: Drop SNPs with MAF <.01
##################################################################################

system("plink --bfile 2.step.annot --maf 0.01 --write-snplist --make-bed --out 3.step.noMAF05")

# check and see what MAF currently is in data after pruning:

system("plink --bfile 3.step.noMAF05 --freq --out MAFdata")

##################################################################################
# Step 3: Drops SNPs with call rates <.05
##################################################################################

system("plink --bfile 3.step.noMAF05 --missing")
snp.miss <- read.table("plink.lmiss",header=TRUE)
bad.fmiss.index <- which(snp.miss$F_MISS > .05) # greater than 5% missingness removed
bad.snp.list <- snp.miss[bad.fmiss.index,2] # create a list a SNPs with call rates <=.95
write.table(bad.snp.list,file="remove.highmiss.list.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

system("plink --bfile 3.step.noMAF05 --exclude remove.highmiss.list.txt --make-bed --out 4.step.snpcall")

##################################################################################
# Step 4:  Prune SNPs out of HWE
##################################################################################
system("plink --bfile 4.step.snpcall --hardy")

HWE.data <- read.table("plink.hwe",header=TRUE)
out.HWE <- which(HWE.data[,9]<.001)# 10^-3
prune.HWE.snps <- HWE.data[out.HWE,2]

write.table(prune.HWE.snps,file="remove.HWE.list.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
#
system("plink --bfile 4.step.snpcall --exclude remove.HWE.list.txt --make-bed --out 5.step.HWE")

##################################################################################
# Step 5: non-random genotyping failure, as inferred by the flanking haplotypic background 
##################################################################################

# drop SNPs with P < 10^-10
system("plink --bfile 5.step.HWE --test-mishap") # generate mishap file

mishap <- read.table("plink.missing.hap",header=TRUE)

fail.index <- which(mishap$P<(10^-10))

drop.hap<- as.data.frame(mishap$SNP[fail.index])

write.table(drop.hap,file="remove.hap.list.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

system("plink --bfile 5.step.HWE --exclude remove.hap.list.txt --make-bed --out 6.step.hap")

##########################################
#****************************************#
##########################################
# QC for individuals:


#################################################################################
# Step 6: Remove individuals who are missing more than 5% of SNP calls
#################################################################################
# First step is generating missingness data with PLINK

system("plink --bfile 6.step.hap --missing")

ind.miss <- read.table("plink.imiss",header=TRUE)
above.five <- which(ind.miss[,6]>.05)
high.missing <- ind.miss[above.five,]

# Create a text file with IID and FID for individuals with missingness above 5%

miss.dat <- cbind(as.character(high.missing$FID),as.character(high.missing$IID))

write.table(miss.dat,file="remove.miss.list.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

system("plink --bfile 6.step.hap --remove remove.miss.list.txt --make-bed --out 7.step.miss")

#################################################################################################################################################################
# Step 7:  evidence of gross non-random plate failure: Use Batch effects inferred from missingness to generate batch covariate to be used in association
################################################################################################################################################################

# read in phenotype data:
risk.scores <- read.table("Risk_Scores.csv",header=FALSE,skip=1,sep=",")
# find indeces of individuals that have a risk score
notna.index <- as.logical((is.na(risk.scores$V6)*-1)+1)
risk.data<- risk.scores[notna.index,] # only individuals with risk scores are included in the data now
names(risk.data)<- c("IID")
merged.data <- merge(risk.data,fam.data,"IID") # merge fam data and risk data

pheno.file <- as.data.frame(cbind(as.character(merged.data[,7]),as.character(merged.data[,1]),as.character(merged.data[,6])))
names(pheno.file) <- c("FID","IID","PHE")
# write risk score data:

write.table(pheno.file,file="risk.score.batch",quote=FALSE,row.names=FALSE,col.names=FALSE)

###########################################
# Generate phenotype data for plink file:
system("nohup plink --bfile 7.step.miss --cluster-missing --out miss.cluster &")

fam.data <- read.table("7.step.miss.fam",header=FALSE)
names(fam.data)<- c("FID","IID")
#### Now for covariate:

# Read in IBM matrix:
IBM <- read.table("miss.cluster.mdist.missing",header=FALSE) # read in individual missingness data

# Use scree plot to determine correct number of eigenvectors to use as covariates

# Now convert data to eigenvectors
dist <- matrix(NA,nrow=nrow(IBM),ncol=ncol(IBM))
dist <- 1-IBM
fit <- cmdscale(dist,eig=TRUE, k=2) # k is the number of dim
# compare predictive effect of dim1 and dim2 to individual missingness rate

covar.fam <- as.data.frame(cbind(as.character(fam.data$FID),as.character(fam.data$IID),as.numeric(as.character(fit$points[,1])),as.numeric(as.character(fit$points[,2]))))
covar.fam[,3] <- as.numeric(as.character(fit$points[,1]))
covar.fam[,4] <- as.numeric(as.character(fit$points[,2]))

covar.file <- covar.fam

# write out covariate files:

write.table(covar.file,file="batch.covar",quote=FALSE,row.names=FALSE,col.names=FALSE)

##################################################################################
# Step 8: Prune data outliers with respect to estimated heterozygosity
##################################################################################

system("plink --bfile 7.step.miss --het --out 8.step.het")
hetchk <- read.table("8.step.het.het",header=TRUE)

hist(hetchk$F[hetchk$F<.5], breaks=100) # plot distribution of heterozygosity

Three.sd <- 3*(sd(hetchk$F)) # 3 sd's 

drop.het.index <- c(which(hetchk$F>Three.sd), which(hetchk$F<(-1*Three.sd))) # create index of individuals outside 3sd of mean het
# read in het file and prune individuals outside of 3 standard deviations from the mean

drop.het.list <- hetchk[drop.het.index,1:2]

write.table(drop.het.list,file="remove.het.list.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

system("plink --bfile 7.step.miss --remove remove.het.list.txt --make-bed --out 8.step.drop.het")

#################################################################################
# Step 9: Remove individuals with discrepencies between reported and genotypic sex
#################################################################################

system("plink --bfile 8.step.drop.het --check-sex --out 9.step.sex")
sexchk <- read.table("9.step.sex.sexcheck",header=TRUE)
hist(sexchk$F[sexchk$F<.5], breaks=100)

problem.sex <-  sexchk[which(sexchk$STATUS!='OK'),]

FID.IID <- cbind(as.character(problem.sex$FID),as.character(problem.sex$IID))

write.table(FID.IID,file="remove.sex.list.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

# use system command to remove discrepant individuals from data set and write new data files

system("plink --bfile 8.step.drop.het --remove remove.sex.list.txt --make-bed --out 9.step.drop.sex")

##################################################################################
# Step 10: Prune data for very closely related indviduals
##################################################################################

#Scan the plink.genome file for any individuals with high PIHAT values (e.g. greater than 0.125 (less related than 2nd degree relatives (3rd mean =.125, 2nd mean=.25 )

system("plink --bfile 9.step.drop.sex --genome --min 0.05 --out 10.close.related") # Create list of closely related individuals

relchk <- read.table("10.close.related.genome",header=TRUE)

d <- which(relchk$PI_HAT>.125)

FID.IID <- as.data.frame(cbind(as.character(relchk$FID1[d]),as.character(relchk$IID1[d])))

write.table(FID.IID,file="remove.rel.list.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

# use system command to remove discrepant individuals from data set and write new data files
#
system("plink --bfile 9.step.drop.sex --make-founders --remove remove.rel.list.txt --make-bed --out 10.step.unrelated")
                    
##################################################################################
# Step 11 Multidimensional scaling analysis with HapMap
##################################################################################
# NEXT STEP:
# merge Framingham with HapMap, THEN prune for LE: 
## convert ss# to rs# so SNPs data can be merged:

SNPlist <- read.table("10.step.unrelated.bim",header=FALSE)# read in list of all SNPs from FHS data

# Read in file with rs# (marker-info file)SNPs:

SNPmi <- read.table("FHSmarker-info.csv",header=FALSE,sep=",")

markers <- SNPmi[,4:5]

# write out SNP names file:

write.table(markers,file="updateSNP.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

# Use --update-name option in plink under data management to rename SNPs so they have rs#'s

system("plink --bfile 10.step.unrelated --update-map updateSNP.txt --update-name --make-bed --out 11.step.update")

# check the new SNP names:

SNPlist <- read.table("11.step.update.bim",header=FALSE)

##################################################################################
# Now merge FHS and Hapmap founders:

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
# only keep SNPs included in FHS

SNPlist <- read.table("11.step.update.bim",header=FALSE) # read in list of FHS SNPs
write.table(SNPlist[,2],file="FHSsnp.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

system("plink --bfile hapmap_r23a --keep hapmap.founders.txt --extract FHSsnp.txt --make-bed --out hapmap_f")

# Also for CEU hapmap:

system("plink --bfile hapmap_CEU_r23a_filtered --keep hapmap.founders.txt --extract FHSsnp.txt --make-bed --out hapmapCEU_f")

###Part C: Merge hapmap founders with FHS data and resolve strand flip and positional issues:

# first get list of detected strand and positional issues from plink

system("plink --bfile hapmap_f --bmerge 11.step.update.bed 11.step.update.bim 11.step.update.fam --make-bed --out FHS.HM.merge")

# read in detected strand flip issues:

strand <- read.table("FHS.HM.merge.missnp",header=FALSE)
posit <- read.table("DiffPos.txt",header=FALSE,fill=TRUE)

droplist <- rbind(strand,posit)

# write list of SNPs to drop

write.table(droplist,file="dropissues.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

# remove bad SNPs from FHS data:

system("plink --bfile 11.step.update --exclude dropissues.txt --make-bed --out 11.step.update")

# remove bad SNPs from HapMap:

system("plink --bfile hapmap_f --exclude dropissues.txt --make-bed --out hapmap_f")

# get MAF from each CEU HapMap data and FHS and drop any outliers as a final check:

system("plink --bfile hapmapCEU_f --freq --out hapMAF")

system("plink --bfile 11.step.update --freq --out fhsMAF")

# read if AF data:
# 
fhsMAF <- read.table("fhsMAF.frq",header=TRUE)
hapMAF <- read.table("hapMAF.frq",header=TRUE)

mergedMAF <- merge(hapMAF,fhsMAF,"SNP")

MAFdiff <- (mergedMAF$MAF.x)-(mergedMAF$MAF.y) # value of MAF differences (nice normal distribution)

MAFwDIFF <- cbind(mergedMAF,MAFdiff)# merge data

drop.index <- c(which(MAFwDIFF$MAFdiff>(2*sd(MAFwDIFF$MAFdiff))),which(MAFwDIFF$MAFdiff<(-2*sd(MAFwDIFF$MAFdiff)))) # create drop index of all SNPs greater than 2 standard deviations from mean of difference distribution

droplist2 <- as.character(MAFwDIFF[drop.index,1])

write.table(droplist2,file="dropissues2.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

# Now drop outlier SNPs (very conservative definition of outliers)
#  NOTE THIS STEP IS DANGEROUS, FOR THE SAKE OF BREVITY, FILES ARE BEING OVERWRITTEN

system("plink --bfile hapmap_f --exclude dropissues2.txt --make-bed --out hapmap_f")

system("plink --bfile 11.step.update --exclude dropissues2.txt --make-bed --out 11.step.update")

# Now merge the data:

system("plink --bfile hapmap_f --bmerge 11.step.update.bed 11.step.update.bim 11.step.update.fam --make-bed --out FHS.HM.merge")

# Get list of SNPs in linkage equilibrium and create independent SNP data set

system("nohup plink --bfile FHS.HM.merge --indep-pairwise 200 5 0.25 &") # this LE value is based on ISC LD pruning value
#
system("plink --bfile FHS.HM.merge --extract plink.prune.in --make-bed --out 11.step.LE &")

# Only SNPs in LE were used to generate the 20 principle dimensions to be used as covariates
# Generate IBS matrix of genetic distance to be used for covariates:
system("nohup plink --bfile 11.step.LE --cluster --mds-plot 20 --out 12.step.MDS &") ##

# 
# read in MDS data:

MDS <- read.table("12.step.MDS.mds",header=TRUE)

# Plot MDS data:

sample.index <- 181:nrow(MDS)
JPT.index <- 1:60
YRI.index <- 121:180
CEU.index <- 61:120

# Plot MDS dimensions 1 & 2 to show genetic relatedness of sample subjects relative to Hapmap
plot(MDS$C1, MDS$C2,main="MDS Genetic Relatedness of Individuals",xlab="Dimension1",ylab="Dimension2")
points(MDS$C1[sample.index],MDS$C2[sample.index],col="blue") # sample points are blue
points(MDS$C1[CEU.index],MDS$C2[CEU.index],col="green") # CEU Points are green
points(MDS$C1[YRI.index],MDS$C2[YRI.index],col="orange") # YRI points are orange
points(MDS$C1[JPT.index],MDS$C2[JPT.index],col="purple") # JPT/CHB points are red
legend(-.15,.08,c('Sample','HapMap CEU','HapMap YRI','HapMap JPT+CHB'),text.col=c('blue','green','orange','purple'))

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
distance <- cbind((as.vector(rdist(center.cau,IND.coord))),MDS[,1:2]) # create a vector of the distance of each individual away from the centroid of the sample cluster merged with FID & IID
names(distance) <- c('distance','FID','IID')

# generate list of individuals who are inside the inner radius to keep:
keep.list <- distance[which(distance$distance<(.1*SJC.dist)),2:3]

write.table(keep.list,file="keepstrat.list.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

system("plink --bfile 11.step.LE --keep keepstrat.list.txt --make-bed --out 11.step.LE")

##################################################################################
# Step 12: Prune data for Stratification & make sure all SNPs are in LE
##################################################################################


fam.data <- read.table("11.step.LE.fam", header=FALSE)
num.subjects <- nrow(fam.data)# each row in file is independent subject
system(paste("nohup plink --bfile 11.step.LE --cluster --neighbour 1 ",(num.subjects-1)," &",sep='')) # compare relatedmness of individuals based on SNPs in linkage equalibrium (MAKE SURE YOU ARE COMPARING NUMEBER OF SUBJECTS - 1, OTHERWISE CRASH BURN OCCURS

 # Check the test file to make sure it worked, but the correct file's name is the one listed below to be read in...
IBS.data <- read.table("plink.nearest",header=TRUE) # IBS data with Z-scores and genetic distance

## Find individual who is most related to everyone else in sample by finding lowest average genetic distance ind

IID.meandist <- tapply(IBS.data$MIN_DST, INDEX=IBS.data$IID,FUN=mean,na.rm=TRUE) # this gives average genetic difference for all individuals, the one with the lowest score is centroid

center.ID <- names(IID.meandist[(which(IID.meandist==min(IID.meandist)))]) # This is the subject ID of the person who is most simililar in terms of IBS to everyone else in the sample

center.index <- which(IBS.data$IID==center.ID) # Indeces in IBS data of center ID rows

center.IBS.data <- IBS.data[center.index,]

hist(center.IBS.data$MIN_DST,10000, xlab="Genetic Distance From Centroid", main="Distribution of IBS relatedness") # histogram of relatedness in sample

## Determine 3 standard deviation cut off distance value and remove bottom percentile from sample:

sort.index <- order(center.IBS.data$MIN_DST)
sorted.IBS.data <- center.IBS.data[sort.index,]
final.IBS.data <- sorted.IBS.data[1:(nrow(sorted.IBS.data)*.997),]
# Now Check with Histogram:info

hist(final.IBS.data$MIN_DST,10000, xlab="Genetic Distance From Centroid", main="Distribution of IBS relatedness Post Outlier Pruning") # histogram of relatedness in sample

# Keep List of individuals who pass stratification test:

keep.strat.list <- rbind(cbind(as.character(final.IBS.data[1,1]),as.character(final.IBS.data[1,2])),cbind(as.character(final.IBS.data[,6]),as.character(final.IBS.data[,7]))) # list of FID and IID for individuals to keep

write.table(keep.strat.list,file="keep.strat.list.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

system("plink --bfile 11.step.LE --keep keep.strat.list.txt --make-bed --out 12.step.LEstrat")
system("plink --bfile 10.step.unrelated --keep keep.strat.list.txt --make-bed --out 12.step.LDstrat") # remove individuals from full LD data set that don't pass stratification

#####################################################################
########## ***** END OF ALL DATA CLEANING PROCEDUERS ***** ##########
#####################################################################



           
#####################################################################
########## *****STANDARD GWAS ANALYSIS OF ENTIRE DATA SET*****   ####
#####################################################################
# Generate alternate phenotype file:
# Format: 3 columns FID, IID, PHE

fam.data <- read.table("12.step.LDstrat.fam",header=FALSE) # Individuals with SNPs still in LD
names(fam.data)<- c("FID","IID")

# read in phenotype data:
risk.scores <- read.table("Risk_Scores.csv",header=FALSE,skip=1,sep=",")
# find indeces of individuals that have a risk score
notna.index <- as.logical((is.na(risk.scores$V6)*-1)+1)
risk.data<- risk.scores[notna.index,] # only individuals with risk scores are included in the data now
names(risk.data)<- c("IID")
merged.data <- merge(risk.data,fam.data,"IID") # merge fam data and risk data

pheno.file <- as.data.frame(cbind(as.character(merged.data[,7]),as.character(merged.data[,1]),as.character(merged.data[,6])))
names(pheno.file) <- c("FID","IID","PHE")
# write risk score data:

write.table(pheno.file,file="risk.score.pheno",quote=FALSE,row.names=FALSE,col.names=FALSE)

######### Generate and Write out file with log transformed risk data:

log.pheno.pre <- pheno.file
e <- (1+(1/10000000))^10000000 # define constant 'e'
log.pheno.data <- log(as.numeric(as.character(pheno.file[,3])),base=e)
hist(log.pheno.data,100) # look at histogram to confirm normal distribution
log.pheno.pre[,3] <- log.pheno.data
log.pheno.file <- log.pheno.pre

write.table(log.pheno.file,file="log.score.pheno",quote=FALSE,row.names=FALSE,col.names=FALSE)


# write list of individuals to keep:

keep <- cbind(as.character(pheno.file[,1]),as.character(pheno.file[,2]))
write.table(keep,file="GWASkeep",quote=FALSE,row.names=FALSE,col.names=FALSE)

# Generate covariate file: 
# Format: Principal components and batch effects per IID

batch <- read.table("batch.covar",header=FALSE) # read in batch data
names(batch) <- c('FID','IID','BatchD1','BatchD2') # give batch col's names

MDS.data <- read.table("12.step.MDS.mds",header=TRUE) # read in dimension data
sample.index <- 181:nrow(MDS.data)

4cv <- merge(MDS.data[sample.index,],batch[,2:4],"IID")

# ADD AGE AND BATCH AS COVARIATES
cv.pre.file <- as.data.frame(merge(cv,merged.data,"IID"))

cv.file <- cbind(cv.pre.file[,2],cv.pre.file[,1],cv.pre.file[,29],cv.pre.file[,4:25])

namelist <- names(cv.file) # get col names
namelist[1] <- c('FID')
namelist[2] <- c('IID')
namelist[3] <- c('AGE')
names(cv.file) <- namelist # Fix variable names

# write file:

write.table(cv.file,file="cv.covar",quote=FALSE,row.names=FALSE,col.names=TRUE)

###################################################################################################################
# Step 13: Run GWAS on Full Data Set controlling for 20 most significant Principal Components and Batch effects 
###################################################################################################################
# NOTE: USE ALL SNPS IN LD FOR ACTUAL ASSOCIATION
                     
#NEXT STEP: Run GWAS while controlling for PComponents as covariates

# Run basic Quantitative association with SNPs in LD:
                     
system("nohup plink --bfile 12.step.LDstrat --make-founders --linear --sex --covar cv.covar --pheno risk.score.pheno --out 13.GWAS &") # non-transformed

# Now Transformed:

system("nohup plink --bfile 12.step.LDstrat --make-founders --linear --sex --covar cv.covar --pheno log.score.pheno --out 13.log.GWAS &")


# Read in GWAS results:
#gwas.data <-  read.table("13.GWAS.assoc.linear",header=TRUE)

# If transformed phenotype:
gwas.data <-  read.table("13.log.GWAS.assoc.linear",header=TRUE)

add.index <- which(gwas.data$TEST=='ADD') # only include those rows that are SNPs from regression
gwas.data <- gwas.data[add.index,]
#
o.index <- order(gwas.data$P)
gwas.data <- gwas.data[o.index,] # sorted with lowest P-vals at top

# Now remove SNPs that were due to bad calls that in the top 30 lowest p-value range
keep.index <- c(2:nrow(gwas.data))
gwas.data <- gwas.data[keep.index,]

# Now order things correctly based on Chromosome and Base Position:

order.index <- order(gwas.data$CHR, gwas.data$BP) # get index to sort by chromosome, then baseposition within chromosome
gwas.data <- gwas.data[order.index,]
#

# write out list of SNPs to keep after removing bad calls from top lowest 30 p-values
good.snps <- as.data.frame(gwas.data$SNP) 
write.table(good.snps,file="good.snps",quote=FALSE,row.names=FALSE,col.names=TRUE)

# FIX NEXT 4 LINES:
data <- as.data.frame(cbind(as.numeric(as.character(gwas.data$CHR)),as.numeric(as.character(gwas.data$BP)),as.numeric(as.character(gwas.data[,9]))))# unadjusted p
names(data) <- c('CHR','BP','P')

### Look at Plot of GWAS P values and QQ plot:
# Note: exclude X,Y, 0 (unplaced), and XY (psuedo-autosomal)

# GWAS PLOTS:
# NOTE: BE SURE TO GIVE SITATION

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ALSO: REMOVE BAD CALL SNPS BEFORE PLOTTING

library(ggplot2)
library(gap)

# NOW PLOT THE DATA FROM GWAS NOW THAT FUNCTIONS ARE INSTALLED

rC0.index <- which(data$CHR!=0)# remove CHR 0
data <- data[rC0.index,]

rC23.index <- which(data$CHR!=23)# remove CHR 23
data <- data[rC23.index,]

rC24.index <- which(data$CHR!=24)# remove CHR 24
data <- data[rC24.index,]

rC25.index <- which(data$CHR!=25)# remove CHR 25
data <- data[rC25.index,]

# GENERATE PLOT:

color <- rep(c("blue", "red"), 11) 
par(las = 2, xpd = TRUE, cex.axis = 1.8, cex = 0.4)
ops <- mht.control(colors = color, yline = 1.5, xline = 3)
mhtplot(data, ops, pch = 25) # manhattan plot
axis(2, pos = 2, at = 1:16)
#abline(h=8,col=c('yellow'))
title("Manhattan Plot Of Heart Disease Risk P-Values", cex.main = 4)
#
qqunif(data$P) # qq-plot
title("QQ-Plot Of Heart Disease Risk P-Values", cex.main = 2)

##########################################################
########## ***** END OF FULL DATASET GWAS ***** ##########
##########################################################




#####################################################################
########## ***** PREDICTIVE POLYGENIC ANALYSIS *****       ##########
#####################################################################

########################################################################
# Step 14: split sample into training set and test set split, be sure each training and test set are balanced for distribution of phenotypic risk
########################################################################

# read in individual list:
# read in phenotype files:
       
phenotype <- read.table("risk.score.pheno",header=FALSE) 
log.phenotype <- read.table("log.score.pheno",header=FALSE)

#############Plots of data:
num.ids <- nrow(phenotype) # total number of subjects
n.crossvals <- 10 # number of crossval groups

plot(log.phenotype$V3, main='Log Transformed Phenotype Data',xlab='Test Sets = Non-Redundant Subjects',ylab='Rank')
abt.val <-seq(1,num.ids,by=(num.ids/n.crossvals))
abline(v=abt.val,col='red')
#

plot(phenotype$V3, main='Untransformed Phenotype Data',xlab='Test Sets =  Non-Redundant Subjects',ylab='Risk Score')
abt.val <-seq(1,num.ids,by=(num.ids/n.crossvals))
abline(v=abt.val,col='red')
#
###############
# generate the 10 test sets, each composed of a different 1/10 of the sample
n.crossvals <- 10 # number of crossval groups

for(i in 1:n.crossvals){
  eval(parse(text=paste('test.data',i,' <- phenotype[(((num.ids/(n.crossvals)*',i,')-(num.ids/n.crossvals)+1):((num.ids/(n.crossvals)*',i,'))),]',sep="")))
}


# now write out the Test subject files
for (i in 1:n.crossvals){
eval(parse(text=paste('write.table(test.data',i,'[,1:2],file="TESTkeep',i,'",quote=FALSE,row.names=FALSE,col.names=FALSE)',sep="")))
}


########################################################################
# Step 15: Run second GWAS on training set series using cross validation
########################################################################
# first remove listed TESTkeep's and generate training set

# First for untransformed
for (i in (1:10)){
eval(parse(text=paste('system("nohup plink --bfile 12.step.LDstrat --make-founders --extract good.snps --linear --sex --covar cv.covar --pheno risk.score.pheno --exclude TESTkeep',i,' --out 15.GWAS.',i,' &")',sep="")))
}

# Second for transformed

for (i in (1:10)){
eval(parse(text=paste('system("nohup plink --bfile 12.step.LDstrat --make-founders --extract good.snps --linear --sex --covar cv.covar --pheno log.score.pheno --exclude TESTkeep',i,' --out 15.log.GWAS.',i,' &")',sep="")))
}
# 
#
#################################################################################################
# Step 16: Generate Score data from Training Set: 
#################################################################################################
# first generate .raw file from beta's.
# read in GWAS results from all 10 training sets:
n.crossvals <- 10
for (i in 1:n.crossvals){
eval(parse(text=paste('gwas.data.',i,' <-  read.table("15.log.GWAS.',i,'.assoc.linear",header=TRUE)',sep='')))
print(paste('data set ',i,' read in',sep=''))
eval(parse(text=paste('add.index <- which(gwas.data.',i,'$TEST=="ADD")',sep=''))) # only include those rows that are SNPs from regression
eval(parse(text=paste('gwas.data.',i,' <- gwas.data.',i,'[add.index,]',sep='')))
}

#
# make sure all numeric data is in correct format
# create updated bed file for next steps:
system("plink --bfile 12.step.LDstrat --make-bed --out 16.step") # create updated bed file for this step

# read in score and phenotype files:
       
phenotype <- read.table("log.score.pheno",header=FALSE) # no FID col
phe.o.index <- order(phenotype[,3])
phenotype <- phenotype[phe.o.index,]
phenotype <- cbind(phenotype,c(1:nrow(phenotype)))
names(phenotype) <- c('FID','IID','PHE')
num.ids <- length(1:nrow(phenotype))


# This loop generates all remaining file for cross validation:                                        
for (i in 1:n.crossvals){
  
eval(parse(text=paste('gwas.data <- gwas.data.',i,'',sep='')))# assign gwas data from current iteration of cross validation
print(paste('gwas data ',i,' read in',sep=''))
# sort raw data by P-value and create seperate files:
o.index <- order(gwas.data$P)
gwas.data <- gwas.data[o.index,] # sorted with lowest P-vals at top
raw <- gwas.data

# generate score file:
snpscore.dat <- as.data.frame(cbind(as.character(raw$SNP),as.character(raw$A1),raw$P))
write.table(snpscore.dat,file='snpprofile.raw',quote=FALSE,row.names=FALSE,col.names=FALSE)

# generate SNPval file:
snpval.dat <- as.data.frame(cbind(as.character(raw$SNP),raw$P))
write.table(snpval.dat,file='snpval.dat',quote=FALSE,row.names=FALSE,col.names=FALSE)

# generate score ranges file:

range.dat <- as.data.frame(matrix(c('S1','S2','S3','S4','S5','S6','S7','S8','S9','S10','S11',0.00,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.00,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1),nrow=11,byrow=FALSE))
write.table(range.dat,file='q.ranges',quote=FALSE,row.names=FALSE,col.names=FALSE)

# All SNPs:

# First for K 1:10 of cross validation
print(paste('performing cross validation ',i,'',sep=''))
eval(parse(text=paste('system("plink --bfile 16.step --keep TESTkeep',i,' --score snpprofile.raw --q-score-file snpval.dat --q-score-range q.ranges --out 16.step.V',i,'")',sep='')))
print(paste('cross validation ',i,' complete!',sep=''))
}

# 
# Combine Data For Models:
for (i in 1:10){
for (j in 1:11){       
eval(parse(text=paste('profile.1 <- read.table("16.step.V',i,'.S',j,'.profile",header=TRUE)',sep='')))
eval(parse(text=paste('top.datV',i,'S',j,' <- merge(profile.1,phenotype,"IID")',sep='')))   
}
} # Finish loop and examine cross-validation models

# merge data together for different score ranges:

total.range1 <- rbind(top.datV1S1,top.datV2S1,top.datV3S1,top.datV4S1,top.datV5S1,top.datV6S1,top.datV7S1,top.datV8S1,top.datV9S1,top.datV10S1)
range1 <- summary(lm(PHE~SCORE,data=total.range1))

total.range2 <- rbind(top.datV1S2,top.datV2S2,top.datV3S2,top.datV4S2,top.datV5S2,top.datV6S2,top.datV7S2,top.datV8S2,top.datV9S2,top.datV10S2)
range2 <- summary(lm(PHE~SCORE,data=total.range2))

total.range3 <- rbind(top.datV1S3,top.datV2S3,top.datV3S3,top.datV4S3,top.datV5S3,top.datV6S3,top.datV7S3,top.datV8S3,top.datV9S3,top.datV10S3)
range3 <- summary(lm(PHE~SCORE,data=total.range3))

total.range4 <- rbind(top.datV1S4,top.datV2S4,top.datV3S4,top.datV4S4,top.datV5S4,top.datV6S4,top.datV7S4,top.datV8S4,top.datV9S4,top.datV10S4)
range4 <- summary(lm(PHE~SCORE,data=total.range4))

total.range5 <- rbind(top.datV1S5,top.datV2S5,top.datV3S5,top.datV4S5,top.datV5S5,top.datV6S5,top.datV7S5,top.datV8S5,top.datV9S5,top.datV10S5)
range5 <- summary(lm(PHE~SCORE,data=total.range5))

total.range6 <- rbind(top.datV1S6,top.datV2S6,top.datV3S6,top.datV4S6,top.datV5S6,top.datV6S6,top.datV7S6,top.datV8S6,top.datV9S6,top.datV10S6)
range6 <- summary(lm(PHE~SCORE,data=total.range6))

total.range7 <- rbind(top.datV1S7,top.datV2S7,top.datV3S7,top.datV4S7,top.datV5S7,top.datV6S7,top.datV7S7,top.datV8S7,top.datV9S7,top.datV10S7)
range7 <- summary(lm(PHE~SCORE,data=total.range7))

total.range8 <- rbind(top.datV1S8,top.datV2S8,top.datV3S8,top.datV4S8,top.datV5S8,top.datV6S8,top.datV7S8,top.datV8S8,top.datV9S8,top.datV10S8)
range8 <- summary(lm(PHE~SCORE,data=total.range8))

total.range9 <- rbind(top.datV1S9,top.datV2S9,top.datV3S9,top.datV4S9,top.datV5S9,top.datV6S9,top.datV7S9,top.datV8S9,top.datV9S9,top.datV10S9)
range9 <- summary(lm(PHE~SCORE,data=total.range9))

total.range10 <- rbind(top.datV1S10,top.datV2S10,top.datV3S10,top.datV4S10,top.datV5S10,top.datV6S10,top.datV7S10,top.datV8S10,top.datV9S10,top.datV10S10)
range10 <- summary(lm(PHE~SCORE,data=total.range10))


############################################################################
# Generate Boxplot:
r.sq <- c(0,0,0.004177,0,0.002056,0.001383,0,0.0002382,0,0)
sets <- c('>0-.1','>.1-.2','>.2-.3','>.3-.4','>.4-.5','>.5-.6','>.6-.7','>.7-.8','>.8-.9','>.9-1')
barplot(r.sq, ylab='Adjusted R-Square',xlab='P-value Discovery Set',names.arg=sets, main='Framingham Heart Disease',col=c("darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","lightgreen","lightgreen","lightgreen"))
colorz<- c('darkgreen','lightgreen')
legend("topright", inset=.05,c('Significant','Null'),fill=colorz)
#

##############################################################################
########## ***** HERITABILITY POLYGENIC CLEANING AND ANALYSIS ***** ##########
##############################################################################                     

# DUE TO CLEANING PROCEDURES BEING ESSENTIALLY THE SAME, YET MORE STRINGENT, SAME FILES FROM ABOVE WERE USED THEN MORE STRINGENT THRESHOLDS WERE APPLIED WHERE NEEDED

##################################################################################
# Step 17:  Prune SNPs out of HWE
##################################################################################
system("plink --bfile 4.step.snpcall --hardy")

HWE.data <- read.table("plink.hwe",header=TRUE)
out.HWE <- which(HWE.data[,9]<.05)# <.05
prune.HWE.snps <- HWE.data[out.HWE,2]

write.table(prune.HWE.snps,file="17.step.remove.HWE.list.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

system("plink --bfile 4.step.snpcall --exclude 17.step.remove.HWE.list.txt --make-bed --out 17.step.HWE")

##################################################################################
# Step 18: non-random genotyping failure, as inferred by the flanking haplotypic background 
##################################################################################

# drop SNPs with P < 10^-10
system("plink --bfile 17.step.HWE --test-mishap --out 18.step") # generate mishap file

mishap <- read.table("18.step.missing.hap",header=FALSE,fill=TRUE,skip=TRUE)

fail.index <- which(as.numeric(as.character(mishap$V8))<(10^-10))

drop.hap<- as.data.frame(mishap$V1[fail.index])

write.table(drop.hap,file="18.step.remove.hap.list.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

system("plink --bfile 17.step.HWE --exclude 18.step.remove.hap.list.txt --make-bed --out 19.step.hap")

##########################################
#****************************************#
##########################################
# QC for individuals:


#################################################################################
# Step 19: Remove individuals who are missing more than 1% of SNP calls
#################################################################################
# First step is generating missingness data with PLINK

system("plink --bfile 19.step.hap --missing --out 19.step")

ind.miss <- read.table("19.step.imiss",header=TRUE)
above.one <- which(ind.miss[,6]>.01)
high.missing <- ind.miss[above.one,]

# Create a text file with IID and FID for individuals with missingness above 1%

miss.dat <- cbind(as.character(high.missing$FID),as.character(high.missing$IID))

write.table(miss.dat,file="19.step.remove.miss.list.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

system("plink --bfile 6.step.hap --remove 19.step.remove.miss.list.txt --make-bed --out 20.step.miss")

#################################################################################################################################################################
# Step 20:  evidence of gross non-random plate failure: Use Batch effects inferred from missingness to generate missingness covariate to be used in association
################################################################################################################################################################
# Generate phenotype data for plink file:
# ALREADY DONE ABOVE

#### Now for covariate:
ind.miss <- read.table("19.step.imiss",header=TRUE) # read in individual missingness data

merged.covar <- merge(merged.data,ind.miss,"IID") # merge fam data and covar

covar.file <- as.data.frame(cbind(as.character(merged.covar[,7]),as.character(merged.covar[,1]),as.character(merged.covar[,14])))
names(covar.file) <- c("FID","IID","N_MISS")
# write out covariate files:

write.table(covar.file,file="batch2.covar",quote=FALSE,row.names=FALSE,col.names=FALSE)

##################################################################################
# Step 21: Prune data outliers with respect to estimated heterozygosity
##################################################################################

system("plink --bfile 20.step.miss --het --out 21.step.het")
hetchk <- read.table("21.step.het.het",header=TRUE)

hist(hetchk$F[hetchk$F<.5], breaks=100) # plot distribution of heterozygosity

Three.sd <- 3*(sd(hetchk$F)) # 3 sd's 

drop.het.index <- c(which(hetchk$F>Three.sd), which(hetchk$F<(-1*Three.sd))) # create index of individuals outside 3sd of mean het
# read in het file and prune individuals outside of 3 standard deviations from the mean

drop.het.list <- hetchk[drop.het.index,1:2]

write.table(drop.het.list,file="remove.het.list.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

system("plink --bfile 20.step.miss --remove remove.het.list.txt --make-bed --out 21.step.drop.het")

#################################################################################
# Step 22: Remove individuals with discrepencies between reported and genotypic sex
#################################################################################
# SAME INDIVIDUALS AS ABOVE
       
system("plink --bfile 21.step.drop.het --remove remove.sex.list.txt --make-bed --out 22.step.drop.sex")

##################################################################################
# Step 23: Prune data for very closely related indviduals
##################################################################################
# First prune SNPs in LD                      

# Get list of SNPs in linkage equilibrium and create independent SNP data set

system("plink --bfile 22.step.drop.sex --out L.equalib --indep-pairwise 200 5 0.25") # this LE value is based on ISC LD pruning value

system("plink --bfile 22.step.drop.sex --extract L.equalib.prune.in --make-bed --out 23.step.LE")
                      
#Scan the plink.genome file for any individuals with high PIHAT values (.025 for Polygenic Heritability Analysis)

system("nohup plink --bfile 23.step.LE --genome --min 0.05 --out 23.close.relatedLE &") # Create list of closely related individuals
                      
relchk <- read.table("23.close.related.genome",header=TRUE)

d <- which(relchk$PI_HAT>.025)

FID.IID <- as.data.frame(cbind(as.character(relchk$FID1[d]),as.character(relchk$IID1[d])))

write.table(FID.IID,file="23remove.rel.list.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

# use system command to remove discrepant individuals from data set and write new data files

system("plink --bfile 22.step.drop.sex --make-founders --remove 23remove.rel.list.txt --make-bed --out 23.step.unrelated")# use file generated from this line after list of individuals post stratification cleaning is generated

                      
system("plink --bfile 23.step.LE --make-founders --remove 23remove.rel.list.txt --make-bed --out 23.step.unrelatedLE")# this line is for speedy stratification
                      
##################################################################################
# Step 24.a Multidimensional scaling analysis with HapMap (First part of Stratification cleaning)
##################################################################################
# NEXT STEP:
# merge Framingham with HapMap, THEN prune for LE: 
## convert ss# to rs# so SNPs data can be merged:

#                     
SNPlist <- read.table("23.step.unrelated.bim",header=FALSE)# read in list of all SNPs from FHS data

# Read in file with rs# (marker-info file)SNPs:

SNPmi <- read.table("FHSmarker-info.csv",header=FALSE,sep=",")

markers <- SNPmi[,4:5]

# write out SNP names file:

write.table(markers,file="update24SNP.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

# Use --update-name option in plink under data management to rename SNPs so they have rs#'s

system("plink --bfile 23.step.unrelated --update-map update24SNP.txt --update-name --make-bed --out 24.step.update")

# check the new SNP names:

SNPlist <- read.table("24.step.update.bim",header=FALSE)

##################################################################################
# Now merge FHS and Hapmap founders:

### Merge hapmap founders with FHS data and resolve strand flip and positional issues:

# first get list of detected strand and positional issues from plink

system("plink --bfile hapmap_f --bmerge 24.step.update.bed 24.step.update.bim 24.step.update.fam --make-bed --out 24.FHS.HM.merge")

# read in detected strand flip issues:

strand <- read.table("24.FHS.HM.merge.missnp",header=FALSE)
posit <- read.table("DiffPos.txt",header=FALSE,fill=TRUE)

droplist <- rbind(strand,posit)

# write list of SNPs to drop

write.table(droplist,file="dropissues24.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

# remove bad SNPs from FHS data:

system("plink --bfile 24.step.update --exclude dropissues24.txt --make-bed --out 24.step.update")

# remove bad SNPs from HapMap:

system("plink --bfile hapmap_f --exclude dropissues.txt --make-bed --out hapmap_f")

# get MAF from each CEU HapMap data and FHS and drop any outliers as a final check:

system("plink --bfile hapmapCEU_f --freq --out hapMAF")

system("plink --bfile 24.step.update --freq --out 24fhsMAF")

# read if AF data:
# 
fhsMAF <- read.table("24fhsMAF.frq",header=TRUE)
hapMAF <- read.table("hapMAF.frq",header=TRUE)

mergedMAF <- merge(hapMAF,fhsMAF,"SNP")

MAFdiff <- (mergedMAF$MAF.x)-(mergedMAF$MAF.y) # value of MAF differences (nice normal distribution)

MAFwDIFF <- cbind(mergedMAF,MAFdiff)# merge data

drop.index <- c(which(MAFwDIFF$MAFdiff>(2*sd(MAFwDIFF$MAFdiff))),which(MAFwDIFF$MAFdiff<(-2*sd(MAFwDIFF$MAFdiff)))) # create drop index of all SNPs greater than 2 standard deviations from mean of difference distribution

droplist2 <- as.character(MAFwDIFF[drop.index,1])

write.table(droplist2,file="24dropissues.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

# Now drop outlier SNPs (very conservative definition of outliers)
#  NOTE THIS STEP IS DANGEROUS, FOR THE SAKE OF BREVITY, FILES ARE BEING OVERWRITTEN

system("plink --bfile hapmap_f --exclude 24dropissues.txt --make-bed --out hapmap_f")

system("plink --bfile 24.step.update --exclude 24dropissues.txt --make-bed --out 24.step.update")

# Now merge the data:

system("plink --bfile hapmap_f --bmerge 24.step.update.bed 24.step.update.bim 24.step.update.fam --make-bed --out 24.FHS.HM.merge")

# Get list of SNPs in linkage equilibrium and create independent SNP data set

system("nohup plink --bfile 24.FHS.HM.merge --indep-pairwise 200 5 0.25 --out 24 &") # this LE value is based on ISC LD pruning value
                      
# 
system("plink --bfile 24.FHS.HM.merge --extract 24.prune.in --make-bed --out 24.step.LE")

# Only SNPs in LE were used to generate the 20 principle dimensions to be used as covariates
# Generate IBS matrix of genetic distance to be used for covariates:
system("nohup plink --bfile 24.step.LE --cluster --mds-plot 20 --out 24.step.MDS &") ##

#
# read in MDS data:

MDS <- read.table("24.step.MDS.mds",header=TRUE)

# Plot MDS data:

sample.index <- 181:nrow(MDS)
JPT.index <- 1:60
YRI.index <- 121:180
CEU.index <- 61:120

# Plot MDS dimensions 1 & 2 to show genetic relatedness of sample subjects relative to Hapmap
plot(MDS$C1, MDS$C2,main="MDS Genetic Relatedness of Individuals",xlab="Dimension1",ylab="Dimension2")
points(MDS$C1[sample.index],MDS$C2[sample.index],col="blue") # sample points are blue
points(MDS$C1[CEU.index],MDS$C2[CEU.index],col="green") # CEU Points are green
points(MDS$C1[YRI.index],MDS$C2[YRI.index],col="orange") # YRI points are orange
points(MDS$C1[JPT.index],MDS$C2[JPT.index],col="purple") # JPT/CHB points are red
legend(-.15,.08,c('Sample','HapMap CEU','HapMap YRI','HapMap JPT+CHB'),text.col=c('blue','green','orange','purple'))

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
distance <- cbind((as.vector(rdist(center.cau,IND.coord))),MDS[,1:2]) # create a vector of the distance of each individual away from the centroid of the sample cluster merged with FID & IID
names(distance) <- c('distance','FID','IID')

# generate list of individuals who are inside the inner radius to keep:
keep.list <- distance[which(distance$distance<(.1*SJC.dist)),2:3]

write.table(keep.list,file="24keepstrat.list.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

system("plink --bfile 24.step.LE --keep 24keepstrat.list.txt --make-bed --out 24.step.LE")

                      
##################################################################################
# Step 24: Prune data for Stratification (Second Part)
##################################################################################
#                    
fam.data <- read.table("24.step.LE.fam", header=FALSE)
num.subjects <- nrow(fam.data)# each row in file is independent subject                       

system(paste("plink --bfile 24.step.LE --cluster --neighbour 1 ",(num.subjects-1)," --out 24",sep='')) # compare relatedmness of individuals based on SNPs in linkage equalibrium (MAKE SURE YOU ARE COMPARING NUMEBER OF SUBJECTS - 1, OTHERWISE CRASH BURN OCCURS
                  
 # Check the test file to make sure it worked, but the correct file's name is the one listed below to be read in...
IBS.data <- read.table("24.nearest",header=TRUE) # IBS data with Z-scores and genetic distance

## Find individual who is most related to everyone else in sample by finding lowest average genetic distance ind

IID.meandist <- tapply(IBS.data$MIN_DST, INDEX=IBS.data$IID,FUN=mean,na.rm=TRUE) # this gives average genetic difference for all individuals, the one with the lowest score is centroid

center.ID <- names(IID.meandist[(which(IID.meandist==min(IID.meandist)))]) # This is the subject ID of the person who is most simililar in terms of IBS to everyone else in the sample

center.index <- which(IBS.data$IID==center.ID) # Indeces in IBS data of center ID rows

center.IBS.data <- IBS.data[center.index,]

hist(center.IBS.data$MIN_DST,1000, xlab="Genetic Distance From Centroid", main="Distribution of IBS relatedness") # histogram of relatedness in sample

## Determine 3 standard deviation cut off distance value and remove bottom percentile from sample:

sort.index <- order(center.IBS.data$MIN_DST)
sorted.IBS.data <- center.IBS.data[sort.index,]
final.IBS.data <- sorted.IBS.data[1:(nrow(sorted.IBS.data)*.997),]
# Now Check with Histogram:info

hist(final.IBS.data$MIN_DST,1000, xlab="Genetic Distance From Centroid", main="Distribution of IBS relatedness Post Outlier Pruning") # histogram of relatedness in sample

# Keep List of individuals who pass stratification test:

keep.strat.list <- rbind(cbind(as.character(final.IBS.data[1,1]),as.character(final.IBS.data[1,2])),cbind(as.character(final.IBS.data[,6]),as.character(final.IBS.data[,7]))) # list of FID and IID for individuals to keep

write.table(keep.strat.list,file="24keep.strat.list.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

# Use input file with SNPs in LD to generate file that contains correct individuals and correct SNPs
system("plink --bfile 23.step.unrelated --keep 24keep.strat.list.txt --make-bed --out 24.step.strat")


#####################################################################
########## ***** END OF ALL DATA CLEANING PROCEDUERS ***** ##########
#####################################################################       

       
#####################################################################
########## ***** HERITABILITY POLYGENIC ANALYSIS *****     ##########
#####################################################################

# Use GCTA (installed in current working directory)
                  
##########################################################################################
# Step 25: Generate final file with correct IID's and SNPs
##########################################################################################

# generate phenotype file for applicable subjects:

fam.data <- read.table("24.step.strat.fam",header=FALSE)# contains FID, IID, Sex for current subjects
names(fam.data) <- c('FID','IID')                      
# read in phenotype data:
risk.scores <- read.table("Risk_Scores.csv",header=FALSE,skip=1,sep=",")
# find indeces of individuals that have a risk score
notna.index <- as.logical((is.na(risk.scores$V6)*-1)+1)
risk.data<- risk.scores[notna.index,] # only individuals with risk scores are included in the data now
names(risk.data)<- c("IID")
merged.data <- merge(risk.data,fam.data,"IID") # merge fam data and risk data

pheno.file <- as.data.frame(cbind(as.character(merged.data[,7]),as.character(merged.data[,1]),as.character(merged.data[,6])))
names(pheno.file) <- c("FID","IID","PHE")
                      
# write list of individuals to keep:

keep <- cbind(as.character(pheno.file[,1]),as.character(pheno.file[,2]))
write.table(keep,file="26.GWASkeep",quote=FALSE,row.names=FALSE,col.names=FALSE)                      

system("plink --bfile 24.step.strat --keep 26.GWASkeep --make-bed --out 26.step.Var")
                      
# write risk score data:

write.table(pheno.file,file="26risk.score.pheno",quote=FALSE,row.names=FALSE,col.names=FALSE)

# write sex covariate file:  (must be seperate covariate file because sex is a discrete variable)                    
sex.file <- as.data.frame(cbind(as.character(merged.data[,7]),as.character(merged.data[,1]),as.character(merged.data[,4])))

names(sex.file) <- c("FID","IID","SEX")
                      
write.table(sex.file,file="26sex.covar",quote=FALSE,row.names=FALSE,col.names=FALSE)

# Generate continuous covariate file: 
# Format: Principal components and batch effects per IID


batch <- read.table("batch.covar",header=FALSE) # read in batch data
names(batch) <- c('FID','IID','BatchD1','BatchD2') # give batch col's names

MDS.data <- read.table("24.step.MDS.mds",header=TRUE) # read in dimension data
sample.index <- 181:nrow(MDS.data)

cv.pre <- merge(MDS.data[sample.index,],batch[,2:4],"IID")
cv.pre2 <- merge(cv.pre,sex.file,"IID")# ensure correct individuals are listed
cv <- cv.pre2[,1:(ncol(cv.pre2)-2)]                      

# ADD AGE AND BATCH AS COVARIATES
cv.pre.file <- as.data.frame(merge(cv,merged.data,"IID"))

cv.file <- cbind(cv.pre.file[,2],cv.pre.file[,1],cv.pre.file[,29],cv.pre.file[,4:25])

namelist <- names(cv.file) # get col names
namelist[1] <- c('FID')
namelist[2] <- c('IID')
namelist[3] <- c('AGE')
names(cv.file) <- namelist # Fix variable names

# write file:

write.table(cv.file,file="26cv.covar",quote=FALSE,row.names=FALSE,col.names=TRUE)

##########################################################################################
# Step 26: Estimation of the genetic relationship matrix (GRM) from all the autosomal SNPs
##########################################################################################

system("./gcta64 --bfile 26.step.Var --autosome --make-grm --out 26.step.GRMat")
                      
##########################################################################################
# Step 27: Estimation of the variance explained by the SNPs Using all 3 REML (restricted maximum likelihood) methods
##########################################################################################
                      
system("nohup ./gcta64 --grm 26.step.GRMat --reml-maxit 10000 --reml-alg 2 --reml --pheno 26risk.score.pheno --qcovar 26cv.covar --covar 26sex.covar --out 27.step.variance.EM &")

system("nohup ./gcta64 --grm 26.step.GRMat --reml-maxit 10000 --reml-alg 1 --reml --pheno 26risk.score.pheno --qcovar 26cv.covar --covar 26sex.covar --out 27.step.variance.Fisher &")

system("nohup ./gcta64 --grm 26.step.GRMat --reml-maxit 10000 --reml-alg 0 --reml --pheno 26risk.score.pheno --qcovar 26cv.covar --covar 26sex.covar --out 27.step.variance.AI &")
#                      

#****************************************************************************************************************#
###################################################################################################################
# Step 28:Estimation of the variance explained by the SNPs Using Haesman-Elston dropping fewer individuals using R:
###################################################################################################################

#First perform new cleaning procedures:


##################################################################################
# Step 29.a Multidimensional scaling analysis with HapMap (First part of Stratification cleaning)
##################################################################################
# Generate Covariates from MDS for regression model

#                     
SNPlist <- read.table("22.step.drop.sex.bim",header=FALSE)# read in list of all SNPs from FHS data

# Read in file with rs# (marker-info file)SNPs:

SNPmi <- read.table("FHSmarker-info.csv",header=FALSE,sep=",")

markers <- SNPmi[,4:5]

# write out SNP names file:

write.table(markers,file="update24SNP.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

# Use --update-name option in plink under data management to rename SNPs so they have rs#'s

system("plink --bfile 22.step.drop.sex --update-map update24SNP.txt --update-name --make-bed --out 29.step.update")

# check the new SNP names:

SNPlist <- read.table("29.step.update.bim",header=FALSE)

##################################################################################
# Now merge FHS and Hapmap founders:

### Merge hapmap founders with FHS data and resolve strand flip and positional issues:

# remove bad SNPs from FHS data:

system("plink --bfile 29.step.update --exclude dropissues24.txt --make-bed --out 29.step.update")

system("plink --bfile 29.step.update --exclude 24dropissues.txt --make-bed --out 29.step.update")

# Now merge the data:

system("plink --bfile hapmap_f --bmerge 29.step.update.bed 29.step.update.bim 29.step.update.fam --make-bed --out 29.FHS.HM.merge")

# This is genetic relatedness of sample individuals and hapmap:
system("nohup gcta --bfile 29.FHS.HM.merge --make-grm --autosome --out 29.step.MDS &")

################## Correct for Stratification: 
# read in data for MDS:
w.hap <- read.table('29.step.MDS.grm',header=FALSE)

grm.matrix.MDS <- matrix(NA,nrow=length(unique(w.hap$V1)),ncol=length(unique(w.hap$V1)))

# First Bottom Half of Matrix
for (i in 1:nrow(w.hap)){
grm.matrix.MDS[w.hap$V1[i],w.hap$V2[i]] <- w.hap$V4[i]   
}#

# Now Top Half of Matrix
for (i in 1:nrow(w.hap)){
grm.matrix.MDS[w.hap$V2[i],w.hap$V1[i]] <- w.hap$V4[i]   
}#

dist.MDS <- max(grm.matrix.MDS,na.rm=TRUE)-grm.matrix.MDS

# generate MDS data:

etest<- cmdscale(dist.MDS,k=20,eig=TRUE)

# place eigenvectors in matrix and write out

evect <- as.matrix(etest$points)#
# now merge with FID and IID
MDS.fam <- read.table("29.FHS.HM.merge.fam",header=FALSE)
evect2 <- cbind(MDS.fam[,c(1:2)],evect)

write.table(evect2,file="MDS.20.HE",quote=FALSE,row.names=FALSE,col.names=FALSE)

# examine with scree plot:
plot(etest$eig,type='b',main='Scree Plot of Eigenvectors',xlab='dimension',ylab='relative variance explained')

##################################################################################
# Step 24: Prune data for Stratification (Second Part)
##################################################################################
#                    
fam.data <- read.table("29.step.LE.fam", header=FALSE)
num.subjects <- nrow(fam.data)# each row in file is independent subject                       

system(paste("nohup plink --bfile 29.step.LE --cluster --neighbour 1 ",(num.subjects-1)," --out 29 &",sep='')) # compare relatedmness of individuals based on SNPs in linkage equalibrium (MAKE SURE YOU ARE COMPARING NUMEBER OF SUBJECTS - 1, OTHERWISE CRASH BURN OCCURS
                  
 # Check the test file to make sure it worked, but the correct file's name is the one listed below to be read in...
IBS.data <- read.table("29.nearest",header=TRUE) # IBS data with Z-scores and genetic distance

## Find individual who is most related to everyone else in sample by finding lowest average genetic distance ind

IID.meandist <- tapply(IBS.data$MIN_DST, INDEX=IBS.data$IID,FUN=mean,na.rm=TRUE) # this gives average genetic difference for all individuals, the one with the lowest score is centroid

center.ID <- names(IID.meandist[(which(IID.meandist==min(IID.meandist)))]) # This is the subject ID of the person who is most simililar in terms of IBS to everyone else in the sample

center.index <- which(IBS.data$IID==center.ID) # Indeces in IBS data of center ID rows

center.IBS.data <- IBS.data[center.index,]

hist(center.IBS.data$MIN_DST,1000, xlab="Genetic Distance From Centroid", main="Distribution of IBS relatedness") # histogram of relatedness in sample

## Determine 3 standard deviation cut off distance value and remove bottom percentile from sample:

sort.index <- order(center.IBS.data$MIN_DST)
sorted.IBS.data <- center.IBS.data[sort.index,]
final.IBS.data <- sorted.IBS.data[1:(nrow(sorted.IBS.data)*.997),]

#######################################################################
# Now Perform again to find correct center after pruning and re-iterate:

## Find individual who is most related to everyone else in sample by finding lowest average genetic distance ind
IBS.data <- final.IBS.data

IID.meandist <- tapply(IBS.data$MIN_DST, INDEX=IBS.data$IID,FUN=mean,na.rm=TRUE) # this gives average genetic difference for all individuals, the one with the lowest score is centroid

center.ID <- names(IID.meandist[(which(IID.meandist==min(IID.meandist)))]) # This is the subject ID of the person who is most simililar in terms of IBS to everyone else in the sample

center.index <- which(IBS.data$IID==center.ID) # Indeces in IBS data of center ID rows

center.IBS.data <- IBS.data[center.index,]

hist(center.IBS.data$MIN_DST,1000, xlab="Genetic Distance From Centroid", main="Distribution of IBS relatedness") # histogram of relatedness in sample

## Determine 3 standard deviation cut off distance value and remove bottom percentile from sample:

sort.index <- order(center.IBS.data$MIN_DST)
sorted.IBS.data <- center.IBS.data[sort.index,]
final.IBS.data <- sorted.IBS.data[1:(nrow(sorted.IBS.data)*.997),]

hist(final.IBS.data$MIN_DST,1000, xlab="Genetic Distance From Centroid", main="Distribution of IBS relatedness") # histogram of relatedness in sample

up.bound <- mean(final.IBS.data$MIN_DST)+(3*sd(final.IBS.data$MIN_DST))
low.bound <- mean(final.IBS.data$MIN_DST)-(3*sd(final.IBS.data$MIN_DST))

final.drop.index <- which((final.IBS.data$MIN_DST>low.bound) & (final.IBS.data$MIN_DST<up.bound))
# Inspect visually:
hist(final.IBS.data$MIN_DST[final.drop.index],1000, xlab="Genetic Distance From Centroid", main="Distribution of IBS relatedness") # histogram of relatedness in sample

Final.IBS.data <- final.IBS.data[final.drop.index,]

# Keep List of individuals who pass stratification test:

keep.strat.list <- rbind(cbind(as.character(Final.IBS.data[1,1]),as.character(Final.IBS.data[1,2])),cbind(as.character(Final.IBS.data[,6]),as.character(Final.IBS.data[,7]))) # list of FID and IID for individuals to keep

write.table(keep.strat.list,file="29keep.strat.list.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

# Use input file with SNPs in LD to generate file that contains correct individuals and correct SNPs
system("plink --bfile 22.step.drop.sex --keep 29keep.strat.list.txt --make-bed --out 29.step.strat")

################################################################################
#############################################Generate genetic relatedness matrix
# This is genetic relatedness of just individuals in the sample that pass stratification and have phenotype data:
system("plink --bfile 29.step.strat --keep ps.keep --make-bed --out 29.step.ps")

system("nohup gcta --bfile 29.step.ps --make-grm --autosome --out 29.step.HE &")

###############################################
# read in data for HE regression:
w.HE <- read.table('29.step.HE.grm',header=FALSE)

grm.matrix.HE <- matrix(NA,nrow=length(unique(w.HE$V1)),ncol=length(unique(w.HE$V1)))

# Generate Bottom Half of Matrix
for (i in 1:nrow(w.HE)){
grm.matrix.HE[w.HE$V1[i],w.HE$V2[i]] <- w.HE$V4[i]   
}#

# remove values in GR matrix greater than .025 to ensure only unrelated CELLS remain
remove.index <- which(grm.matrix.HE>.025)
grm.matrix.HE[remove.index] <- NA
#
###############################################
################## Generate Phenotype Matrix
fam.data <- read.table("29.step.strat.fam",header=FALSE)
names(fam.data) <- c('FID','IID','V3','V4','SEX','NA')
# read in phenotype data:
risk.scores <- read.table("Risk_Scores.csv",header=FALSE,skip=1,sep=",")
names(risk.scores)<- c("IID")
merged.data <- merge(risk.scores,fam.data,"IID") # merge fam data and risk data
# Now merge with list of individuals that passed stratification test, to ensure individuals that have phenotype data and pass strat test are in list:
strat.dat <- as.data.frame(read.table("29keep.strat.list.txt",header=FALSE))
names(strat.dat) <- c('FID','IID')
ps.merge <- merge(merged.data,strat.dat,"IID")

#############################################
# extract phenotype data:
PHE.pre <- ps.merge[,c(7,1,6)]
names(PHE.pre) <- c('FID','IID','RAW.PHE')

# Now run natural log transform on phenotype
e <- (1+(1/100000))^100000
PHE.pre[,3] <- log(PHE.pre[,3],base=e)
# remove single outlier:
PHE.pre[which(PHE.pre[,3]< (-10)),3] <- mean(PHE.pre[,3])
# assign to PHE (phenotype variable)
PHE <- PHE.pre
names(PHE) <- c('FID','IID','log.PHE')

####### Now run linear model with covariates being controlled for and use residuals as PHE for HE regression:

# Step 1: read in covariates and merge with data: (Age, Batch, Sex, Stratification)

# read in batch:

batch <- read.table("batch.covar",header=FALSE) # read in batch data
names(batch) <- c('FID','IID','BatchD1','BatchD2') # give batch col's names

# read in sex:

sex <- read.table("29.step.ps.fam",header=FALSE)[,c(1,2,5)] # read in batch data
names(sex) <- c('FID','IID','sex') # give sex col's names

# extract age:

age <- ps.merge[,c(7,1,5)]
names(age) <- c('FID','IID','age')

# read in stratification data generated using GCTA and R:

strat <- read.table("MDS.20.HE",header=FALSE)
names(strat) <- c('FID','IID','E1','E2','E3','E4','E5','E6','E7','E8','E9','E10','E11','E12','E13','E14','E15','E16','E17','E18','E19','E20')

# Step 2: merge Phenotype and covariates and run regression, use residuals in HE regression:
PHE.batch <- merge(PHE,batch,'IID')[,c(1,3,5,6)]
PHE.b.sex <- merge(PHE.batch,sex,'IID')[,c(1:4,6)]
PHE.b.s.age <- merge(PHE.b.sex,age,'IID')[,c(1:5,7)]
PHE.b.s.a.strat <- merge(PHE.b.s.age,strat,'IID')
PHE.tot <- PHE.b.s.a.strat

# remove NA's in columns and fill in with average of that column:

for (i in 1:ncol(PHE.tot)){
  index.NA <- which(is.na(PHE.tot[,i])==TRUE)
  PHE.tot[index.NA,i] <- mean(PHE.tot[,i],na.rm=TRUE)
}

# Now run model to control for covariates:
PHE.covar.model <- lm(log.PHE ~ BatchD1 + BatchD2 + sex + age + E1 + E2 + E3 + E4 + E5 + E6 + E7 + E8 + E9 + E10 + E11 + E12 + E13 + E14 + E15 + E16 + E17 + E18 + E19 + E20, data=PHE.tot)

PHE.r <- resid(PHE.covar.model) # get residuals from model with covariates

PHE.data.final <- as.data.frame(cbind(PHE.tot[,1],PHE.r))
names(PHE.data.final) <- c('IID','PHE.r')

###Now write out data that will be used as phenotype data

write.table(PHE.data.final,file="PHE.data.HE",quote=FALSE,row.names=FALSE,col.names=FALSE)
PHE.r <- read.table("PHE.data.HE",header=FALSE)[,2]

PHE.table <- cbind(PHE.r,PHE.r)
# generate covariance matrix from phenotype data:

PHE.matrix <- as.matrix(dist(PHE.table))
PHE2.matrix <- PHE.matrix^2 # make sure to square phenotypic difference as is done in HE regression

# Now remove top half of Matrix because only single pairwise comparison's are desired
for (i in 1:nrow(PHE2.matrix)){
PHE2.matrix[w.HE$V2[i],w.HE$V1[i]] <- NA # this uses convenient index provided by GCTA .grm input file    
}#

PHE.G <- cbind(as.matrix(as.vector(PHE2.matrix)),as.matrix(as.vector(grm.matrix.HE))) # transform matrices into a single rows for binding and regression analysis

# change mode of data to single instead of double so half as much memory is used:

mode(PHE.G) <- c('single')

# Now perform HE regression: (also time it)

HE.model <- lm(PHE.G[,1]~PHE.G[,2]) # h^2 =  -Beta1/2

Intercepts <- rep(0,1000)
Betas <- rep(0,1000) # generate vector to be filled with resampled betas
# Now resample with replacement
for (i in 1:1000){
index.sample <- sample(nrow(PHE2.matrix),replace=TRUE)# generate sampling index            
PHE.G.resample <- cbind(as.matrix(as.vector(PHE2.matrix[index.sample,])),as.matrix(as.vector(grm.matrix.HE[index.sample,]))) # transform matrices into a single rows for binding and regression analysis
HE.model.resample <- lm(PHE.G.resample[,1]~PHE.G.resample[,2])
Intercepts[i] <- HE.model.resample$coefficients[1]
Betas[i] <- HE.model.resample$coefficients[2]
} # each loop takes 70.8 seconds; total time for 1k iterations roughly 20 hrs

#  write out Va, Vp, and H2 data:
I.frame <- as.data.frame(Intercepts)
B.frame <- as.data.frame(Betas)
H2.frame <- as.data.frame(Betas/(-1*Intercepts)) # create data frame of heritability
names(H2.frame) <- c('H2')

write.table(B.frame,file="VA.dist",quote=FALSE,row.names=FALSE,col.names=TRUE)
write.table((-1*I.frame),file="VP.dist",quote=FALSE,row.names=FALSE,col.names=TRUE)
write.table(H2.frame,file="H2.dist",quote=FALSE,row.names=FALSE,col.names=TRUE)

# Based on results, empirical p=0.026
# Standard Error is just sd(Betas)
#compute 95% CI=Stat +- 1.96*(sd(Stat)) 
#*************** add related individuals to see the additional effect of rare variants
PHE.v <- as.matrix(read.table("VP.dist",header=TRUE))
G.v <- as.matrix(read.table("VA.dist",header=TRUE))
H2.v <- as.matrix(read.table("H2.dist",header=TRUE))

# Take standard deviation of parameter estimates (PHE.v, G.V, H2.v) to find empirical standard error


########################################### Experimental Garbage Below This Point:
##################################################################################
# Test effect of using PCA's of grm matrix as predictors of phenotype:

# generate grm matrix with top and bottom:

full.grm.matrix <- matrix(0,nrow=length(unique(w.HE$V1)),ncol=length(unique(w.HE$V1)))

# Generate Bottom Half of Matrix
for (i in 1:nrow(w.HE)){
full.grm.matrix[w.HE$V1[i],w.HE$V2[i]] <- w.HE$V4[i]   
}#
# Generate Top Half of Matrix
for (i in 1:nrow(w.HE)){
full.grm.matrix[w.HE$V2[i],w.HE$V1[i]] <- w.HE$V4[i]   
}#

PHE.dist <- 1-PHE.matrix # pairwise difference in phenotype matrix

PHE.e <- cmdscale(PHE.dist,k=10,eig=TRUE)
plot(PHE.e$eig,type='b',main='Scree Plot of Eigenvectors',xlab='dimension',ylab='relative variance explained')

eig.PHE1 <- PHE.e$points[,1]

# how well does eigenvector 1 of PHE similarity predict Phenotype after having controlled for all covariates?

eig.model <- lm(PHE.r ~ eig.PHE1) # 90% variance explained in phenotype when predicted by first eigenvector of pairwise phenotypic similarity matrix

grm.dist <- max(full.grm.matrix,na.rm=TRUE)-full.grm.matrix # genetic relatedness pairwise distance

GEN.e <- cmdscale(grm.dist,k=20,eig=TRUE)

G <- GEN.e$points

# look at scree plot of effect for each eigen vector:
plot(GEN.e$eig,type='b',main='Scree Plot of Eigenvectors',xlab='dimension',ylab='relative variance explained')

# Examine predictive ability of eigenvectors:

pred.model <- lm(phe.r ~ G[,1] + G[,2] + G[,3] + G[,4] + G[,5] + G[,6] + G[,7] + G[,8] + G[,9] + G[,10] + G[,11] + G[,12] + G[,13] + G[,14] + G[,15] + G[,16] + G[,17] + G[,18] + G[,19] + G[,20])

pred.model2 <- lm(eig.PHE1 ~ G[,6]) # barely significant and explains basically no variance
