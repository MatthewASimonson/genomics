

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

MDS <- read.table("ft5.shellfish_MDS.tevecs",header=FALSE)
names(MDS) <- c('C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14','C15','C16','C17','C18','C19','C20')
sample <- read.table("ft5.shellfish_MDS.sample",header=TRUE)
names(sample) <- c('IID','IID2')

MDS.samp <- cbind.data.frame(sample[2:nrow(sample),],MDS)

samp.used <-  read.table("MESA.clean.FINAL.fam",header=FALSE)
names(samp.used) <- c('FID','IID')
MDS.used <- merge(MDS.samp,samp.used,by='IID')

batch.covar <- read.table("batch.covar",header=FALSE)
names(batch.covar) <- c('IID','SEX','AGE','BATCH')
total.covars <- merge(batch.covar,MDS.used,by='IID')

covar.file <- cbind.data.frame(total.covars$IID,total.covars$SEX,total.covars$AGE,total.covars$BATCH,total.covars[,9:28])
names(covar.file) <-  c('IID','SEX','AGE','BATCH','C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14','C15','C16','C17','C18','C19','C20')
# write out covariate files:

write.table(covar.file,file="MESA.clean.FINAL.covar",quote=FALSE,row.names=FALSE,col.names=FALSE)

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

# First convert SNPs to rs format and flip strand of SNPs listed as negative:

bim <- read.table("5.step.unrelated.bim",header=FALSE)
snp.annot <- read.table("GenomeWideSNP_6.na32.annot.csv",header=TRUE,fill=TRUE,sep=",")
dup.drops <- as.data.frame(snp.annot$Probe.Set.ID[which(duplicated(snp.annot$dbSNP.RS.ID)==TRUE)])
rs.conversion <- cbind.data.frame(snp.annot$Probe.Set.ID,snp.annot$dbSNP.RS.ID)

write.table(rs.conversion,file="rs.conversion",row.names=FALSE,quote=FALSE,col.names=FALSE)
write.table(dup.drops,file="dup.drops",row.names=FALSE,quote=FALSE,col.names=FALSE)

strand.info <- cbind.data.frame(snp.annot$dbSNP.RS.ID,snp.annot$Strand,snp.annot$Strand.Versus.dbSNP)

flip.index <- which(snp.annot$Strand.Versus.dbSNP=="reverse")
flip.snps <- snp.annot$Probe.Set.ID[flip.index]

write.table(flip.snps,file="flip.snps",row.names=FALSE,quote=FALSE,col.names=FALSE)

system("plink --bfile 5.step.unrelated --flip flip.snps --exclude dup.drops --make-bed --out 6.step.flpd") # first flip strand
system("plink --bfile 6.step.flpd --update-map rs.conversion --update-name --make-bed --out 6.step.rs") # now convert to rs

sample.snps <- read.table("6.step.rs.bim",header=FALSE)
write.table(sample.snps$V2,file="sample.snps",row.names=FALSE,quote=FALSE,col.names=FALSE)

# only include SNPs that are common to hapmap and data set to expedite merge:
system("plink --bfile hapmap_270_r22_b36_fwd --extract sample.snps --make-bed --out hapmap_sample")

common.snps <- read.table("hapmap_sample.map",header=FALSE)
write.table(common.snps$V2,file="common.snps",row.names=FALSE,quote=FALSE,col.names=FALSE)

system("plink --bfile 6.step.rs --extract common.snps --make-bed --out 6.step.common")

# Now merge with HapMap and detect strand problems and positional issues:
system("plink --bfile hapmap_sample --bmerge 6.step.common.bed 6.step.common.bim 6.step.common.fam --make-bed --out strand.check1")

# flip hapmap SNPs:

system("plink --bfile hapmap_sample --flip strand.check1.missnp --make-bed --out hapmap_sample_flipped1")


# Run GWAS with hapmap CEU and Sample CEU ancestry to detect final strand issues:
system("plink --bfile hapmap_sample_flipped1 --bmerge 6.step.common.bed 6.step.common.bim 6.step.common.fam --make-bed --out flip.detect")
system("plink --bfile hapmap_sample_flipped1 --exclude flip.detect.missnp --make-bed --out hapmap_sample_flipped2")
system("plink --bfile 6.step.common --exclude flip.detect.missnp --make-bed --out 6.step.common2")
system("plink --bfile hapmap_sample_flipped2 --bmerge 6.step.common2.bed 6.step.common2.bim 6.step.common2.fam --make-bed --out flip.gwas")

fam <- read.table("flip.gwas.fam",header=FALSE)
ids <- fam[,c(1,2)]
names(ids) <- c('FID','IID')
race <- read.table("MESA_race.txt",header=TRUE,sep="\t")
hap.race <- read.table("hapmap.pop",header=FALSE)
names(hap.race) <- c('FID','IID','RACE')
total.race <- rbind(race,hap.race[,c(2,3)])
ids.race <- merge(ids,total.race,by="IID") 

white.index <- which((ids.race$RACE=="CEU") | (ids.race$RACE=="1"))
ids.white <- ids.race[white.index,c(2,1,3)]

hap.pheno <- ids.white
names(hap.pheno) <- c('FID','IID','HAPMAP')
hap.pheno[which(hap.pheno$HAPMAP=="CEU"),3] <- 2
write.table(hap.pheno,file="hap.pheno",row.names=FALSE,quote=FALSE,col.names=FALSE)
system("plink --bfile flip.gwas --pheno hap.pheno --assoc --out flip.gwas.assoc") 

# Read in association results and flip SNPs with significant associations:

flip.assoc <- read.table("flip.gwas.assoc.assoc",header=TRUE)
drop.snps <- as.data.frame(as.character(flip.assoc$SNP[which(flip.assoc$P<0.05)]))
write.table(drop.snps,file="hap.drop",row.names=FALSE,quote=FALSE,col.names=FALSE)
system("plink --bfile flip.gwas --exclude hap.drop --recode --out 7_step_shellfish")

# NOTE: TRY AND CONVERT FROM BINARY TO PED/MAP USING GTOOL B4 RUNNING
system("gtool -P --ped 7_step_shellfish.ped --map 7_step_shellfish.map --og shellfish_MDS.gen --os shellfish_MDS.sample")

# Generate IBS matrix of genetic distance to be used for covariates:
system("shellfish.py --pca --numpcs 20 --maxprocs 17 --file shellfish_MDS --out 7_shell_mds") ##

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

# Read in MDS data:

MDS <- read.table("7_shell_mds.tevecs",header=FALSE)
names(MDS) <- c('C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14','C15','C16','C17','C18','C19','C20')

# Plot MDS data:

sample.index <- 181:nrow(MDS)
JPT.index <- 1:60
YRI.index <- 121:180
CEU.index <- 61:120

# Plot MDS dimensions 1 & 2 to show genetic relatedness of sample subjects relative to Hapmap
plot(MDS$C1[c(sample.index,YRI.index,JPT.index)], MDS$C2[c(sample.index,YRI.index,JPT.index)],main="MDS Genetic Relatedness of Individuals",xlab="Dimension1",ylab="Dimension2")
points(MDS$C1[sample.index],MDS$C2[sample.index],col="blue") # sample points are blue
points(MDS$C1[CEU.index],MDS$C2[CEU.index],col="green") # CEU Points are green
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
IND.dat <- read.table("6.step.rs.fam",header=FALSE) # read and merge IID's
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

## Merging Data-Sets:

# read in SNPs from dataset with lowest count:

bim <- read.table("MESA.clean.FINAL.bim",header=FALSE)
write.table(bim[,2],file="common.snps",quote=FALSE,row.names=FALSE,col.names=FALSE)

# Run the following code on all data sets so all have same SNPs:
system("plink --bfile 11.step.hap --extract common.snps --make-bed --out CARDIA.clean.FINAL")

# Merge all datasets:
CARDIA.clean.FINAL.bed CARDIA.clean.FINAL.bim CARDIA.clean.FINAL.fam
MESA.clean.FINAL.bed MESA.clean.FINAL.bim MESA.clean.FINAL.fam 

system("plink --bfile ARIC.clean.FINAL --merge-list merge.txt --make-bed --out MERGE.clean.FINAL")


