# Framingham GWAS

setwd("/Net/statgen2/home/shared/km/Framingham/Merged_data/")

#########################################################################
# Step 1:
# Remove individuals with discrepencies between reported and genotypic sex

system("plink --bfile merge.NHLBI --check-sex --out 1.step.checksex")
sexchk <- read.table("1.step.checksex",header=TRUE)
hist(sexchk$F[sexchk$F<.5], breaks=20)

highF <-  sexchk[which(sexchk$STATUS!='OK'),]
discrepant <- highF[which(highF$F>.3),]

FID.IID <- cbind(as.character(discrepant$FID),as.character(discrepant$FID))

write.table(FID.IID,file="remove.list.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

# use system command to remove discrepant individuals from data set and write new data files

system("plink --bfile merge.NHLBI --remove remove.list.txt --make-bed --out 1.1.step.sex")

#################################################################################
#### Step 2:

# Remove individuals who are missing more than 5% of SNP calls
# First step is generating missingness data with PLINK

system("plink --bfile 1.1.step.sex --missing")

ind.miss <- read.table("plink.imiss",header=TRUE)
above.five <- which(ind.miss[,6]>.05)
high.missing <- ind.miss[above.five,]

# Create a text file with IID and FID for individuals with missingness above 5%

miss.dat <- cbind(as.character(high.missing$FID),as.character(high.missing$FID))

write.table(miss.dat,file="remove.miss.list.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

system("plink --bfile 1.1.step.sex --remove remove.miss.list.txt --make-bed --out 3.step.ind_miss")


#################################################################################
# Step 3: Make sure physical positions of SNPs is correct

#Check position of SNPs in Illumina annotation file and update position of all SNPs

# read in annotated affymetrix data and compare with previous SNP data:
#part 1: Sort data by chromosome

affy.annot <- read.table("GenomeWideSNP_5.na30.annot.csv", header=TRUE,sep=",", colClasses=c(rep("character",5),rep("NULL",22)))

sorted.affy.index <-order(affy.annot$Chromosome,affy.annot$Physical.Position)

sorted.affy <- affy.annot[sorted.affy.index,] 

# remove SNPs with no chromosome or physical marker info
keep.rows <- which(sorted.affy$Chromosome!='---')
sorted.affy <- sorted.affy[keep.rows,1:5]# remove after column 7
affy.snps <- sorted.affy[,1]
# read in Map file from original data

map <- read.table("merge.NHLBI.bim", header=FALSE)
snps <- as.vector(map[,2])

# match snps in old data with affy list

snps.2.compare <- match(snps,affy.snps)
relevent.affy <- sorted.affy[snps.2.compare,] # only relevent snps
diff.position <- which(map[,4!=relevent.affy[,5]) # indeces of snps with changed base postions
snp.diff <- abs(as.numeric(relevent.affy[diff.position,5])-map[diff.position,4]) # how large is the difference
large.diff.index <- which(snp.diff>100)
remove.snps <- as.data.frame(map[diff.position,2][large.diff.index])
write.table(remove.snps,file="diff.position.txt", quote=FALSE,sep="",row.names=FALSE,col.names=FALSE) # write list of snps to remove that have significantly different base positions

# Write new data file that excludes specified SNPs
system("plink --bfile merge.NHLBI.bim --exclude diff.position.txt --make-bed --out 4.step.annot")

##################################################################################
# Step 4: Drop SNPs with MAF <.05

system("plink --bfile 4.step.annot --maf 0.05 --write-snplist")

system("plink --bfile 4.step.annot --extract plink.snplist --make-bed --out 5.step.noMAF05")
##################################################################################
# Step 5: Drops SNPs with call rates <.98


system("plink --bfile 5.step.noMAF05 --missing")
snp.miss <- read.table("plink.lmiss",header=TRUE)
bad.fmiss.index <- which(snp.miss$F_MISS > .02)
bad.snp.list <- snp.miss[bad.fmiss.index,2] # create a list a SNPs with call rates <=.98
write.table(bad.snp.list,file="remove.highmiss.list.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

system("plink --bfile 5.step.noMAF05 --exclude remove.highmiss.list.txt --make-bed --out 6.step.snpcall")                       
##################################################################################
# Step 6:  Prune SNPs out of HWE

system("plink --bfile 6.step.snpcall --hardy")

HWE.data <- read.table("plink.hwe",header=TRUE)
out.HWE <- which(HWE.data[,9]<.0001)
prune.HWE.snps <- HWE.data[out.HWE,2]

write.table(prune.HWE.snps,file="remove.HWE.list.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

system("plink --bfile 6.step.snpcall --exclude remove.HWE.list.txt --make-bed --out 7.step.HWE")            
                                             
system("plink --bfile 7.step.HWE --remove drop.list.txt --keep keep.IND.nohap.txt --recode --out 7.step.pruned.geno")

##################################################################################
# Step: 7 Prune data for admixture
##################################################################################
#NEXT STEP:
                       
#fam.data <- read.table("merge.NHLBI.fam", header=FALSE) # read in fam file to determine number of subjects
#num.subjects <- nrow(fam.data)# each row in file is independent subject

#system(paste("plink --bfile merge.NHLBI --cluster --neighbour 1 ",num.subjects," &",sep='')) # takes a long time to run, about a week for 10000 subjects

# More needs to be added here to examine the z-scores and prune outliers
#
