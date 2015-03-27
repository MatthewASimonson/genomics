# Flipping Strands:

################################################
# Convert 1000genomes data into binary format: #
################################################
system('/parkinsons/1000genome')

 for(l in 1:22){ # START LOOP l
        eval(parse(text=paste('system("plink --file EUR.20100804.chr',l,' --make-bed --out EUR.20100804.chr',l,'")',sep='')))
      } # END LOOP l

#######################
# CIDR:               #
#######################
setwd('/projects/KellerLab/simonson/parkinsons/geno_calls/CIDR/bim_bed_fam/phg000022.CIDR_PD.genotype-calls.HumanCNV370v1.v1.p1.c1.GRU.matrixfmt.genotype')

# Step 1: Update position of all SNPs so same as HapMap and flip SNPs based on Illumina strand list:
#########

map.data <- read.table("hapmap_r23a.bim",header=FALSE)
map.dat.sub <- map.data[,c(2,4)]

write.table(map.dat.sub,file="hapmap23a_pos.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)# write file with hapmap positions of SNPs

system("plink --bfile CIDR_CNV370_clean --update-map hapmap23a_pos.txt --flip CIDR.Illu.flip.txt --make-bed --out CIDR_CNV370_clean_update")

# Step 2: read .fam data and create alternate phenotype that represents hapmap or dataset
#########

hap.fam <- read.table("hapmapCEU_f.fam",header=FALSE) # only CEU hapmap will be used because LD is most similiar to current sample
dat.fam <- read.table("CIDR_CNV370_clean.fam",header=FALSE)

# generate data frame with 3 col's FID,IID,PHE; PHE refers to hapmap or not in this case

hap.dat <- cbind(as.character(hap.fam$V1),as.character(hap.fam$V2),rep('1',nrow(hap.fam)))
dat.dat <- cbind(as.character(dat.fam$V1),as.character(dat.fam$V2),rep('2',nrow(dat.fam)))

cidr.hap <- as.data.frame(rbind(hap.dat,dat.dat))

# write out phenotype file:

write.table(cidr.hap,file="cidr.hap.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')
 
# Step 3: merge hapmap CEU and CIDR data:
#########

system("plink --bfile hapmap_CEU_r23a --bmerge CIDR_CNV370_clean_update.bed CIDR_CNV370_clean_update.bim CIDR_CNV370_clean_update.fam --make-bed --out strand.CIDR.HM.merge")

# now flip: WARNING THIS OVERWRITES A PRVIOUS FILE
system("plink --bfile CIDR_CNV370_clean_update --flip strand.CIDR.HM.merge.missnp --make-bed --out CIDR_CNV370_clean_update")

# Merge attempt 2:

system("plink --bfile hapmap_CEU_r23a --bmerge CIDR_CNV370_clean_update.bed CIDR_CNV370_clean_update.bim CIDR_CNV370_clean_update.fam --make-bed --out strand.CIDR.HM.merge")

# Step 4: Use linkage-disequilibrium data to detect any remaining strand issues:
# NOTE: BE SURE TO USE '--pheno' COMMAND TO SPECIFY SAMPLE VS. HAPMAP SPLIT
#########

# First keep only those SNPs from CIDR

keep.snps <- read.table("CIDR_CNV370_clean_update.bim",header=FALSE)

write.table(keep.snps[,c(2)],file="snp.extract.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

system("plink --bfile strand.CIDR.HM.merge --extract snp.extract.txt --make-bed --out CIDR.HM.flip.extract")

# Second generate LD information:

system("plink --bfile CIDR.HM.flip.extract --pheno cidr.hap.txt --flip-scan --out CIDR.flip.detect")

#  Read in Flip data: 

flip.data <- read.table("CIDR.flip.detect.flipscan",header=TRUE,fill=TRUE)

#################################################
#KEY FOR FLIP-DATA:
#     CHR     Chromosome
#     SNP     SNP identifier for index SNP
#     BP      Base-pair position
#     A1      Minor allele code
#     A2      Major allele code
#     F       Allele frequency (A1 allele)
#     POS     Number of positive LD matches
#     R_POS   Average correlation of these 
#     NEG     Number of negative LD matches
#     R_NEG   Average correlation of these
#     NEGSNPS The SNPs showing negative correlation

# Examine how many SNPs have negative LD matches, if less than 5000, just drop these SNPs:

hist(flip.data$NEG[which(flip.data$NEG>0)],30,main='Number of Negative LD Correlation SNPs')
table(flip.data$NEG)

drop.snps <- as.data.frame(flip.data$SNP[which(flip.data$NEG>0)])

write.table(drop.snps,file="drop.flipscan.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

# Third: Drop or flip SNPs with uncertain strand direction (assumes not many):

system("plink --bfile CIDR.HM.flip.extract --exclude drop.flipscan.txt --make-bed --out CIDR.FINAL")

###############################################################################################
#######################
# NINDS_250:          #
#######################
setwd('/projects/KellerLab/simonson/parkinsons/geno_calls/NINDS/bim_bed_fam/250s')

# Step 1: Update position of all SNPs so same as HapMap and flip SNPs based on Illumina strand list:
#########

map.data <- read.table("hapmap_r23a.bim",header=FALSE)
map.dat.sub <- map.data[,c(2,4)]

write.table(map.dat.sub,file="hapmap23a_pos.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)# write file with hapmap positions of SNPs

system("plink --bfile NINDS_250clean --update-map hapmap23a_pos.txt --make-bed --out NINDS_250clean_update") # NOTE: DUE TO THE USE OF A CUSTOM GENOTYPING CHIP NO SNPS ARE PREEMPTIVELY FLIPPED BASED ON ILLUMINA DATA UNLIKE ALL OTHER DATA SETS USED

# Step 2: read .fam data and create alternate phenotype that represents hapmap or dataset
#########

hap.fam <- read.table("hapmapCEU_f.fam",header=FALSE) # only CEU hapmap will be used because LD is most similiar to current sample
dat.fam <- read.table("NINDS_250clean.fam",header=FALSE)

# generate data frame with 3 col's FID,IID,PHE; PHE refers to hapmap or not in this case

hap.dat <- cbind(as.character(hap.fam$V1),as.character(hap.fam$V2),rep('1',nrow(hap.fam)))
dat.dat <- cbind(as.character(dat.fam$V1),as.character(dat.fam$V2),rep('2',nrow(dat.fam)))

ninds250.hap <- as.data.frame(rbind(hap.dat,dat.dat))

# write out phenotype file:

write.table(ninds250.hap,file="ninds250.hap.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')

# Step 3: merge hapmap CEU and CIDR data:
#########

system("plink --bfile hapmap_CEU_r23a --bmerge NINDS_250clean_update.bed NINDS_250clean_update.bim NINDS_250clean_update.fam --make-bed --out strand.NINDS.HM.merge")

# now flip: WARNING THIS OVERWRITES A PRVIOUS FILE
system("plink --bfile NINDS_250clean_update --flip strand.NINDS.HM.merge.missnp --make-bed --out NINDS_250clean_update")

# Merge attempt 2:

system("plink --bfile hapmap_CEU_r23a --bmerge NINDS_250clean_update.bed NINDS_250clean_update.bim NINDS_250clean_update.fam --make-bed --out strand.NINDS.HM.merge")

# Step 4: Use linkage-disequilibrium data to detect any remaining strand issues:
# NOTE: BE SURE TO USE '--pheno' COMMAND TO SPECIFY SAMPLE VS. HAPMAP SPLIT
#########

# First keep only those SNPs from NINDS and exclude those with different chromosome positions:

keep.snps <- read.table("NINDS_250clean_update.bim",header=FALSE)

write.table(keep.snps[,c(2)],file="snp.extract.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

system("plink --bfile strand.NINDS.HM.merge --extract snp.extract.txt --make-bed --out NINDS250.HM.flip.extract")

system("plink --bfile NINDS250.HM.flip.extract --exclude NINDS_250.dc.exclude --make-bed --out NINDS250.HM.flip.extract") # Now exlcude chromosome location mistake SNPs

# Second generate LD information:
#START HERE
system("plink --bfile NINDS250.HM.flip.extract --pheno ninds250.hap.txt --flip-scan --out NINDS250.flip.detect")

#  Read in Flip data: 

flip.data <- read.table("NINDS250.flip.detect.flipscan",header=TRUE,fill=TRUE)

#################################################
#KEY FOR FLIP-DATA:
#     CHR     Chromosome
#     SNP     SNP identifier for index SNP
#     BP      Base-pair position
#     A1      Minor allele code
#     A2      Major allele code
#     F       Allele frequency (A1 allele)
#     POS     Number of positive LD matches
#     R_POS   Average correlation of these 
#     NEG     Number of negative LD matches
#     R_NEG   Average correlation of these
#     NEGSNPS The SNPs showing negative correlation

# Examine how many SNPs have negative LD matches, if less than 5000, just drop these SNPs:

hist(flip.data$NEG[which(flip.data$NEG>0)],30,main='Number of Negative LD Correlation SNPs')
table(flip.data$NEG)

drop.snps <- as.data.frame(flip.data$SNP[which(flip.data$NEG>0)])

write.table(drop.snps,file="drop.flipscan.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

# Third: Drop or flip SNPs with uncertain strand direction (assumes not many):

system("plink --bfile NINDS250.HM.flip.extract --exclude drop.flipscan.txt --make-bed --out NINDS250.FINAL")

###############################################################################################
#######################
# NINDS300v1:         #
#######################
setwd('/projects/KellerLab/simonson/parkinsons/geno_calls/NINDS/bim_bed_fam/300v1')

# Step 1: Update position of all SNPs so same as HapMap and flip SNPs based on Illumina strand list:
#########

map.data <- read.table("hapmap_r23a.bim",header=FALSE)
map.dat.sub <- map.data[,c(2,4)]

write.table(map.dat.sub,file="hapmap23a_pos.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)# write file with hapmap positions of SNPs

system("plink --bfile NINDS_300clean --update-map hapmap23a_pos.txt --make-bed --out NINDS_300clean_clean_update")

# Step 2: read .fam data and create alternate phenotype that represents hapmap or dataset
#########

hap.fam <- read.table("hapmapCEU_f.fam",header=FALSE) # only CEU hapmap will be used because LD is most similiar to current sample
dat.fam <- read.table("NINDS_300clean.fam",header=FALSE)

# generate data frame with 3 col's FID,IID,PHE; PHE refers to hapmap or not in this case

hap.dat <- cbind(as.character(hap.fam$V1),as.character(hap.fam$V2),rep('1',nrow(hap.fam)))
dat.dat <- cbind(as.character(dat.fam$V1),as.character(dat.fam$V2),rep('2',nrow(dat.fam)))

ninds3.hap <- as.data.frame(rbind(hap.dat,dat.dat))

# write out phenotype file:

write.table(ninds3.hap,file="ninds3.hap.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')

# Step 3: merge hapmap CEU and CIDR data:
#########

system("plink --bfile hapmap_CEU_r23a --bmerge NINDS_300clean_clean_update.bed NINDS_300clean_clean_update.bim NINDS_300clean_clean_update.fam --make-bed --out strand.NINDS3.HM.merge")

# now flip: WARNING THIS OVERWRITES A PRVIOUS FILE
system("plink --bfile NINDS_300clean_clean_update --flip strand.NINDS3.HM.merge.missnp --make-bed --out NINDS_300clean_clean_update")

# Merge attempt 2:

system("plink --bfile hapmap_CEU_r23a --bmerge NINDS_300clean_clean_update.bed NINDS_300clean_clean_update.bim NINDS_300clean_clean_update.fam --make-bed --out strand.NINDS3.HM.merge")

# Step 4: Use linkage-disequilibrium data to detect any remaining strand issues:
# NOTE: BE SURE TO USE '--pheno' COMMAND TO SPECIFY SAMPLE VS. HAPMAP SPLIT
#########

# First keep only those SNPs from CIDR

keep.snps <- read.table("NINDS_300clean_clean_update.bim",header=FALSE)

write.table(keep.snps[,c(2)],file="snp.extract.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

system("plink --bfile strand.NINDS3.HM.merge --extract snp.extract.txt --make-bed --out NINDS3.HM.flip.extract")

system("plink --bfile NINDS3.HM.flip.extract --exclude NINDS3.dc.exclude --make-bed --out NINDS3.HM.flip.extract") # Now exlcude chromosome location mistake SNPs

# Second generate LD information:

system("plink --bfile NINDS3.HM.flip.extract --pheno ninds3.hap.txt --flip-scan --out NINDS3.flip.detect")

#  Read in Flip data: 

flip.data <- read.table("NINDS3.flip.detect.flipscan",header=TRUE,fill=TRUE)

#################################################
#KEY FOR FLIP-DATA:
#     CHR     Chromosome
#     SNP     SNP identifier for index SNP
#     BP      Base-pair position
#     A1      Minor allele code
#     A2      Major allele code
#     F       Allele frequency (A1 allele)
#     POS     Number of positive LD matches
#     R_POS   Average correlation of these 
#     NEG     Number of negative LD matches
#     R_NEG   Average correlation of these
#     NEGSNPS The SNPs showing negative correlation

# Examine how many SNPs have negative LD matches, if less than 5000, just drop these SNPs:

hist(flip.data$NEG[which(flip.data$NEG>0)],30,main='Number of Negative LD Correlation SNPs')
table(flip.data$NEG)

drop.snps <- as.data.frame(flip.data$SNP[which(flip.data$NEG>0)])

write.table(drop.snps,file="drop.flipscan.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

# Third: Drop or flip SNPs with uncertain strand direction (assumes not many):

system("plink --bfile NINDS3.HM.flip.extract --exclude drop.flipscan.txt --make-bed --out NINDS3.FINAL")

###############################################################################################
#######################
# NINDS_CORIEL:       #
#######################
setwd('/projects/KellerLab/simonson/parkinsons/geno_calls/NINDS/bim_bed_fam/lng_coriell_pd')

# Step 1: Update position of all SNPs so same as HapMap and flip SNPs based on Illumina strand list:
#########

map.data <- read.table("hapmap_r23a.bim",header=FALSE)
map.dat.sub <- map.data[,c(2,4)]

write.table(map.dat.sub,file="hapmap23a_pos.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)# write file with hapmap positions of SNPs

system("plink --bfile NINDS_CORIELclean --update-map hapmap23a_pos.txt --make-bed --out NINDS_CORIELclean_clean_update")

# Step 2: read .fam data and create alternate phenotype that represents hapmap or dataset
#########

hap.fam <- read.table("hapmapCEU_f.fam",header=FALSE) # only CEU hapmap will be used because LD is most similiar to current sample
dat.fam <- read.table("NINDS_CORIELclean.fam",header=FALSE)

# generate data frame with 3 col's FID,IID,PHE; PHE refers to hapmap or not in this case

hap.dat <- cbind(as.character(hap.fam$V1),as.character(hap.fam$V2),rep('1',nrow(hap.fam)))
dat.dat <- cbind(as.character(dat.fam$V1),as.character(dat.fam$V2),rep('2',nrow(dat.fam)))

ninds3.hap <- as.data.frame(rbind(hap.dat,dat.dat))

# write out phenotype file:

write.table(ninds3.hap,file="ninds3.hap.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')

# Step 3: merge hapmap CEU and CIDR data:
#########

system("plink --bfile hapmap_CEU_r23a --bmerge NINDS_CORIELclean_clean_update.bed NINDS_CORIELclean_clean_update.bim NINDS_CORIELclean_clean_update.fam --make-bed --out strand.NINDS3.HM.merge")

# now flip: WARNING THIS OVERWRITES A PRVIOUS FILE
system("plink --bfile NINDS_CORIELclean_clean_update --flip strand.NINDS3.HM.merge.missnp --make-bed --out NINDS_CORIELclean_clean_update")

# Merge attempt 2:

system("plink --bfile hapmap_CEU_r23a --bmerge NINDS_CORIELclean_clean_update.bed NINDS_CORIELclean_clean_update.bim NINDS_CORIELclean_clean_update.fam --make-bed --out strand.NINDS3.HM.merge")

# Step 4: Use linkage-disequilibrium data to detect any remaining strand issues:
# NOTE: BE SURE TO USE '--pheno' COMMAND TO SPECIFY SAMPLE VS. HAPMAP SPLIT
#########

# First keep only those SNPs from CIDR

keep.snps <- read.table("NINDS_CORIELclean_clean_update.bim",header=FALSE)

write.table(keep.snps[,c(2)],file="snp.extract.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

system("plink --bfile strand.NINDS3.HM.merge --extract snp.extract.txt --make-bed --out NINDS3.HM.flip.extract")
#
#system("plink --bfile NINDS3.HM.flip.extract --exclude NINDS3.dc.exclude --make-bed --out NINDS3.HM.flip.extract") # Now exlcude chromosome location mistake SNPs

# Second generate LD information:

system("plink --bfile NINDS3.HM.flip.extract --pheno ninds3.hap.txt --flip-scan --out NINDS3.flip.detect")

#  Read in Flip data: 

flip.data <- read.table("NINDS3.flip.detect.flipscan",header=TRUE,fill=TRUE)

#################################################
#KEY FOR FLIP-DATA:
#     CHR     Chromosome
#     SNP     SNP identifier for index SNP
#     BP      Base-pair position
#     A1      Minor allele code
#     A2      Major allele code
#     F       Allele frequency (A1 allele)
#     POS     Number of positive LD matches
#     R_POS   Average correlation of these 
#     NEG     Number of negative LD matches
#     R_NEG   Average correlation of these
#     NEGSNPS The SNPs showing negative correlation

# Examine how many SNPs have negative LD matches, if less than 5000, just drop these SNPs:

hist(flip.data$NEG[which(flip.data$NEG>0)],30,main='Number of Negative LD Correlation SNPs')
table(flip.data$NEG)

drop.snps <- as.data.frame(flip.data$SNP[which(flip.data$NEG>0)])

write.table(drop.snps,file="drop.flipscan.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

# Third: Drop or flip SNPs with uncertain strand direction (assumes not many):

system("plink --bfile NINDS3.HM.flip.extract --exclude drop.flipscan.txt --make-bed --out NINDS3_CORIEL.FINAL")
