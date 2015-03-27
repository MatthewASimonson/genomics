# Merge plink format files from ARIC, MESA, WHI, CARDIA, GENEVA_diabetes and MGS:

######################################################
# Make sure same SNPs are used across all data sets: #
######################################################
setwd("/home/simonsom/ROH_pathway/TOTAL/GENEVA_diab/TOTAL")
GENEVA.snps <- read.table("GENEVA.total.final.bim",header=FALSE) # GENEVA # 5445 subjects

setwd("/home/simonsom/ROH_pathway/TOTAL/WHI")
WHI.snps <- read.table("WHI.clean.final.bim",header=FALSE) # WHI # 3167 subjects

setwd("/home/simonsom/ROH_pathway/TOTAL")
AMC.snps <- read.table("MERGE.clean.FINAL.bim",header=FALSE) # ARIC, MESA, CARDIA # 11883

setwd("/home/simonsom/ROH_pathway/TOTAL/MGS/mgs_ea_recalled")
MGS.snps <- read.table("MGS.final.bim",header=FALSE) # MGS controls # 1443 subjects

# total: 21938
# split:
# discovery = ARIC, MESA, CARDIA, WHI = 15050 subjects
# replication = GENEVA, MGS = 6888 subjects

wg.snps <- merge(GENEVA.snps,WHI.snps,by="V2") 
