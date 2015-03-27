
#GOAL: Now that the genotype files are ready, this script prunes the data and runs ROH analyses

# 4 levels used:
# 1) Moderate LD pruning - 35 SNP threshold - No het allowance
# 2) Moderate LD pruning - 45 SNP threshold - No het allowance
# 3) Moderate LD pruning - 50 SNP threshold - No het allowance
# 4) Light LD pruning - 65 SNP threshold - No het allowance


#LD pruning results will be stored in the excel file QC_simonsonROH_2012.xls 
#NOTE: Script accidently overwrote pruning results log file with pruning log file


#==================================
#====== Merged ARIC, MESA, CARDIA =
#==================================

# On STATGEN
setwd("/home/simonsom/ROH_pathway/TOTAL/WHI")

#Variables
data_set <- "WHI.clean.final"

#Remove MAF < 0.05
system(paste("plink --bfile ",data_set," --maf .05 --make-bed --out ",data_set,"ROH_maf",sep=""))

#Light LD-pruning
system(paste("plink --bfile ",data_set,"ROH_maf --indep 50 5 10 --out ",data_set,"ROH_lite",sep=""))
system(paste("plink --bfile ",data_set,"ROH_maf --extract ",data_set,"ROH_lite.prune.in --make-bed --out ",data_set,"ROH_lite",sep=""))
system(paste("cp ",data_set,"ROH_lite.log ",data_set,"ROH_lite.log",sep=""))

# --- Genotype Info in pruned data

#light LD-pruning
system(paste("plink --bfile ",data_set,"ROH_lite --missing --out ",data_set,"ROH_lite",sep="")) #missingness
system(paste("plink --bfile ",data_set,"ROH_lite --het --out ",data_set,"ROH_lite",sep="")) #Inb. coeff
system(paste("plink --bfile ",data_set,"ROH_lite --freq --out ",data_set,"ROH_lite",sep=""))
frq <- read.table("WHI.clean.finalROH_lite.frq",header=TRUE)
list <- frq[,c(2,3)]
one <- rep(1,nrow(frq))
list <- cbind(list,one)
write.table(list,"list",row.names=FALSE,col.names=FALSE,quote=FALSE)
system(paste("plink --bfile ",data_set,"ROH_lite --score list --out ",data_set,"ROH_lite",sep="")) #minor allele load

# 4) Lite LD pruning - 65 SNP threshold - No het allowance (ROH_lite_snp65)

# Parameters:
# --homozyg-window-snp 65
# --homozyg-snp 65
# --homozyg-window-missing 3
# --homozyg-window-het 0
# --homozyg-window-threshold 0.03

# --homozyg-window-kb 0
# --homozyg-kb 0
# --homozyg-gap 2000
# --homozyg-density 300

# --homozyg-group
# --homozyg-match .95 
  system(paste("nohup plink --bfile ",data_set,"ROH_lite --homozyg-window-snp 65 --homozyg-snp 65 --homozyg-window-missing 3 --homozyg-window-het 0 --homozyg-window-threshold 0.03 --homozyg-window-kb 0 --homozyg-kb 0 --homozyg-gap 2000 --homozyg-density 300 --homozyg-group --homozyg-match .95 --out ",data_set,"ROH_lite_snp65 &",sep=""))









