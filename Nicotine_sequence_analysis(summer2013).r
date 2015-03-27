# Created by Matthew A. Simonson 6/20/2012:

###########################################################################
# Examine if number of variants is predictive of substance abuse phenotype#
###########################################################################

#########################
# Step 1: Organize Data #
#########################
# Load list of data:

setwd("/home/simonsom/SOLiD_new")

system("ls *.csv > files")
files <- as.matrix(read.table("files",header=FALSE))

SS <- list()
SS.names <- vector()
for(i in 1:nrow(files)){
  SS[[i]] <- read.csv(files[i,1],header=TRUE)
  SS.names[i] <- strsplit(files[i,1],"\\.")$V1[1]
  print(i)
}
names(SS) <- SS.names  


#install.packages("SKAT")
library("SKAT")

##############
# Run tests: #
##############
##############################################
# Using functions from SKAT library:   #
##############################################

# SKAT
# Notes:
#  To run SKAT-O, please set method=”optimal.adj” in the SKAT function,
# which provides a better type I error control in tail areas then method=”optimal”
#
#  You can use SKAT_CommonRare function to run a test for the combined effect of common
# and rare variants described in Ionita-Laza et al (-C and -A methods).
#
# THIS FUCNTION CAN ALSO BE USED FOR CONTINUOUS PHENOTYPES; LOOK INTO DETAILS
#
# Variable formatting:
#
# Z = a numeric genotype matrix of individuals and SNPs. Each row represents a different individual, and each column represents a different SNP marker.
# X a numeric matrix of covariates
# y.b a numeric vector of binary phenotypes
#
#

# Compute the P-value of SKAT with default Beta(1,25) Weights
# dichotomous trait with covariates:

X <- as.matrix(SS[[14]][2:ncol(SS[[14]])])
y.b <- as.matrix(SS[[8]][,2])

null <- SKAT_Null_Model(y.b ~ X, out_type="D") # calculate null
SKAT(Z, obj.b)$p.value

##############################################
# Using functions from AssotesteR library:   #
##############################################

#install.packages("AssotesteR")
library("AssotesteR")

# C-ALPHA
for(i in 
calp = CALPHA(CHRNA4.phe, CHRNA4.gen, perm=10000)
