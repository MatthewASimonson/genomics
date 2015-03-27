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
