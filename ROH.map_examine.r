# Examine map of autozygosity data:

setwd('/home/simonsom/PGC.ROH/map') # Set working directory

############################
# Read in ROH data files:  #
############################

# Read in values for un-permuted data:

runs.data <- read.table('ROH.map.data',header=FALSE, fill=TRUE)

permute.data <- read.table('ROH.map.permute',header=FALSE, fill=TRUE) 
################################
# Examine distribution of OR:  #
################################

B1.permute.index <- which(permute.data[,1]=='Beta1')
B1.permute <- permute.data[B1.permute.index,]

hist(B1.permute[,2],100)

############### Plot initial state of case/control status:
resample.index.pre <- 1:nrow(ROH.FINAL) # generate vector of sample size

###### Plot Sex:
plot(resample.index.pre,cex=.001, main='Resampling Ranges based on Dataset/Sex')

male.location <- which(as.numeric(ROH.FINAL[,2])==2)
male.points <- points(male.location,male.location,col='red',cex=.001)

female.location <- which(as.numeric(ROH.FINAL[,2])==3)
female.points <- points(female.location,female.location,col='blue',cex=.001)

abline(v=c(1,984,1420,2324,3204,3764,4340,6106,7752,8365,8554,9106,9502,9947,10634,11046,11286,14419,16620,17067,17408,17700,17970,18144,18305,18651,18924,19243,19541,20121,20556,21398,21898,22117))

abline(v=c(1,1420,3204,4340,7752,8554,9502,10634,11286,16620,17408,17970,18305,18924,19541,20556,21898,22279),col='green')

###### Plot Affection:
plot(resample.index.pre,cex=.001, main='Resampling Ranges based on Dataset/Sex of Case/Control status')

case.location <- which(as.numeric(ROH.FINAL[,3])==2)
case.points <- points(case.location,case.location,col='red',cex=.001)

control.location <- which(as.numeric(ROH.FINAL[,3])==1)
control.points <- points(control.location,control.location,col='blue',cex=.001)

abline(v=c(1,984,1420,2324,3204,3764,4340,6106,7752,8365,8554,9106,9502,9947,10634,11046,11286,14419,16620,17067,17408,17700,17970,18144,18305,18651,18924,19243,19541,20121,20556,21398,21898,22117))

abline(v=c(1,1420,3204,4340,7752,8554,9502,10634,11286,16620,17408,17970,18305,18924,19541,20556,21898,22279),col='green')

############# Re-order correctly:

### FIRST SORT DATA MATRIX BY DATA SET THEN BY SEX WITHIN EACH DATA SET:
order.index <- order(ROH.FINAL$DATA, ROH.FINAL$SEX)
ROH.FINAL <- ROH.FINAL[order.index,]

# Compare data at strange index to normal:

# Normal Permutation:

plot(ROH.FINAL[,25],cex=.3)

control.location <- which(as.numeric(ROH.FINAL[,3])==1)
control.points <- points(control.location,ROH.FINAL[control.location,25],col='blue',cex=.3)

case.location <- which(as.numeric(ROH.FINAL[,3])==2)
case.points <- points(case.location,ROH.FINAL[case.location,25],col='red',cex=.3)

abline(v=c(1,984,1420,2324,3204,3764,4340,6106,7752,8365,8554,9106,9502,9947,10634,11046,11286,14419,16620,17067,17408,17700,17970,18144,18305,18651,18924,19243,19541,20121,20556,21398,21898,22117))

abline(v=c(1,1420,3204,4340,7752,8554,9502,10634,11286,16620,17408,17970,18305,18924,19541,20556,21898,22279),col='green')

# Re-orderm just as in script

### FIRST SORT DATA MATRIX BY DATA SET THEN BY SEX WITHIN EACH DATA SET:
order.index <- order(ROH.FINAL$DATA, ROH.FINAL$SEX)
ROH.FINAL <- ROH.FINAL[order.index,]

# Weird Permutation:

weird.perm <- as.numeric(read.table('temp.perm105')) # this is a strange permutation

# order things using the weird permutation:

ROH.FINAL <- ROH.FINAL[weird.perm,]

plot(ROH.FINAL[,25],cex=.3)

control.location <- which(as.numeric(ROH.FINAL[,3])==1)
control.points <- points(control.location,ROH.FINAL[control.location,25],col='blue',cex=.3)

case.location <- which(as.numeric(ROH.FINAL[,3])==2)
case.points <- points(case.location,ROH.FINAL[case.location,25],col='red',cex=.3)

abline(v=c(1,984,1420,2324,3204,3764,4340,6106,7752,8365,8554,9106,9502,9947,10634,11046,11286,14419,16620,17067,17408,17700,17970,18144,18305,18651,18924,19243,19541,20121,20556,21398,21898,22117))

abline(v=c(1,1420,3204,4340,7752,8554,9502,10634,11286,16620,17408,17970,18305,18924,19541,20556,21898,22279),col='green')
#
