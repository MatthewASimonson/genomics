

ROH.matrix <- function(data.sets){# FUNCTION IS USED TO GENERATE A MAP OF ROH ACROSS GENOME IN MEGA-BASE BINS MERGED WITH PHENOTYPE DATA AND COVARIATES
# NOTE: all files relevent to any data-set must have same name before extension
# such as: 'x.hom', 'x.fam' and 'x.mds'

#############################################################
# Read in data-set .mds file and store data for covariates: #
#############################################################

mds.append <- matrix(NA,nrow=1,ncol=23) # generate matrix to store and append .mds files that are read in

for (i in 1:length(data.sets)){ # read in each .fam file and append to variable
temp.mds <- read.table(paste('',data.sets[i],'.mds',sep=""),header=TRUE)
temp.mds$IID <- paste(temp.mds[,1],temp.mds[,2],sep="")
temp.mds <- as.matrix(temp.mds)             
mds.append <- rbind(mds.append,temp.mds)
}# end loop i
mds.append <- mds.append[2:nrow(mds.append),] # remove empty top row of merged .fam files
mds.data <- as.data.frame(mds.append[,c(2,4:23)])

###############################################
# Read in data-set .fam files and store data: #
###############################################

fam.append <- matrix(NA,nrow=1,ncol=7) # generate matrix to store and append .fam files that are read in

for (i in 1:length(data.sets)){ # read in each .fam file and append to variable
temp.fam <- read.table(paste('',data.sets[i],'.fam',sep=""),header=FALSE)
temp.fam[,2] <- paste(temp.fam[,1],temp.fam[,2],sep="")
temp.fam <- as.matrix(temp.fam)
set.col <- matrix(rep(data.sets[i],nrow(temp.fam)),nrow=nrow(temp.fam),ncol=1)
temp.fam <- cbind(temp.fam,set.col) # merge
fam.append <- rbind(fam.append,temp.fam)
}# end loop i
fam.append <- fam.append[2:nrow(fam.append),] # remove empty top row of merged .fam files

keep.data <- as.data.frame(fam.append[,c(2,5,6,7)]) # keep relevent data from .fam files
names(keep.data) <- c('IID','SEX','PHE','DATA') # label IID, sex, and phenotype columns

###############################################
# Read in data-set .hom files and store data: #
###############################################

hom.append <- matrix(NA,nrow=1,ncol=13) # generate matrix to store and append .fam files that are read in
names(hom.append) <- c('FID','IID','PHE','CHR','SNP1','SNP2','POS1','POS2','KB','NSNP','DENSITY','PHOM','PHET')

for (i in 1:length(data.sets)){ # read in each .hom file and append to variable
temp.hom <- read.table(paste('',data.sets[i],'.hom',sep=""),header=TRUE)
temp.hom$IID <- paste(temp.hom[,1],temp.hom[,2],sep="")
temp.hom <- as.matrix(temp.hom)
hom.append <- rbind(hom.append,temp.hom)
}# end loop i
hom.append <- as.data.frame(hom.append[2:nrow(hom.append),]) # remove empty top row of merged .hom files and convert to data frame for sorting

order.index <- order(hom.append[,2],hom.append$CHR,hom.append$POS1) # sort .hom file by indivual, chromosome, BP start of run
hom.append <- hom.append[order.index,]

##############################################
# Generate storage matrix with 3200 Columns: #
##############################################

##################################
# MAP KEY: each bin containts 1Mb
#------------------------
# CHROMOSOME  |   BINS  |
#-------------+----------
# CHR 1:       1:247   
# CHR 2:       248:492
# CHR 3:       493:693
# CHR 4:       694:886
# CHR 5:       887:1068
# CHR 6:       1069:1240
# CHR 7:       1241:1400
# CHR 8:       1401:1548
# CHR 9:       1549:1686
# CHR 10:      1687:1822
# CHR 11:      1823:1958
# CHR 12:      1959:2091
# CHR 13:      2092:2205
# CHR 14:      2206:2312
# CHR 15:      2313:2414
# CHR 16:      2415:2505
# CHR 17:      2506:2588
# CHR 18:      2589:2666
# CHR 19:      2667:2731
# CHR 20:      2732:2796
# CHR 21:      2797:2844
# CHR 22:      2845:2895
#################################

num_inds <- length(unique(as.character(hom.append$IID))) # total number of subjects with runs
map.matrix <- matrix(0, nrow=num_inds,ncol=2895)

inds <- unique(as.character(hom.append$IID))
chrom.constant <- c(0,248,493,694,887,1069,1241,1401,1549,1687,1823,1959,2092,2206,2313,2415,2506,2589,2667,2732,2797,2845)

# Fill map matrix:

 for(i in 1:num_inds){ # execute for each subject
  for(h in 1:22){ # execute for each chromosome
  ind.rows <- which(hom.append$IID==inds[i]) # get index of rows for specific individual from hom.append
    for (j in ind.rows){# fill bins for specific indivual in individuals' row in map.matrix using 'if' switches for specific chromosome
      if (as.numeric(as.character(hom.append$CHR[j]))==h){ # SET SPECIFIC BIN RANGE FOR CHROMOSOME
        column.location1 <- round(as.numeric(as.character(hom.append[j,7]))/1000000) + chrom.constant[h]
        column.location2 <- round(as.numeric(as.character(hom.append[j,8]))/1000000) + chrom.constant[h]
        for(k in (column.location1:column.location2)){
         if (map.matrix[i,k]==0){
          map.matrix[i,k] <- map.matrix[i,k] +1 # ensuring max value of 1 in all bins
        } # end if statement  
           }# end loop k
          }# end if statement 1
        } # end loop j
      } # end loop h
    } # end loop i
#

#####################################################
# Merge ROH data with individual ID and Covariates: #
#####################################################

IID.matrix <- matrix(inds,nrow=num_inds,ncol=1)
IID.w.runs <- as.data.frame(cbind(IID.matrix,map.matrix))
names(IID.w.runs) <- c('IID')

keep.data2 <- merge(keep.data,mds.data,by='IID') # merge covariates with .fam data
ROH.matrix <- merge(keep.data2,IID.w.runs,by='IID') # This ROH.matrix excludes individuals that had no runs, they must be merged now:

no_run.index <- which(is.na(match(keep.data2[,1],IID.w.runs[,1]))) # index in keep.data2[,1] of subjects who have no runs
no_runs.subjects <- keep.data2[no_run.index,] # list of subjects with no runs, bind them to the bottom of IID.w.runs
fill.matrix <- matrix(0,nrow=nrow(no_runs.subjects),ncol=2895)

no.run.matrix <- cbind(no_runs.subjects,fill.matrix)
names(ROH.matrix) <- names(no.run.matrix)

ROH.FINAL <- rbind(ROH.matrix,no.run.matrix) # merge rows of subjects with and without runs (This is complete data)
return(ROH.FINAL)
} # end function ROH.matrix
#

###################################################################################################
# Generate Runs information for data sets seperately then merge data together and write out file: #
###################################################################################################

# Process data sets
data.sets <- c('ab.Final','bon.Final','bulg.Final','carwtc.Final','cat2.Final','dk.Final','dub.Final','edi.Final','mgs2.Final','muc.Final','port.Final','sw1.Final','sw2.Final','top3.Final','ucl.Final','ucla.Final','zhh.Final')

ROH.FINAL <- ROH.matrix(data.sets) # generate ROH.matrix using function from specified data sets

write.table(ROH.FINAL,file="ROH.matrix.data",quote=FALSE,row.names=FALSE,col.names=TRUE)#

##############################################################################
# FUNCTION TO GENERATE LOGISTIC REGRESSION DATA AND RUN PERMUTATION TESTING: #
##############################################################################

# Inputs include names of data sets in char vector, data matrix with runs, and number of cores to use
# Ouputs include 2 files: ROH.map.data and ROH.map.permute, containing results of oringinal model in first file and permutation results in second

ROH.permute <- function(data.sets,data.matrix.name,PROCESSES,permutes){
#################
# read in data: #
#################

ROH.FINAL <- read.table(data.matrix.name,header=TRUE,fill=TRUE,colClasses=c(rep('factor',4),rep('numeric',2915)))

#######################################################################
# Generate indeces of sex within dataset in ROH.FINAL for permutation #
#######################################################################

### GENERATE DATA STRUCTURE TO START AND STOP LOCATIONS OF SEX WITHIN DATA SET
data_set.sex.index <- matrix(NA, nrow=2,ncol=2*length(data.sets))
rownames(data_set.sex.index) <- c('male','female')
colnames(data_set.sex.index) <- c(paste(data.sets[1:length(data.sets)],'_start',sep=""),paste(data.sets[1:length(data.sets)],'_stop',sep=""))

### FIRST SORT DATA MATRIX BY DATA SET THEN BY SEX WITHIN EACH DATA SET:
order.index <- order(ROH.FINAL$DATA, ROH.FINAL$SEX)
ROH.FINAL <- ROH.FINAL[order.index,]

### SECOND FIND START AND STOP INDECES OF EACH SEX NESTED WITHIN DATA SET:
for (i in 1:length(data.sets)){
  set.index <- which(ROH.FINAL$DATA==data.sets[i]) # get index of data set
  data_set.sex.index[1,i] <- min(set.index[which(ROH.FINAL$SEX[set.index]==1)]) # start index for males in dataset i
  data_set.sex.index[1,i+length(data.sets)] <- max(set.index[which(ROH.FINAL$SEX[set.index]==1)]) # stop index for males in dataset i
  data_set.sex.index[2,i] <- min(set.index[which(ROH.FINAL$SEX[set.index]==2)]) # start index for females in dataset i
  data_set.sex.index[2,i+length(data.sets)] <- max(set.index[which(ROH.FINAL$SEX[set.index]==2)]) # stop index for females in dataset i
} # end loop i

### RE-ARRANGE DATA STRUCTURE SO MORE EASILY VIEWED:
order.index2 <- order(colnames(data_set.sex.index))
data_set.sex.index <- data_set.sex.index[,order.index2]

####################################################################
# Loop through permutations of Case/Control Status Using Milticore #
####################################################################
#Load packages & important functions
# NOTE: THIS IS CURRENTLY SET FOR MULTI-NODE
require(foreach)
require(doMPI) # This must be loaded for 'foreach' to work
cl <- startMPIcluster(PROCESSES) 
registerDoMPI(cl) # DoMC must be registered so the foreach tasks are executed in parallel instead of sequentially; NUMBER OF CORES CAN BE SET HERE

corecheck <-getDoParWorkers() # assign number of cores being used by foreach
write.table(corecheck,file="corecheck",quote=FALSE,row.names=FALSE,col.names=FALSE)# write out number of cores being used by foreach
loop <- foreach(z=1:permutes) %dopar% {# Start permutation loop z (first iteration writes un-permuted output)

############################################################################
# Re-arrange data structure for each permutation based on nested structure #
############################################################################

resample.index <- 1:nrow(ROH.FINAL) # generate vector of sample size
for (i in seq(from=1,to=length(data.sets)*2,2)){
  resample.index[data_set.sex.index[1,i]:data_set.sex.index[1,i+1]] <- sample(resample.index[data_set.sex.index[1,i]:data_set.sex.index[1,i+1]],replace=TRUE) # reshuffle male index for given data set
  resample.index[data_set.sex.index[2,i]:data_set.sex.index[2,i+1]] <- sample(resample.index[data_set.sex.index[2,i]:data_set.sex.index[2,i+1]],replace=TRUE) # reshuffle female index for given data set
  }  # end loop i        
 
#####################################################
# Generate matrix to fill with relevent information #
#####################################################

Beta0 <- rep(0,(ncol(ROH.FINAL)-24))
Beta1 <- rep(0,(ncol(ROH.FINAL)-24))
Beta0_P_val <- rep(0,(ncol(ROH.FINAL)-24))
Beta1_P_val <- rep(0,(ncol(ROH.FINAL)-24))
Beta0_StandardErr <- rep(0,(ncol(ROH.FINAL)-24))
Beta1_StandardErr <- rep(0,(ncol(ROH.FINAL)-24))
Tot_resid_deviance <- rep(0,(ncol(ROH.FINAL)-24))
df <- rep(0,(ncol(ROH.FINAL)-24))

logistic.matrix <- rbind(Beta0,Beta1,Beta0_P_val,Beta1_P_val,Beta0_StandardErr,Beta1_StandardErr,Tot_resid_deviance,df)

###########################################################
# loop through each Mb window and run logistic regression #
###########################################################

# First loop doesn't permute:
if (z>1){
  ROH.FINAL <- ROH.FINAL[resample.index,] # turns resampling on at all values of z loop greater than 1
} # end if statement turning re-sampling on/off

for( j in 25:27){# (ncol(ROH.FINAL))
logistic.model <- glm(PHE ~ as.numeric(ROH.FINAL[,j]) + SEX + DATA + as.numeric(C1) + as.numeric(C2) + as.numeric(C3) + as.numeric(C4) + as.numeric(C5) + as.numeric(C6) + as.numeric(C7) + as.numeric(C8) + as.numeric(C9) + as.numeric(C10), family=binomial("logit"),data=ROH.FINAL)

logistic.matrix[1,j-24] <- as.numeric(logistic.model$coefficients[1]) # Beta0 
logistic.matrix[2,j-24] <- as.numeric(logistic.model$coefficients[2]) # Beta1 
logistic.matrix[3,j-24] <- as.numeric(summary(logistic.model)$coefficients[1,4])# Beta0_P_val
logistic.matrix[4,j-24] <- as.numeric(summary(logistic.model)$coefficients[2,4])# Beta1_P_val
logistic.matrix[5,j-24] <- as.numeric(summary(logistic.model)$coefficients[1,2])# Beta0_StandardErr
logistic.matrix[6,j-24] <- as.numeric(summary(logistic.model)$coefficients[2,2])# Beta1_StandardErr 
logistic.matrix[7,j-24] <- logistic.model$deviance # Tot_resid_deviance
logistic.matrix[8,j-24] <- summary(logistic.model)$df[2] # df
} # end loop j

if (z<2){
write.table(logistic.matrix,file="ROH.map.data",quote=FALSE,row.names=TRUE,col.names=FALSE)# write out un-permuted results
}else{
write.table(logistic.matrix,file="ROH.map.perm",quote=FALSE,row.names=TRUE,col.names=FALSE,append=TRUE)# append permuted results to file with columns labeled so grep can seperate iterations
}

} # end permutation loop z
} # end function ROH.permute

##############################################################################
# GENERATE LOGISTIC REGRESSION DATA AND RUN PERMUTATION TESTING: #
##############################################################################
# NOTE: MAKE SURE TO CHANGE LOOP z to 1001 AND FIX LOOP J

# Function Inputs:
 data.sets <- c('ab.Final','bon.Final','bulg.Final','carwtc.Final','cat2.Final','dk.Final','dub.Final','edi.Final','mgs2.Final','muc.Final','port.Final','sw1.Final','sw2.Final','top3.Final','ucl.Final','ucla.Final','zhh.Final')

data.matrix.name <- c('ROH.matrix.data') # make sure this file is in directory

cores <- 4 #Number of cores to use

permutations <- 1000 # number of permutations

ROH.permute(data.sets,data.matrix.name,cores,permutations)
