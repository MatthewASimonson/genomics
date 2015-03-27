##############################################################################
# FUNCTION TO GENERATE LOGISTIC REGRESSION DATA AND RUN PERMUTATION TESTING: #
##############################################################################

# INPUT: 1.) names of data sets in char vector, 2.) data matrix with runs, 3.) number of cores to use, 4.) number of permutations to run

# OUTPUT: 1.) A file is written: ROH.map.data, containing the initial results of logistic regression, with one row for 8 examined variables, 2.) 2 files are written containing data from specified number of permutations for Beta0 and Beta1

ROH.permute <- function(data.sets,data.matrix.name,PROCESSES,permutes){
#################
# read in data: #
#################

ROH.FINAL <- read.table(data.matrix.name,header=TRUE,fill=TRUE,colClasses=c(rep('factor',4),rep('integer',2875)))

#######################################################################
# Generate indeces of sex within dataset in ROH.FINAL for permutation #
#######################################################################

### GENERATE DATA STRUCTURE TO STORE START AND STOP LOCATIONS OF SEX WITHIN DATA SET
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
# Loop through using Milticore #
####################################################################
#Load packages & important functions
# NOTE: THIS IS CURRENTLY SET FOR MULTI-NODE
#FOR JANUS: NEXT 4 LINES

require(foreach)
require(doMPI) # This must be loaded for 'foreach' to work
cl <- startMPIcluster(PROCESSES) 
registerDoMPI(cl) # DoMC must be registered so the foreach tasks are executed in parallel instead of sequentially; NUMBER OF CORES CAN BE SET HERE

loop <- foreach(z=1:(permutes+1)) %dopar% {# Start permutation loop z (first iteration writes un-permuted output)
  set.seed(z) # SEED IS SET TO ITERATOR
##################################################################################################################
# Generate index to re-arrange Phenotype column of data structure for each permutation based on nested structure #
##################################################################################################################

resample.index <- 1:nrow(ROH.FINAL) # generate vector of sample size
for (i in seq(from=1,to=length(data.sets)*2,2)){
  resample.index[data_set.sex.index[1,i]:data_set.sex.index[1,i+1]] <- sample(resample.index[data_set.sex.index[1,i]:data_set.sex.index[1,i+1]],replace=FALSE) # reshuffle male index for given data set
  resample.index[data_set.sex.index[2,i]:data_set.sex.index[2,i+1]] <- sample(resample.index[data_set.sex.index[2,i]:data_set.sex.index[2,i+1]],replace=FALSE) # reshuffle female index for given data set
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
Adjusted_Rsquared <- rep(0,(ncol(ROH.FINAL)-24))
DF <- rep(0,(ncol(ROH.FINAL)-24))

logistic.matrix <- rbind(Beta0,Beta1,Beta0_P_val,Beta1_P_val,Beta0_StandardErr,Beta1_StandardErr,Adjusted_Rsquared,DF)

###########################################################
# loop through each Mb window and run logistic regression #
###########################################################
ROH.LOGISTIC <- ROH.FINAL # assign data to temporary frame

# First loop doesn't permute:
if (z>1){ 
  ROH.LOGISTIC[,3] <- ROH.LOGISTIC[resample.index,3] # turns resampling on at all values of z loop greater than 1 in the PHE COLUMN ONLY
} # end 'if' statement turning re-sampling on/off

for( j in 25:(ncol(ROH.LOGISTIC))){# loop through each megabase | (ncol(ROH.LOGISTIC))
logistic.model <- lm(as.numeric(PHE) ~ ROH.LOGISTIC[,j] + SEX + DATA + as.numeric(C1) + as.numeric(C2) + as.numeric(C3) + as.numeric(C4) + as.numeric(C5) + as.numeric(C6) + as.numeric(C7) + as.numeric(C8) + as.numeric(C9) + as.numeric(C10), data=ROH.LOGISTIC)

logistic.matrix[1,j-24] <- as.numeric(logistic.model$coefficients[1]) # Beta0 
logistic.matrix[2,j-24] <- as.numeric(logistic.model$coefficients[2]) # Beta1 
logistic.matrix[3,j-24] <- as.numeric(summary(logistic.model)$coefficients[1,4])# Beta0_P_val
logistic.matrix[4,j-24] <- as.numeric(summary(logistic.model)$coefficients[2,4])# Beta1_P_val
logistic.matrix[5,j-24] <- as.numeric(summary(logistic.model)$coefficients[1,2])# Beta0_StandardErr
logistic.matrix[6,j-24] <- as.numeric(summary(logistic.model)$coefficients[2,2])# Beta1_StandardErr 
logistic.matrix[7,j-24] <- as.numeric(summary(logistic.model)$adj.r.squared) # Adjusted R^2
logistic.matrix[8,j-24] <- as.numeric(summary(logistic.model)$df[2]) # degrees of freedom
} # end loop j

# 
if (z==1){
write.table(logistic.matrix,file="ROH.map.data",quote=FALSE,row.names=TRUE,col.names=FALSE)# write out un-permuted results
}else{
write.table(logistic.matrix,file=paste("temp.ROH",z-1,"",sep=""),quote=FALSE,row.names=TRUE,col.names=FALSE)# write each iteration to a temporary file
}

} # end permutation loop z; end multicore

###################################################
# Read in data-set temp.ROH files and merge data #
###################################################

roh.append <- matrix(NA,nrow=1,ncol=2896) # generate matrix to store and append .mds files that are read in

for (k in 1:permutes){ # read in each .fam file and append to variable
temp.roh <- read.table(paste("temp.ROH",k,"",sep=""),header=FALSE)
temp.roh <- as.matrix(temp.roh)             
roh.append <- rbind(roh.append,temp.roh)
system(paste("rm temp.ROH",k,"",sep="")) # remove temp file after reading in 
}# end loop k

roh.append <-roh.append[2:nrow(roh.append),] # remove empty top row of merged .fam files
write.table(roh.append,file="ROH.map.permute",quote=FALSE,row.names=FALSE,col.names=FALSE,append=FALSE)#
#

} # end function ROH.permute

##############################################################################
# GENERATE LOGISTIC REGRESSION DATA AND RUN PERMUTATION TESTING: #
##############################################################################

# Function Inputs:
data.sets <- c('ab','bon','bulg','carwtc','cat','dk','dub','edi','mgs','muc','port','sw1','sw2','top3','ucl','ucla','zhh')

data.matrix.name <- c('ROH.matrix.data') # make sure this file is in directory

cores <- 100 #Number of cores to use

permutations <- 1000 # number of permutations

ROH.permute(data.sets,data.matrix.name,cores,permutations)

