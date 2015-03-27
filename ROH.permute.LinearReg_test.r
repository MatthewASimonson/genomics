############################################################
# FOLLOWING WERE WRITTEN BY MATTHEW A. SIMONSON 2011-2012 #
############################################################

##############################################################################
# FUNCTION TO GENERATE LINEAR REGRESSION DATA AND RUN PERMUTATION TESTING:   #
##############################################################################

# INPUT: 1.) names of data sets in char vector, 2.) data matrix with runs in 1MB bins, 3.) number of cores to use, 4.) number of permutations to run

# OUTPUT: 1.) A file is written: ROH.map.data, containing the initial results of logistic regression, with one row for 8 examined variables, 2.) 2 files are written containing data from specified number of permutations for Beta0 and Beta1

ROH.permute <- function(data.sets,data.matrix.name,PROCESSES,permutes){
#################
# read in data: #
#################

ROH.FINAL <- read.table(data.matrix.name,header=TRUE,fill=TRUE,colClasses=c(rep('character',4),rep('integer',2877)))

#################################################################
# Generate indeces of each dataset in ROH.FINAL for permutation #
#################################################################

### GENERATE DATA STRUCTURE TO STORE START AND STOP LOCATIONS OF SEX WITHIN DATA SET
data_set.index <- matrix(NA, nrow=1,ncol=2*length(data.sets))
colnames(data_set.index) <- c(paste(data.sets[1:length(data.sets)],'_start',sep=""),paste(data.sets[1:length(data.sets)],'_stop',sep=""))

### FIRST SORT DATA MATRIX BY DATA SET THEN BY SEX WITHIN EACH DATA SET:
order.index <- order(ROH.FINAL$DATA)
ROH.FINAL <- ROH.FINAL[order.index,]

### SECOND FIND START AND STOP INDECES OF EACH SEX NESTED WITHIN DATA SET:
for (i in 1:length(data.sets)){
  set.index <- which(ROH.FINAL$DATA==data.sets[i]) # get index of data set
  data_set.index[1,i] <- min(set.index) # start index for dataset i
  data_set.index[1,i+length(data.sets)] <- max(set.index) # stop index for dataset i
} # end loop i

################################
# Loop through using Milticore #
################################
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
  sample.cols <- ncol(data_set.index)/2
for (i in 1:sample.cols){
  resample.index[data_set.index[1,i]:data_set.index[1,i+sample.cols]] <- sample(resample.index[data_set.index[1,i]:data_set.index[1,i+sample.cols]],replace=FALSE) # reshuffle index for given data set
  }  # end loop i        

#####################################################
# Generate matrix to fill with relevent information #
#####################################################

Beta0 <- rep(0,(ncol(ROH.FINAL)-4))
Beta1 <- rep(0,(ncol(ROH.FINAL)-4))
Beta0_P_val <- rep(0,(ncol(ROH.FINAL)-4))
Beta1_P_val <- rep(0,(ncol(ROH.FINAL)-4))
Beta0_StandardErr <- rep(0,(ncol(ROH.FINAL)-4))
Beta1_StandardErr <- rep(0,(ncol(ROH.FINAL)-4))
Adjusted_Rsquared <- rep(0,(ncol(ROH.FINAL)-4))
DF <- rep(0,(ncol(ROH.FINAL)-4))

linear.matrix <- rbind(Beta0,Beta1,Beta0_P_val,Beta1_P_val,Beta0_StandardErr,Beta1_StandardErr,Adjusted_Rsquared,DF)

###########################################################
# loop through each Mb window and run linear regression #
###########################################################
ROH.LINEAR <- ROH.FINAL # assign data to temporary frame

# First loop doesn't permute:
if (z>1){ 
  ROH.LINEAR$PHE <- ROH.LINEAR$PHE[resample.index] # turns resampling on at all values of z loop greater than 1 in the PHE COLUMN ONLY
} # end 'if' statement turning re-sampling on/off

for( j in 5:(ncol(ROH.LINEAR))){# loop through each megabase | 5:(ncol(ROH.LINEAR))
linear.model <- lm(as.numeric(as.character(ROH.LINEAR$PHE)) ~ ROH.LINEAR[,j])

write.table(j,file="status",quote=FALSE,row.names=FALSE,col.names=FALSE)# Write out permutation status

linear.matrix[1,j-4] <- as.numeric(linear.model$coefficients[1]) # Beta0 
linear.matrix[2,j-4] <- as.numeric(linear.model$coefficients[2]) # Beta1 
linear.matrix[3,j-4] <- as.numeric(summary(linear.model)$coefficients[1,4])# Beta0_P_val

linear.matrix[5,j-4] <- as.numeric(summary(linear.model)$coefficients[1,2])# Beta0_StandardErr

if(is.na(as.numeric(linear.model$coefficients[2]))==TRUE){
linear.matrix[4,j-4] <- NA # Beta1_P_val
linear.matrix[6,j-4] <- NA # Beta1_StandardErr
}else{ # in case NA for Beta 1
linear.matrix[4,j-4] <- as.numeric(summary(linear.model)$coefficients[2,4])# Beta1_P_val
linear.matrix[6,j-4] <- as.numeric(summary(linear.model)$coefficients[2,2])# Beta1_StandardErr
}
linear.matrix[7,j-4] <- as.numeric(summary(linear.model)$adj.r.squared) # Adjusted R^2
linear.matrix[8,j-4] <- as.numeric(summary(linear.model)$df[2]) # degrees of freedom
} # end loop j

# 
if (z==1){
write.table(linear.matrix,file="ROH.map.data",quote=FALSE,row.names=TRUE,col.names=FALSE)# write out un-permuted results
}else{
write.table(linear.matrix,file=paste("temp.ROH",z-1,"",sep=""),quote=FALSE,row.names=TRUE,col.names=FALSE)# write each iteration to a temporary file
}

} # end permutation loop z; end multicore

###################################################
# Read in data-set temp.ROH files and merge data #
###################################################

roh.append <- matrix(NA,nrow=1,ncol=(ncol(ROH.FINAL)-3)) 

for (k in 1:permutes){ # read in each temp file
temp.roh <- matrix(scan(paste("temp.ROH",k,"",sep=""),what='character'),nrow=8,ncol=(ncol(ROH.FINAL)-3),byrow=TRUE)
temp.roh <- as.matrix(temp.roh)             
roh.append <- rbind(roh.append,temp.roh)
#system(paste("rm temp.ROH",k,"",sep="")) # remove temp file after reading in 
}# end loop k

roh.append <-roh.append[2:nrow(roh.append),] # remove empty top row of merged .fam files
write.table(roh.append,file="ROH.map.permute",quote=FALSE,row.names=FALSE,col.names=FALSE,append=FALSE)#
#

} # end function ROH.permute

##############################################################################
# GENERATE LINEAR REGRESSION DATA AND RUN PERMUTATION TESTING: #
##############################################################################

# Function Inputs:
#data.sets <- c('ROH_IQ_rd5')

data.matrix.name <- c('ROH.matrix.data') # make sure this file is in directory

cores <- 1 #Number of cores to use

permutations <- 0 # number of permutations

ROH.permute(data.sets,data.matrix.name,cores,permutations)
           
