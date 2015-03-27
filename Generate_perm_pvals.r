######################################################################
# Generate empirical p-values for ROH data based on permutation test #
######################################################################

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
setwd('/STATGEN/home/simonsom/PGC.ROH/map/Logistic_test')

# Read in data with runs and covariate information:
runs.data <- read.table("ROH.matrix.data",header=TRUE,fill=TRUE,colClasses=c(rep('factor',4),rep('numeric',2915)))# runs start at col 25
runs.only <- runs.data[,25:ncol(runs.data)] # generate matrix without covariates

# Read in map and permute data:
map.data <- read.table("ROH.map.data",header=FALSE) # read in un-permuted results
permute.data <- read.table("ROH.map.permute",header=FALSE) # read in results of permutation

# get count of runs at each Mb position:
run.count <- rep(NA,2895)

for (i in 1:2895){
  run.count[i] <- sum(runs.only[,i]) # count number of runs in each bin
                    }

exclude.index <- which(run.count<2) # index of where fewer than 2 runs exist in the genome

# Fill indeces to be excluded with NA's that have too few runs

map.data[,exclude.index+1] <- NA
permute.data[,exclude.index+1] <- NA

# Store all Beta1's in structure:
Beta1.index <- which(permute.data$V1=='Beta1')
Beta1.data <- permute.data[Beta1.index,]

# Store all P-values for Beta1's in structure:
pvals.index <- which(permute.data$V1=='Beta1_P_val')
pval.data <- permute.data[pvals.index,]

# Generate objects to store lowest p-value indeces and selected Beta1's:

low_p.index <- rep(NA, 1000)
Beta1s <- rep(NA, 1000)

# get index of lowest p-value Beta from each iteration across genome (exclude NA bins)

for (i in 1:1000){
  # low_p.index[i] <- which(as.numeric(pval.data[i,2:ncol(pval.data)])==min(as.numeric(pval.data[i,2:ncol(pval.data)]),na.rm=TRUE)) # BUG IN 'WHICH' FUNCTION RETURNS MULTIPLE VALUES WHEN ONLY A SINGLE VALUE MATCHES THE MIN!!!!
  low_p.index[i] <- which.min(as.numeric(pval.data[i,2:ncol(pval.data)]))
  print(i)
}

# Double check to make sure the numbr of runs where lowest p-values are from looks correct:

run.count[low_p.index]

# Assign Beta's to structure:

for (i in 1:1000){
  Beta1s[i] <- Beta1.data[i,low_p.index[i]+1] # all the '+1's in indexing are due to the first column being the field name
}

# Examine 95% CI's

hist(abs(Beta1s),100)

one.tail <- sort(abs(Beta1s))

AUC <- sum(one.tail) # total area under curve
AUC.percentile <- rep(NA,1000)

for (i in 1:1000){
  AUC.percentile[i] <- 1-(sum(one.tail[1:i])/AUC)
}

plot(one.tail,AUC.percentile, xlab="Absolute Beta Value",ylab="Probability of Value",main="Probability Distribution of Beta's")

# Find 95% cutoff:

value.index <- which.max(which(AUC.percentile<.95))
one.tail.95 <- one.tail[value.index]

# Plot data:

plot(abs(as.numeric(map.data[2,2:2895])),type='l')
