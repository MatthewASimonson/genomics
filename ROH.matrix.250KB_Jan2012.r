# INPUTS:
#  character vector containing list of data sets, the functions expects a .hom, .fam, and .mds file for each data set all with the same name minus different extensions.
# WHOLE GENOME DATA IS EXPECTED!!!


# OUTPUTS: a data frame containing rows of subjects and columns for each megabase across all autosomal regions of the genome (28950 total); if a subject containts a run within a given 100KB bin (contained in each column), their respective row will contain a '1', else a '0' for that megabase bin. The first 24 columns contain ID, sex, phenotype, dataset, and 20 PCA's.

# A key file is also generated specifying chromosome and base-position for each bin


ROH.matrix <- function(data.sets){# FUNCTION IS USED TO GENERATE A MAP OF ROH ACROSS GENOME IN MEGA-BASE BINS MERGED WITH PHENOTYPE DATA AND COVARIATES
# NOTE: all files relevent to any data-set must have same name before extension
# such as: 'x.hom', 'x.fam' and 'x.mds'

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

#################################################
# Get base position indeces of each chromosome: #
#################################################

hom.data <- hom.append
order.hom <- order(hom.data$CHR,hom.data$POS2) # order by chromosome and stop indeces
hom.data <- hom.data[order.hom,]

bp.chrom.count <- rep(NA,22)
bp.chrom.pre <- rep(NA,22)
bp.chrom.indeces <- rep(NA,22) # This object stores the last base position for each chromosome after generated below

for(i in 1:22){
bp.chrom.count[i]  <- max(which(as.numeric(hom.data$CHR)==i))
bp.chrom.pre[i] <- as.numeric(as.character(hom.data$POS2[bp.chrom.count[i]]))
bp.chrom.indeces[i] <- sum(bp.chrom.pre[1:i])
}

chrom.constant <- rep(0,22) # generate object to store stop indeces of each chromosome
base.positions.per.bin <- 250000 # set base positions per bin here

for (i in 2:22){
chrom.constant[i] <- round((bp.chrom.indeces[i-1]/base.positions.per.bin)) + 5
}

# USEFUL FOR PLOTTING START AND STOP INDECES:


for(i in 1:22){
bp.chrom.count[i]  <- min(which(as.numeric(hom.data$CHR)==i)) # index or row with last position of chromosome
bp.chrom.pre[i] <- as.numeric(as.character(hom.data$POS2[bp.chrom.count[i]])) # first base position of each chromosome
}

chrom.constant.start <- rep(NA,22)
for (i in 1:22){
chrom.constant.start[i] <- round((bp.chrom.pre[i]/base.positions.per.bin))+chrom.constant[i] # first bin that has a run in each chromosome
}

chrom.constant.end <- c(chrom.constant[2:22],(round((bp.chrom.indeces[22]/base.positions.per.bin)) + 5)) # last bin with run in each chromosome

#############################
# Generate storage matrix : #
#############################

num_inds <- length(unique(as.character(hom.append$IID))) # total number of subjects with runs
map.matrix <- matrix(0, nrow=num_inds,ncol=(round((bp.chrom.indeces[22]/base.positions.per.bin))+10)) # generate storage matrix to be filled

inds <- unique(as.character(hom.append$IID))

chrom.label <- rep(NA,(round((bp.chrom.indeces[22]/base.positions.per.bin)) + 10))
chrom.bp <- rep(NA,(round((bp.chrom.indeces[22]/base.positions.per.bin)) + 10))

map.matrix.key <- rbind(chrom.label,chrom.bp)

map.matrix.key.fill <- matrix(0, nrow=num_inds,ncol=round((bp.chrom.indeces[22]/base.positions.per.bin)) + 10) # generate storage matrix to be filled
map.matrix.chrom.fill <- matrix(0, nrow=num_inds,ncol=round((bp.chrom.indeces[22]/base.positions.per.bin)) + 10) # generate storage matrix to be filled
# Fill map matrix and generate map key:

 for(i in 1:num_inds){ # execute for each subject
   print(i) # output subject being started
  for(h in 1:22){ # execute for each chromosome
  ind.rows <- which(hom.append$IID==inds[i]) # get index of rows for specific individual from hom.append
    for (j in ind.rows){# fill bins for specific indivual in individuals' row in map.matrix using 'if' switches for specific chromosome
      if (as.numeric(as.character(hom.append$CHR[j]))==h){ # SET SPECIFIC BIN RANGE FOR CHROMOSOME (if statement 1)
        column.location1 <- round(as.numeric(as.character(hom.append[j,7]))/base.positions.per.bin) + chrom.constant[h] 
        column.location2 <- round(as.numeric(as.character(hom.append[j,8]))/base.positions.per.bin) + chrom.constant[h]
        # Fill Chrom Matrix:
        map.matrix.chrom.fill[i,(column.location1:column.location2)] <- as.numeric(as.character(hom.append[j,4]))      
        # Fill Key Matrix:
         bp1 <- as.numeric(as.character(hom.append[j,7]))
         bp2 <- as.numeric(as.character(hom.append[j,8]))
         map.matrix.key.fill[i,(column.location1:column.location2)] <- c(seq(from=bp1,to=bp2,by=base.positions.per.bin),bp2)[1:((column.location2-column.location1)+1)]
        # Fill Runs Matrix:
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

# Find Base Position for each window and place in map.matrix.key:

for (i in 1:ncol(map.matrix.key.fill)){
map.matrix.key[2,i] <- max(map.matrix.key.fill[,i],na.rm=TRUE) # max value from last run position across all runs in each bin
print(i)
}

# Find Chromosome position for each window and place in map.matrix.key:

for (i in 1:ncol(map.matrix.chrom.fill)){
map.matrix.key[1,i] <- max(map.matrix.chrom.fill[,i],na.rm=TRUE)
print(i)
}

###########################################
# Fill in gaps at chromosome junctions:   #
###########################################

# Fill residual indeces:
   res.ind <- matrix(NA,nrow=2,ncol=22)
# FIX FUNCTION:
chrom.junction.fill <- function(map.matrix.key,res.ind,base.positions.per.bin){
for(i in 1:22){
# fill matrix that holds index for empty junctions between chromosomes to be filled
  res.ind[1,i] <- min(which(map.matrix.key[1,]==i))
  res.ind[2,i] <- max(which(map.matrix.key[1,]==i))
}# end res.ind fill loop

for (i in 1:21){
# fill empty cells before first run starts on chromosome 1:
  if ((res.ind[1,1])>1){ # if start gaps, fill them by subtracting from known start
    for(j in 1:((res.ind[1,1])-1)){
      map.matrix.key[1,(res.ind[1,1])-j] <- map.matrix.key[1,(res.ind[1,1])] # fix chrom row in chrom 1
      map.matrix.key[2,(res.ind[1,1])-j] <- (map.matrix.key[2,(res.ind[1,1])]-base.positions.per.bin) # fix bp row in chrom 1
    } # end loop j
  }# end if statement for chromosome 1
  
  # detect if a gap exists:
  if((res.ind[1,i+1])>(res.ind[2,i]+1)){ # if first filled bin in next chromosome is not directly after last filled bin of previous then ...
   fill.range <- ((res.ind[2,i])+1):(res.ind[1,i+1]-1)# through the range of empty bins
   # detect if first bin of next chromosome has fewer bp's than minimum; if so, fill range starting at end of last run
      if((map.matrix.key[2,(res.ind[1,i+1])])<base.positions.per.bin){
        map.matrix.key[1,fill.range] <- map.matrix.key[1,(res.ind[2,i])]
        map.matrix.key[2,fill.range] <- seq(from=(map.matrix.key[2,(res.ind[2,i])])+base.positions.per.bin,by=base.positions.per.bin,length.out=length(fill.range))
      }else{ # if first bin of next chromosome has more bp's then minimum:
        bin.count.next <- round(map.matrix.key[2,(res.ind[1,i+1])]/base.positions.per.bin) # get a count of number of bins to fill from next chromosome
        if(bin.count.next>=length(fill.range)){ # detects if all empty bins belong to next chromosome
        map.matrix.key[1,fill.range] <- map.matrix.key[1,(res.ind[1,i+1])]
        map.matrix.key[2,fill.range] <- seq(from=(map.matrix.key[2,(res.ind[1,i+1])])-base.positions.per.bin,by=-base.positions.per.bin,length.out=length(fill.range))
        }# end bin.count.next loop
      }
  }# end if for detect gap
} # end loop i

return(map.matrix.key) 
} # end function chrom.junction.fill

ROH.map.key <- chrom.junction.fill(map.matrix.key,res.ind,base.positions.per.bin)
ROH.map.key <- chrom.junction.fill(ROH.map.key,res.ind,base.positions.per.bin) # RUN FUNCTION A SECOND TIME TO PICK UP ANY MISSED DURING FIRST PASS, WHICH WILL HAVE BEEN MADE DETECTABLE BY FIRST PASS

#############################
# Fill in remaining gaps:   #
#############################

gap.index <- which(as.numeric(ROH.map.key[1,])==0)

for (i in gap.index){
  ROH.map.key[1,i] <- ROH.map.key[1,i-1]
  ROH.map.key[2,i] <- ROH.map.key[2,i-1]+base.positions.per.bin
}

write.table(ROH.map.key,file="ROH.matrix.key",quote=FALSE,row.names=TRUE,col.names=FALSE)# Write out Key

#####################################################
# Merge ROH data with individual ID and Covariates: #
#####################################################

IID.matrix <- matrix(inds,nrow=num_inds,ncol=1)
IID.w.runs <- as.data.frame(cbind(IID.matrix,map.matrix))
names(IID.w.runs) <- c('IID')

keep.data2 <- keep.data
ROH.matrix <- merge(keep.data2,IID.w.runs,by='IID') # This ROH.matrix excludes individuals that had no runs, they must be merged now:
#

no_run.index <- which(is.na(match(keep.data2[,1],IID.w.runs[,1]))) # index in keep.data2[,1] of subjects who have no runs
no_runs.subjects <- keep.data2[no_run.index,] # list of subjects with no runs, bind them to the bottom of IID.w.runs
fill.matrix <- matrix(0,nrow=nrow(no_runs.subjects),ncol=ncol(map.matrix))

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

data.sets <- c("training10")

ROH.FINAL <- ROH.matrix(data.sets) # generate ROH.matrix using function from specified data sets

write.table(ROH.FINAL,file="ROH.matrix.data",quote=FALSE,row.names=FALSE,col.names=TRUE)#
#
