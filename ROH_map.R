# By: Matthew Simonson
# 10/25/10
#
setwd('/STATGEN/home/simonsom/PGC.ROH')

# Step 1.

###############################################
# Define conversion functions:

hom_map <- function(summary_dat)
{
 # convert summary file into map file

map <- as.data.frame(cbind(summary_dat[,1:2],rep(0,nrow(summary_dat)),summary_dat[,3]))
names(map) <- c("V1","V2","V3","V4")

return(map) # return map file
}

hom_ped <- function(hom_dat,indiv_dat,map){

roh <- hom_dat
per.roh <- indiv_dat  

# Create the ROH-by-SNP map

roh$start <- match(roh$POS1,map$V4)
roh$end <- match(roh$POS2,map$V4)

# set matrix up
roh.matrix <- matrix(1,nrow=nrow(per.roh),ncol=nrow(map)) # fill matrix with 1's, defualt genotype = 1 1
pers.id <- as.character(per.roh$IID)
snp.id <- as.character(map$V2)

# make the matrix
for (i in 1:nrow(per.roh)){
nmat <- roh[roh$IID==pers.id[i],]
inner.lp <- nrow(nmat)
if (inner.lp>0){
for (k in 1:inner.lp){roh.matrix[i,nmat$start[k]:nmat$end[k]] <-2} # fill regions with runs with 2's, thus run genotype = 2 2
 }
}

# Generate first 6 cols of .ped and combine with matrix

pre.ped <- as.data.frame(cbind(as.character(per.roh$FID),as.character(per.roh$IID),rep(0,nrow(per.roh)),rep(0,nrow(per.roh)),rep(0,nrow(per.roh)),rep(0,nrow(per.roh))))

# duplicate each column in matrix so as of form .ped file

ped.fill <- matrix(NA,nrow=nrow(roh.matrix),ncol=(2*ncol(roh.matrix)))

even.index <- seq(from=1,to=(2*ncol(roh.matrix)),by=2) 
odd.index <- (even.index+1)

ped.fill[,even.index] <- roh.matrix[,1:ncol(roh.matrix)]
ped.fill[,odd.index] <- roh.matrix[,1:ncol(roh.matrix)]

ped.fdata <- as.data.frame(cbind(pre.ped,ped.fill))

return(ped.fdata) # return ped file
}

######################################
# Run on All data sets:

datasets <- c('qc2report_scz_ab_eur_QC1B_0.02_0.02_0.02.TIME2.send5','qc2report_scz_bon_eur_QC1B_0.02_0.02_0.02.TIME2.send5','qc2report_scz_bulg_eur_QC1B_0.02_0.02_0.02.TIME2.send5','qc2report_scz_carwtc_eur_QC1B_0.02_0.02_0.02.TIME2.send5','qc2report_scz_cat2_eur_QC1B_0.02_0.02_0.02.TIME2.send5','qc2report_scz_dk_eur_QC1B_0.02_0.02_0.02.TIME2.send5','qc2report_scz_dub_eur_QC1B_0.02_0.02_0.02.TIME2.send5','qc2report_scz_edi_eur_QC1B_0.02_0.02_0.02.TIME2.send5','qc2report_scz_mgs2_eur_QC1B_0.02_0.02_0.02.TIME2.send5','qc2report_scz_muc_eur_QC1B_0.02_0.02_0.02.TIME2.send5','qc2report_scz_port_eur_QC1B_0.02_0.02_0.02.TIME2.send5','qc2report_scz_sw1_eur_QC1B_0.02_0.02_0.02.TIME2.send5','qc2report_scz_sw2_eur_QC1B_0.02_0.02_0.02.TIME2.send5','qc2report_scz_top3_eur_QC1B_0.02_0.02_0.02.TIME2.send5','qc2report_scz_ucl_eur_QC1B_0.02_0.02_0.02.TIME2.send5','qc2report_scz_ucla_eur_QC1B_0.02_0.02_0.02.TIME2.send5','qc2report_scz_zhh_eur_QC1B_0.02_0.02_0.02.TIME2.send5')

#dataset Number:
i <-17
  
# enter data set:
dataset <- datasets[i]
print("data set being read...")
# read data
hom_dat <- read.table(paste(dataset,".hom",sep=""),header=T)
indiv_dat <- read.table(paste(dataset,".hom.indiv",sep=""),header=T)
summary_dat <- read.table(paste(dataset,".hom.summary",sep=""),header=T)

# convert data and write out .map and .ped  and phenotype files and generate covariate file:

phe.dat <- as.data.frame(cbind(as.character(indiv_dat$FID),as.character(indiv_dat$IID),as.character(indiv_dat$PHE)))

names(phe.dat) <- c('FID','IID','PHE')

map_dat <- hom_map(summary_dat)

ped.dat <- hom_ped(hom_dat,indiv_dat,map_dat)

value <- substr(datasets[i],14,17)
covar <- cbind(ped.dat[,1:2], rep(value,nrow(ped.dat)))

# writing files:

write.table(phe.dat,file=paste(dataset,".phe",sep=""),quote=FALSE,sep=" ",row.names=FALSE,col.names=FALSE)

write.table(map_dat,file=paste(dataset,".map",sep=""),quote=FALSE,sep=" ",row.names=FALSE,col.names=FALSE)

write.table(ped.dat,file=paste(dataset,".ped",sep=""),quote=FALSE,sep=" ",row.names=FALSE,col.names=FALSE)

write.table(covar,file=paste(dataset,".cov",sep=""),quote=FALSE,sep="",row.names=FALSE,col.names=FALSE)





##########
# Now Merge the data sets into a single file:

                

