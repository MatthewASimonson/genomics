##### Convert Calls to tped file ####
## NOTE: This script produces a tped and tfam file from birdcall files produced by birdsuite and places it in the working directory. Read through the entire script before running it because some files that are read in are modified outside this script, such as certain ped and map files.


plates <- as.character(c('ANNUL','ARDOR','BENCH','BIDES','BLOOD','CHOMP','CORSE','CRAVE','EMEUS','FLIPS','GRAPY','IMAGE','INERT','INGLE','INNED','JOYED','KRAAL','MALTS','MAXES','MBIRA','MOTET','PEELS','PERDU','RIYAL','SAPID','SCHWA','SEELY','SHARD','STAYS','THYME','TOWNY','TREYS','TWAES','VULGO'))

range <- length(plates)

for (z in range){
setwd("/scratch/Bird.Calls/")
setwd(paste("/scratch/Bird.Calls/",plates[z], sep=""))
#### Step1: read in Map file:

file.map <- read.table("phg000013.unfiltered.SARC.plink.bim", header=FALSE) ## The original unfiltered bim file from dbgap must be in the working directory
                      
#### Step2: Read in birdseed calls and annotated summary files

calls.birdseed <- read.table("CALLS.larry_bird_calls", header=TRUE) ## larry_bird calls file
## rename columns of birdseed call files:
names(calls.birdseed)<- substr(as.character(names(calls.birdseed)),start=2,stop=6) # rename birdseed column names to subject ID's
calls.birdseed <- gsub(","," ",as.matrix(calls.birdseed))
# Find all call types from input file:
call.types <- unique(calls.birdseed[,ncol(calls.birdseed)])
# Replace CNV's with equivalent calls

calls.birdseed <- gsub("-1 -1","0 0",calls.birdseed)

## 1
calls.birdseed <- gsub("1 0","2 2",calls.birdseed)
calls.birdseed <- gsub("0 1","1 1",calls.birdseed)
calls.birdseed <- gsub("1 3","2 1",calls.birdseed)
calls.birdseed <- gsub("3 1","2 1",calls.birdseed)
calls.birdseed <- gsub("2 1","2 1",calls.birdseed)
calls.birdseed <- gsub("1 2","2 1",calls.birdseed)
calls.birdseed <- gsub("1 1","2 1",calls.birdseed)
calls.birdseed <- gsub("1 4","2 1",calls.birdseed)
calls.birdseed <- gsub("4 1","2 1",calls.birdseed)
calls.birdseed <- gsub("3 1","2 1",calls.birdseed)
calls.birdseed <- gsub("5 1","2 1",calls.birdseed)
calls.birdseed <- gsub("1 5","2 1",calls.birdseed)

## 2
calls.birdseed <- gsub("2 0","2 2",calls.birdseed)
calls.birdseed <- gsub("0 2","1 1",calls.birdseed)
calls.birdseed <- gsub("2 1","2 1",calls.birdseed)
calls.birdseed <- gsub("1 2","2 1",calls.birdseed)
calls.birdseed <- gsub("2 3","2 1",calls.birdseed)
calls.birdseed <- gsub("3 2","2 1",calls.birdseed)
calls.birdseed <- gsub("2 4","2 1",calls.birdseed)
calls.birdseed <- gsub("4 2","2 1",calls.birdseed)
calls.birdseed <- gsub("5 2","2 1",calls.birdseed)
calls.birdseed <- gsub("2 5","2 1",calls.birdseed)


## 3
calls.birdseed <- gsub("3 0","2 2",calls.birdseed)
calls.birdseed <- gsub("0 3","1 1",calls.birdseed)
calls.birdseed <- gsub("1 3","2 1",calls.birdseed)
calls.birdseed <- gsub("3 1","2 1",calls.birdseed)
calls.birdseed <- gsub("2 3","2 1",calls.birdseed)
calls.birdseed <- gsub("3 2","2 1",calls.birdseed)
calls.birdseed <- gsub("3 3","2 1",calls.birdseed)
calls.birdseed <- gsub("3 4","2 1",calls.birdseed)
calls.birdseed <- gsub("4 3","2 1",calls.birdseed)
calls.birdseed <- gsub("3 5","2 1",calls.birdseed)
calls.birdseed <- gsub("5 3","2 1",calls.birdseed)

## 4
calls.birdseed <- gsub("4 0","2 2",calls.birdseed)
calls.birdseed <- gsub("0 4","1 1",calls.birdseed)
calls.birdseed <- gsub("1 4","2 1",calls.birdseed)
calls.birdseed <- gsub("4 1","2 1",calls.birdseed)
calls.birdseed <- gsub("2 4","2 1",calls.birdseed)
calls.birdseed <- gsub("4 2","2 1",calls.birdseed)
calls.birdseed <- gsub("3 4","2 1",calls.birdseed)
calls.birdseed <- gsub("4 3","2 1",calls.birdseed)
calls.birdseed <- gsub("4 4","2 1",calls.birdseed)
calls.birdseed <- gsub("4 5","2 1",calls.birdseed)
calls.birdseed <- gsub("5 4","2 1",calls.birdseed)

##5
calls.birdseed <- gsub("6 0","2 2",calls.birdseed)
calls.birdseed <- gsub("5 0","2 2",calls.birdseed)
calls.birdseed <- gsub("0 6","1 1",calls.birdseed)
calls.birdseed <- gsub("1 5","2 1",calls.birdseed)
calls.birdseed <- gsub("5 1","2 1",calls.birdseed)
calls.birdseed <- gsub("2 5","2 1",calls.birdseed)
calls.birdseed <- gsub("5 2","2 1",calls.birdseed)
calls.birdseed <- gsub("3 5","2 1",calls.birdseed)
calls.birdseed <- gsub("5 3","2 1",calls.birdseed)
calls.birdseed <- gsub("5 4","2 1",calls.birdseed)
calls.birdseed <- gsub("4 5","2 1",calls.birdseed)
calls.birdseed <- gsub("5 5","2 1",calls.birdseed)



file.summary <- read.table("CALLS.annotated_summary", header=FALSE) ## birdseed summary file which contains format of map file at begining

head(file.summary[which(file.summary[,2]==1),1:3])
### Now remove copy number probe rows

row.snp.index <- grep("SNP_A",file.summary$V1)
snp.file.summary <- file.summary[row.snp.index,]

### remove every other row so only one copy of each SNP is listed

keep.index <- (1:(nrow(snp.file.summary)/2)*2)-1 ## generates an index of all odd numbers
a.snps <- snp.file.summary[keep.index,1:3] ## should double check to make sure no B's in col 4

## object chrom"i"index, is the index for the rows of chromsome i in a.snps; there are 24 chromosomes listed in Map file
for (i in 1:24){
eval(parse(text=paste("chrom",i,"index.map <- which(a.snps$V2==",i,")",sep='')))
}

### Create indeces for the row of each chromosome found in the summary file

## object chrom"i"index, is the index for the rows of chromsome i in a.snps; there are 24 chromosomes listed in summary file
for (i in 1:24){
eval(parse(text=paste("chrom",i,"index.summary <- which(a.snps$V2==",i,")",sep='')))
}

for (i in 1:24){
eval(parse(text=paste("match.index",i," <- match(substr(as.character(a.snps$V1[chrom",i,"index.summary]),start=1,stop=13),as.character(file.map$V2[chrom",i,"index.map]))",sep='')))

eval(parse(text=paste("nonmatches",i,"<- which(match.index",i,"==NA)",sep='')))
eval(parse(text=paste("which(nonmatches",i,">0)",sep=''))) ## Do any of the chromosome SNPs not match?
}
### All SNPs match

#### Step3: Generate tped file
snp.list <- substr(as.character(a.snps$V1),start=1,stop=13)
chrom<-a.snps$V2
zeros <- rep(0,length(chrom))
bp <- as.character(a.snps$V3)
map.format <-data.frame(cbind(chrom,snp.list,zeros,bp)) ## map file format data frame

## Sort calls data so ordered the same as map.format data rows
row.index <- match(as.character(map.format[,2]),as.character(calls.birdseed[,1]),)
sorted.rows.calls <- calls.birdseed[row.index,] ## checked at 500th, 5000th & 50000th row, and it does match

## Now add genotype data to map.format and write tped

tped <- data.frame(cbind(map.format,sorted.rows.calls[,2:ncol(sorted.rows.calls)])) ## Merge data into tped format for entire plate

write.table(tped,paste(plates[z],".tped",sep=""), quote=FALSE, row.names=FALSE,col.names=FALSE) ## write out tped file for plate

##### Create .tfam file ####

### Step1: read in first 5 rows of the original SARC and GRU ped files:
sorted.calls<- as.data.frame(sorted.rows.calls)

ped.data.SARC <- read.table("SARC.ped",header=FALSE)# NOTE: these are modified versions of the original ped files with only the first 6 columns
ped.data.GRU <- read.table("GRU.ped",header=FALSE)

## combine the SARC and GRU ped data:

ped.data <- rbind(ped.data.SARC, ped.data.GRU)

sorted.ped.data <- ped.data[match(ordered.subjects.id,ped.data$V1),] ## Now ped data is sorted in the same order as tped file and only includes individuals from this plate

write.table(sorted.ped.data,paste(plates[z],".tfam",sep=""), quote=FALSE, row.names=FALSE,col.names=FALSE) ## write out fam file for plate
} # end loop z
