############################################ 
# merge .mds files with Red.Imputed files: #
############################################

sets <- c('ab','bon','bulg','carwtc','cat','dk','dub','edi','mgs','muc','port','sw1','sw2','top3','ucl','ucla','zhh')

# Read in .fam files:

a10.fam <- read.table('Red.Imputed4a10.fam',header=FALSE) 
b10.fam <- read.table('Red.Imputed4b10.fam',header=FALSE)
c10.fam <- read.table('Red.Imputed4c10.fam',header=FALSE)
d10.fam <- read.table('Red.Imputed4d10.fam',header=FALSE)

fam.merge <- rbind(a10.fam,b10.fam,c10.fam,d10.fam)

# find index of each data set in .fam file:

ab.index.fam <- grep("ab",fam.merge$V1) #1
bon.index.fam <- grep("bon",fam.merge$V1) #2
bulg.index.fam <- grep("bulg",fam.merge$V1) #3
carwtc.index.fam <- grep("carwtc",fam.merge$V1) #4
cat.index.fam <- grep("cat",fam.merge$V1) #5
dk.index.fam <- grep("dk",fam.merge$V1) #6
dub.index.fam <- grep("dub",fam.merge$V1) #7
edi.index.fam <- grep("edi",fam.merge$V1) #8
mgs.index.fam <- grep("mgs",fam.merge$V1) #9
muc.index.fam <- grep("muc",fam.merge$V1) #10
port.index.fam <- grep("port",fam.merge$V1) #11
sw1.index.fam <- grep("sw1",fam.merge$V1) #12
sw2.index.fam <- grep("sw2",fam.merge$V1) #13
top3.index.fam <- grep("top3",fam.merge$V1) #14
ucl.index.fam <- grep("ucl",fam.merge$V1) #15
ucla.index.fam <- grep("ucla",fam.merge$V1) #16
zhh.index.fam <- grep("zhh",fam.merge$V1) #17

names(fam.merge) <- c('FID','IID')

# Read in .mds files:

ab.mds <- read.table('ab6.mds',header=TRUE) #1
bon.mds <- read.table('bon6.mds',header=TRUE) #2
bulg.mds <- read.table('bulg6.mds',header=TRUE) #3
carwtc.mds <- read.table('carwtc6.mds',header=TRUE) #4
cat2.mds <- read.table('cat26.mds',header=TRUE) #5
dk.mds <- read.table('dk6.mds',header=TRUE) #6
dub.mds <- read.table('dub6.mds',header=TRUE) #7
edi.mds <- read.table('edi6.mds',header=TRUE) #8
mgs2.mds <- read.table('mgs26.mds',header=TRUE) #9
muc.mds <- read.table('muc6.mds',header=TRUE) #10
port.mds <- read.table('port6.mds',header=TRUE) #11
sw1.mds <- read.table('sw16.mds',header=TRUE) #12
sw2.mds <- read.table('sw26.mds',header=TRUE) #13
top3.mds <- read.table('top36.mds',header=TRUE) #14
ucl.mds <- read.table('ucl6.mds',header=TRUE) #15
ucla.mds <- read.table('ucla6.mds',header=TRUE) #16
zhh.mds <- read.table('zhh6.mds',header=TRUE) #17

mds.merge <- rbind(ab.mds,bon.mds,bulg.mds,carwtc.mds,cat2.mds,dk.mds,dub.mds,edi.mds,mgs2.mds,muc.mds,port.mds,sw1.mds,sw2.mds,top3.mds,ucl.mds,ucla.mds,zhh.mds)

# find index of each data set in .mds file:

ab.index.mds <- grep("_ab_",keep.data[,1]) #1
bon.index.mds <- grep("_bon_",keep.data[,1]) #2
bulg.index.mds <- grep("_bulg_",keep.data[,1]) #3
carwtc.index.mds <- grep("_carwtc_",keep.data[,1]) #4
cat.index.mds <- grep("_cat2_",keep.data[,1]) #5
dk.index.mds <- grep("_dk_",keep.data[,1]) #6
dub.index.mds <- grep("_dub_",keep.data[,1]) #7
edi.index.mds <- grep("_edi_",keep.data[,1]) #8
mgs.index.mds <- grep("_mgs2_",keep.data[,1]) #9
muc.index.mds <- grep("_muc_",keep.data[,1]) #10
port.index.mds <- grep("_port_",keep.data[,1]) #11
sw1.index.mds <- grep("_sw1_",keep.data[,1]) #12
sw2.index.mds <- grep("_sw2_",keep.data[,1]) #13
top3.index.mds <- grep("_top3_",keep.data[,1]) #14
ucl.index.mds <- grep("_ucl_",keep.data[,1]) #15
ucla.index.mds <- grep("_ucla_",keep.data[,1]) #16
zhh.index.mds <- grep("_zhh_",keep.data[,1]) #17

sum(length(ab.index.mds)+length(bon.index.mds)+length(bulg.index.mds)+length(carwtc.index.mds)+length(cat.index.mds)+length(dk.index.mds)+length(dub.index.mds)+length(edi.index.mds)+length(mgs.index.mds)+length(muc.index.mds)+length(port.index.mds)+length(sw1.index.mds)+length(sw2.index.mds)+length(top3.index.mds)+length(ucl.index.mds)+length(ucla.index.mds)+length(zhh.index.mds))

indeces <- c(ab.index.mds,bon.index.mds,bulg.index.mds,carwtc.index.mds,cat.index.mds,dk.index.mds,dub.index.mds,edi.index.mds,mgs.index.mds,muc.index.mds,port.index.mds,sw1.index.mds,sw2.index.mds,top3.index.mds,ucl.index.mds,ucla.index.mds,zhh.index.mds)


##############
# .mds IID's have too much information, cut them so they match .fam files:
#

mds.data$IID <- gsub("_NOXLS","",mds.data$IID)

mds.data$IID <- gsub("1_case_scz_car_eur_A500K*","",mds.data$IID,fixed=TRUE)

mds.data$IID <- gsub("0_control_bip_wtc_eur_A500k*","",mds.data$IID,fixed=TRUE)

##############
# CHECK TO MAKE SURE FID AND IID ARE THE SAME IN .mds and .fam files:
#############
# Create new ID column for .mds and .fam that is merge of FID and IID:

mds.names <- paste(mds.merge$FID,mds.merge$IID,sep="")
fam.names <- paste(fam.merge$FID,fam.merge$IID,sep="")

mds.total <- cbind(mds.names,mds.merge) # merge with .mds
fam.total <- cbind(fam.names,fam.merge) # merge with .fam

# now make merged name column have same name in both mds.total and fam.total:

names(mds.total) <- c('NAME',names(mds.total[2:ncol(mds.total)]))
names(fam.total) <- c('NAME',names(fam.total[2:ncol(fam.total)]))

# Check to be sure .mds data and .fam data can now merge:

match.index <- match(as.character(fam.total$NAME),as.character(mds.total$NAME))

nomatch.index <- which(is.na(match.index)) # this should be zero

###########
# Now order rows of mds.merge so same as fam.merge
###########

mds.merge <- mds.merge[match.index,]

###########
# Write out seperate .mds files for each of the 4 .fam files
###########

a7.range <- 1:nrow(a7.fam)
b7.range <- (nrow(a7.fam)+1):(nrow(a7.fam)+nrow(b7.fam))
c7.range <- (nrow(b7.fam)+1):(nrow(b7.fam)+nrow(c7.fam))
d7.range <- (nrow(c7.fam)+1):(nrow(c7.fam)+nrow(d7.fam))

a7.mds <- mds.merge[a7.range,]
b7.mds <- mds.merge[b7.range,]
c7.mds <- mds.merge[c7.range,]
d7.mds <- mds.merge[d7.range,]

write.table(a7.mds,file="Red.Imputed4a7.mds",quote=FALSE,row.names=FALSE,col.names=TRUE)#
write.table(b7.mds,file="Red.Imputed4b7.mds",quote=FALSE,row.names=FALSE,col.names=TRUE)#
write.table(c7.mds,file="Red.Imputed4c7.mds",quote=FALSE,row.names=FALSE,col.names=TRUE)#
write.table(d7.mds,file="Red.Imputed4d7.mds",quote=FALSE,row.names=FALSE,col.names=TRUE)#

/home/simonsom/PGC.ROH/map/Final_ROH .
