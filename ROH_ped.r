# Create ROH ped file:

 a.fam <- read.table("WHI.NPU.commonROH_lite.fam",header=FALSE)
 a.hom <- read.table("WHI.NPU.commonROH_lite_snp65.hom",header=TRUE)
 a.bim <- read.table("WHI.NPU.commonROH_lite.bim",header=FALSE)

chroms <- c(1:22)
for(c in chroms){
# split by chromosome so no errors due to memory allocation:
print(c)
chr <- c
bim <- a.bim[which(a.bim$V1==chr),]
hom <- a.hom[which(a.hom$CHR==chr),]
hom$START <- match(hom$SNP1,bim$V2)
hom$STOP <- match(hom$SNP2,bim$V2)
fam <- a.fam
hom.matrix <- matrix(as.integer(1),nrow=nrow(fam),ncol=(nrow(bim))) # matrix of zeros, each row is subject, every column is SNP
null.matrix <- cbind(hom.matrix,hom.matrix)

evp <- function(char.vector){
  range <- eval(parse(text=char.vector)) # return ranges from each char index
  return(range)
}

for(i in 1:nrow(fam)){ # fill in hom.matrix with allele count at each loci inside a run
  if(nrow(hom[which(fam$V2[i]==hom$IID),])>0){ # skip subjects with no runs
  row.hom <- hom[which(fam$V2[i]==hom$IID),] # find rows in hom matrix for ind
  ind.char <- paste(row.hom$START,":",row.hom$STOP,sep="")
  fill.col.index <-  as.numeric(unlist(sapply(ind.char,evp)))
  hom.matrix[i,fill.col.index] <- 2 # fill the hom.matrix then bind with ids after again
  print(i)
}
}
print("creating .ped file...")
# fill ped file for chromosome
a1.index <- seq(from=1,to=ncol(hom.matrix)*2,by=2)
a2.index <- seq(from=2,to=ncol(hom.matrix)*2,by=2)
null.matrix[,a1.index] <- hom.matrix[,1:ncol(hom.matrix)]
null.matrix[,a2.index] <- hom.matrix[,1:ncol(hom.matrix)]
#


# Write out files:
write.table(bim[,1:4],file=paste("chr",chr,".roh.map",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
write.table(fam,file=paste("chr",chr,".roh.fam",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
write.table(null.matrix,file=paste("chr",chr,".roh.gen",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
#
system(paste("rm chr",chr,".roh.ped",sep=""))
system(paste("paste chr",chr,".roh.fam chr",chr,".roh.gen > chr",chr,".roh.ped",sep="")) # merge two files into .ped file
#
}
#
# merge files:

system("plink --file chr1.roh --merge-list roh.merge --make-bed --out GNHS.roh")


# examine merged data:

d <- read.table("WHI.NPU.test.qassoc",header=TRUE)
d$P[which(is.na(d$P)==TRUE)] <- 1 # p-value of 1 where no runs exist

# create manhattan plot:

source("http://people.virginia.edu/~sdt5z/0STABLE/qqman.r") # load manhattan functions
manhattan(d, colors=c("black","#666666","#CC6600"), pch=20, genomewideline=F, suggestiveline=F,main="WHI.NPU Autozygosity")

plink --bfile WHI.NPU.roh --assoc --pheno disc.resid.pheno --allow-no-sex --out WHI.NPU.test
