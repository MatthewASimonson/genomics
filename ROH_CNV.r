# ROH CNV analysis:
# Matthew A. Simonson, October 2012
############################################
# create .cnv file from .hom.overlap file: #
############################################

#    FID     Family ID
#    IID     Individual ID
#    CHR     Chromosome
#    BP1     Start position (base-pair)
#    BP2     End position (base-pair)
#    TYPE    Type of variant, e.g. 0,1 or 3,4 copies # All runs should be listed as 1
#    SCORE   Confidence score associated with variant # Place zeros here
#    SITES   Number of probes in the variant # Place zeros here

# First extract all runs from .hom.overlap file:

datafile <- c('MERGE.clean.FINALMERGE_ROH_lite_snp65') # Be sure there is also a '.fam' file with this name
system(paste("grep '\\-9' ",datafile,".hom > ",datafile,".runs",sep="")) 

# Read in runs:

runs <- read.table(paste("",datafile,".runs",sep=""),header=FALSE)

# create output file:

output.cnv <- cbind.data.frame(runs$V2,runs$V3,runs$V5,runs$V8,runs$V9,rep(1,nrow(runs)),rep(0,nrow(runs)),rep(0,nrow(runs)))
order.index <- order(runs$V5,runs$V8)
output.cnv <- output.cnv[order.index,]

# write out .cnv file:
#install.packages("gdata")
library("gdata")
write.fwf(output.cnv,file=paste("",datafile,".cnv",sep=""),sep="\t",rownames=FALSE,colnames=FALSE,quote=FALSE)

# Use plink to create map file fore CNV/ROH data:

system(paste("plink --cnv-list ",datafile,".cnv --cnv-make-map --out ",datafile,"",sep=""))

# ROH mapping
system(paste("plink --cfile ",datafile," --pheno MERGE.clean.FINAL.phe --mperm 10000 --out ",datafile," &",sep="")) 

###################################################################
# Merge .mperm file with .bim file to generate .assoc format file #
###################################################################

roh.cnv.map <- read.table("test.cnv.roh.cnv.qt.summary.mperm",header=TRUE)
bim <- read.table("MERGE.clean.FINALMERGE_ROH_lite.bim",header=FALSE)
names(bim) <- c('CHR','SNP','CM','BP')


s <- read.table("two.side.test.cnv.qt.summary",header=TRUE)
sp <- read.table("two.side.test.cnv.qt.summary.mperm",header=TRUE)
t <- cbind(s,sp)

