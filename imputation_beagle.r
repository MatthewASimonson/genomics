# Impute data from all sets before 
# NOTE: include some hapmap offspring

# STEP 1:
##########################################
# Make sure all assemblies are the same: #
##########################################

# beagle reference genome is based on NCBI build 37 (UCSC hg19)

# Read in SNPs from all data sets:

setwd("/home/simonsom/ROH_pathway/TOTAL/GENEVA_diab/TOTAL")
g.bim <- read.table("GENEVA.total.final.bim",header=FALSE)
setwd("/home/simonsom/ROH_pathway/TOTAL/WHI")
w.bim <- read.table("WHI.clean.final.bim",header=FALSE)
setwd("/home/simonsom/ROH_pathway/TOTAL")
amc.bim <- read.table("MERGE.clean.FINAL.bim",header=FALSE)
setwd("/home/simonsom/ROH_pathway/TOTAL/MGS/mgs_plink_raw")
mgs.bim <- read.table("MGS.final.bim",header=FALSE)

ts.1 <- rbind(g.bim,w.bim)
ts.2 <- rbind(amc.bim,mgs.bim)
ts.pre <- rbind(ts.1,ts.2)
ts.pre <- ts.pre[order(ts.pre$V1,ts.pre$V4),]
total.snps <- unique(ts.pre$V2)
setwd("/home/simonsom/ROH_pathway/TOTAL")
write.table(total.snps,file="total.old.snps",row.names=FALSE,col.names=FALSE,quote=FALSE) # write out all snps from across sets

# Use table browser to generate file with new positions and correct strand and remove snps that differ between AFFY verions:

setwd("/STATGEN/home/simonsom/beagle/strand")

g <- read.table("GenomeWideSNP_6.na32.strand.txt",skip=1,header=FALSE) # affy snps strand info
setwd("/home/simonsom/ROH_pathway/TOTAL")
new.snps <- read.table("total.new.snps",header=FALSE,skip=1)

ambi.snps <- as.character(g$V2[duplicated(g$V2)]) # ambigous snps to drop
ambi.snps <- c(ambi.snps,'rs12179084')

names(new.snps) <- c('C','BP1','BP2','SNP')
names(g) <- c('AFFY','SNP','CHR','BP','STRAND')

tot.snps <- merge(new.snps,g,by="SNP")
strand.flip <- tot.snps$SNP[which(tot.snps$STRAND=='-')]
update.pos <- cbind.data.frame(tot.snps$SNP,tot.snps$BP2)
keep.snps <- update.pos[,1]

setwd("/home/simonsom/ROH_pathway/TOTAL/GENEVA_diab/TOTAL")
write.table(strand.flip,file="strand.flip",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(ambi.snps,file="ambi.snps",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(update.pos,file="update.pos",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(keep.snps,file="keep.snps",row.names=FALSE,col.names=FALSE,quote=FALSE)
system("plink --bfile GENEVA.total.final --exclude ambi.snps --make-bed --out GENEVA.preflip")
system("plink --bfile GENEVA.preflip --flip strand.flip --update-map update.pos --extract keep.snps --make-bed --out GENEVA.total.final")

setwd("/home/simonsom/ROH_pathway/TOTAL")
write.table(strand.flip,file="strand.flip",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(ambi.snps,file="ambi.snps",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(update.pos,file="update.pos",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(keep.snps,file="keep.snps",row.names=FALSE,col.names=FALSE,quote=FALSE)
system("plink --bfile MERGE.clean.FINAL --exclude ambi.snps --make-bed --out MERGE.preflip")
system("plink --bfile MERGE.preflip --flip strand.flip --update-map update.pos --extract keep.snps --make-bed --out MERGE.clean.FINAL")

setwd("/home/simonsom/ROH_pathway/TOTAL/MGS/mgs_plink_raw")
write.table(strand.flip,file="strand.flip",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(ambi.snps,file="ambi.snps",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(update.pos,file="update.pos",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(keep.snps,file="keep.snps",row.names=FALSE,col.names=FALSE,quote=FALSE)
system("plink --bfile MGS.final --exclude ambi.snps --make-bed --out MGS.preflip")
system("plink --bfile MGS.preflip --flip strand.flip --update-map update.pos --extract keep.snps --make-bed --out MGS.FINAL")

# Now for Illumina data: (WHI)

setwd("/home/simonsom/ROH_pathway/TOTAL/WHI")
