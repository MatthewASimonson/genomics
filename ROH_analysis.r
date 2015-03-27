
#############################
# Run this on imputed data: #
#############################


# STEP 1:
#############################################
# Run ROH Analysis on Data Sets Seperately: #
#############################################

#Set data set:
data_set <- "GLAUC.clean"
data_set <- "ADDICT.clean"
data_set <- "PREMA.clean"
data_set <- "LUNG.clean"
data_set <- "VTHROM.clean"

#Remove MAF < 0.05
system(paste("plink --bfile ",data_set," --maf .05 --make-bed --out ",data_set,"ROH_maf",sep=""))

#Light LD-pruning
system(paste("plink --bfile ",data_set,"ROH_maf --indep 50 5 10 --out ",data_set,"ROH_lite",sep="")) # set to run in background
system(paste("plink --bfile ",data_set,"ROH_maf --extract ",data_set,"ROH_lite.prune.in --make-bed --out ",data_set,"ROH_lite",sep=""))
system(paste("cp ",data_set,"ROH_lite.log ",data_set,"ROH_lite.log",sep=""))

# --- Genotype Info in pruned data

#light LD-pruning
system(paste("plink --bfile ",data_set,"ROH_lite --missing --out ",data_set,"ROH_lite",sep="")) #missingness
system(paste("plink --bfile ",data_set,"ROH_lite --het --out ",data_set,"ROH_lite",sep="")) #Inb. coeff
system(paste("plink --bfile ",data_set,"ROH_lite --freq --out ",data_set,"ROH_lite",sep=""))
frq <- read.table(paste("",data_set,"ROH_lite.frq",sep=""),header=TRUE)
list <- frq[,c(2,3)]
one <- rep(1,nrow(frq))
list <- cbind(list,one)
write.table(list,file=paste("",data_set,".list",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
system(paste("plink --bfile ",data_set,"ROH_lite --score ",data_set,".list --out ",data_set,"ROH_lite",sep="")) #minor allele load

# 4) Lite LD pruning - 65 SNP threshold - No het allowance (ROH_lite_snp65)

# Parameters:
# --homozyg-window-snp 65
# --homozyg-snp 65
# --homozyg-window-missing 3
# --homozyg-window-het 0
# --homozyg-window-threshold 0.03

# --homozyg-window-kb 0
# --homozyg-kb 0
# --homozyg-gap 2000
# --homozyg-density 300

# --homozyg-group
# --homozyg-match .95

# 5) Run ROH burden analysis on full genome:
setwd("/STATGEN/home/simonsom/ROH_pathway/imputed/GNHS.roh")
msig <- as.matrix(read.table("misgdb.ranges",header=FALSE))


# 6) Calculate total size of each pathway so ROH % can be determined:
path.size <- vector()  
  for(i in 1:nrow(msig)){
temp <- read.table(msig[i,1],header=FALSE)
path.size[i] <- sum(temp$V3-temp$V2)
print(i)
  }

# 7) Find runs that overlap with pathways:

setwd("/STATGEN/home/simonsom/ROH_pathway/imputed/GNHS.roh")
g.hom <- read.table("GNHS.commonROH_lite_snp65.hom",header=TRUE)
setwd("/STATGEN/home/simonsom/ROH_pathway/imputed/ARIC.roh")
a.hom <- read.table("ARIC.commonROH_lite_snp65.hom",header=TRUE)
setwd("/STATGEN/home/simonsom/ROH_pathway/imputed/WHI.GRU.roh")
wg.hom <- read.table("WHI.GRU.commonROH_lite_snp65.hom",header=TRUE)
setwd("/STATGEN/home/simonsom/ROH_pathway/imputed/WHI.NPU.roh")
wn.hom <- read.table("WHI.NPU.commonROH_lite_snp65.hom",header=TRUE)
#
  
rohs <- a.hom
path.runs <- list()  
  for(i in 1:nrow(msig)){
  temp <- read.table(msig[i,1],header=FALSE)
    for(j in 1:nrow(temp)){
   roh.list <- list()   
   in.index <- which((rohs$CHR==temp$V1[j]) & (rohs$POS1>=temp$V2[j]) & (rohs$POS2<=temp$V3[j])) # runs inside gene
   front.index <- which((rohs$CHR==temp$V1[j]) & (rohs$POS1<=temp$V2[j]) & (rohs$POS2>=temp$V2[j])) # runs front overlap gene
   back.index <- which((rohs$CHR==temp$V1[j]) & (rohs$POS1<=temp$V3[j]) & (rohs$POS2>=temp$V3[j])) # runs back overlap gene
   enc.index <- which((rohs$CHR==temp$V1[j]) & (rohs$POS1<=temp$V2[j]) & (rohs$POS2>=temp$V3[j])) # run that encompasses gene
   roh.index <- c(front.index,back.index,enc.index)
   roh.temp <- rohs[roh.index,]
   roh.list[[j]] <- roh.temp
   }
 path.runs[[i]] <- as.data.frame(do.call(rbind, roh.list))
  print(i)
  }
  #
path.len <- vector()
  for(i in 1:880){
    path.len[i] <- nrow(path.runs[[i]])
    print(i)
                        }

  
merge.paths <- list()
    for(i in 1:880){
merge.paths[[i]] <- cbind.data.frame(path.runs[[i]],rep(strsplit(msig[i],"\\.")[[1]][1],path.len[[i]]),rep(path.size[[i]],path.len[[i]]))
names(merge.paths[[i]]) <- c(names(rohs),'PATH','BP')
merge.paths[[i]]$RATE <- (merge.paths[[i]]$KB*1000)/merge.paths[[i]]$BP 
print(i)
    }

id.rate <- list()  
for(i in 1:880){ # list filled with subjects ID and ROH rate in each pathway
id.rate[[i]] <- cbind.data.frame(merge.paths[[i]]$IID,merge.paths[[i]]$RATE)
names(id.rate[[i]]) <- c('IID',strsplit(msig[i],"\\.")[[1]][1])
print(i)
}
 save(id.rate,file="ARIC.roh.Rdata")
#
 disc.phe <- read.table("disc.resid.pheno",header=FALSE)
 names(disc.phe) <- c('FID','IID','PHE')
setwd("/STATGEN/home/simonsom/ROH_pathway/imputed/ARIC.roh")
load("ARIC.roh.Rdata")
aric.id.rate <- id.rate
rm(id.rate)
setwd("/STATGEN/home/simonsom/ROH_pathway/imputed/GNHS.roh")
load("GNHS.roh.Rdata")
gnhs.id.rate <- id.rate
rm(id.rate)
setwd("/STATGEN/home/simonsom/ROH_pathway/imputed/WHI.GRU.roh")
load("WHI.GRU.roh.Rdata")
whi.gru.id.rate <- id.rate
rm(id.rate)
setwd("/STATGEN/home/simonsom/ROH_pathway/imputed")
load("WHI.NPU.roh.Rdata")
whi.npu.id.rate <- id.rate
rm(id.rate)

path.model.matrix <- disc.phe[duplicated(disc.phe$IID)==FALSE,] # merge all pathways with phenotype data for multiple regression
for(i in 1:880){
path.model.matrix <- merge(path.model.matrix,aric.id.rate[[i]],by="IID",all=TRUE)
path.model.matrix <- path.model.matrix[duplicated(path.model.matrix$IID)==FALSE,]
print(i)
}
save(path.model.matrix,file="ARIC.roh.matrix.Rdata")
aric.matrix <- path.model.matrix
#
for(i in 1:880){
path.model.matrix <- merge(path.model.matrix,gnhs.id.rate[[i]],by="IID",all=TRUE)
path.model.matrix <- path.model.matrix[duplicated(path.model.matrix$IID)==FALSE,]
print(i)
}
save(path.model.matrix,file="GNHS.roh.matrix.Rdata")
gnhs.matrix <- path.model.matrix
#
for(i in 1:880){
path.model.matrix <- merge(path.model.matrix,whi.gru.id.rate[[i]],by="IID",all=TRUE)
path.model.matrix <- path.model.matrix[duplicated(path.model.matrix$IID)==FALSE,]
print(i)
}
save(path.model.matrix,file="WHI.GUR.roh.matrix.Rdata")
whi.gru.matrix <- path.model.matrix
#
for(i in 1:880){
path.model.matrix <- merge(path.model.matrix,whi.npu.id.rate[[i]],by="IID",all=TRUE)
path.model.matrix <- path.model.matrix[duplicated(path.model.matrix$IID)==FALSE,]
print(i)
}
save(path.model.matrix,file="WHI.NPU.roh.matrix.Rdata")
whi.npu.matrix <- path.model.matrix
#

pmp1 <- rbind(aric.matrix,gnhs.matrix)
pmp2 <- rbind(whi.gru.matrix,whi.npu.matrix)
path.model.matrix <- rbind(pmp1,pmp2)

na.counter <- function(x){
  return(length(which(is.na(x)==TRUE)))}

row.nas <- apply(path.model.matrix,1,na.counter)
keep.rows <- c(1:9155,35152:37702,45977:50290,64365:68679)
path.model.matrix <- path.model.matrix[keep.rows,]
save(path.model.matrix,file="disc.model.matrix.Rdata")

path.model.matrix <- as.matrix(path.model.matrix)
  for(i in 1:(ncol(path.model.matrix)-3)){
 path.model.matrix[is.na(path.model.matrix[,(3+i)]),(3+i)] <- 0 # fill NA's with zero
 print(i)
 }
#
path.roh.count <- vector() # number of runs in each pathway
for(i in 1:(ncol(path.model.matrix)-3)){
path.roh.count[i] <- sum(path.model.matrix[,i+3])
                         print(i)
}

path.keep <- c(1:3,(which(path.roh.count>1)+3))
path.model.matrix <- as.data.frame(path.model.matrix)
path.model.matrix <- path.model.matrix[,path.keep]
path.model.matrix <- path.model.matrix[,-c(which(names(path.model.matrix)=='REACTOME_VITAMIN_B5_(PANTOTHENATE)_METABOLISM'),which(names(path.model.matrix)=='REACTOME_CDC6_ASSOCIATION_WITH_THE_ORC:ORIGIN_COMPLEX'))]
kegg.paths <- path.model.matrix[214:387]


# Just examine KEGG pathways
betas <- vector()
ses <- vector()
tvals <- vector()
pvals <- vector()

for(i in 1:length(names(path.model.matrix)[214:387])){
  ROH.variable <- names(path.model.matrix)[214:387][i] # 
  model.paths <- eval(parse(text=(paste("lm(as.numeric(PHE) ~ ",ROH.variable,", data=path.model.matrix)",sep="")))) 
   betas[i] <- summary(model.paths)$coefficients[2,1]
   ses[i] <- summary(model.paths)$coefficients[2,2]
   tvals[i] <- summary(model.paths)$coefficients[2,3]
   pvals[i] <- summary(model.paths)$coefficients[2,4]
  print(i)
}

mod.matrix <- cbind.data.frame(names(path.model.matrix)[214:387],betas,ses,tvals,pvals)    

# permutation:

perm.dat <- matrix(0,nrow=length(tvals),ncol=1000)
perm.model.matrix <- path.model.matrix
PHEs <- as.numeric(perm.model.matrix[,3])
for(j in 1:1000){
for(i in 1:length(names(path.model.matrix)[214:387])){
  PHEs <- sample(PHEs,replace=FALSE)
  ROH.variable <- names(path.model.matrix)[214:387][i] # 
  model.paths <- eval(parse(text=(paste("lm(PHEs ~ ",ROH.variable,", data=perm.model.matrix)",sep="")))) 
  perm.dat[i,j] <- summary(model.paths)$coefficients[2,3]
  print(i)
  print(j)
}
}

perm.path.srt <- apply(perm.dat,1,sort)

lower.tail <- perm.path.srt[1,]
upper.tail <- perm.path.srt[nrow(perm.path.srt),]

adj.p <- vector()
for(i in 1:nrow(perm.dat)){
if((tvals[i]<mean(lower.tail))&(tvals[i]<0)){  
pval <- pnorm(-abs((tvals[i]-mean(lower.tail))/sd(lower.tail)))
}
if((tvals[i]>mean(lower.tail))&(tvals[i]<0)){
pval <- pnorm(-abs((tvals[i]-mean(lower.tail))/(sd(lower.tail))),lower.tail=FALSE)
}
if((tvals[i]>mean(upper.tail))&(tvals[i]>0)){  
pval <- pnorm(-abs((tvals[i]-mean(upper.tail))/(sd(upper.tail))))
}
if((tvals[i]<mean(upper.tail))&(tvals[i]>0)){
pval <- pnorm(-abs((tvals[i]-mean(upper.tail))/(sd(upper.tail)/(sd(upper.tail)))),lower.tail=FALSE)
}   
adj.p[i] <- pval
print(i)
}

mod.matrix$adj.p <- adj.p
mm.order <- mod.matrix[order(mod.matrix$pvals),]
save(mm.order,file="discover.roh.results.Rdata")


##############################
# Examine genes in pathways:
##############################
g <- read.table("genes.files",header=FALSE)
sig.path.genes <- list()
for(i in 1:nrow(g)){
  sig.path.genes[[i]] <- read.table(as.character(g[i,1]),header=FALSE)
  print(i)
}
u.genes <- unique(unlist(sig.path.genes)

                  
count.matrix <- matrix(0,nrow=nrow(g),ncol=nrow(g))
pct.matrix <- matrix(0,nrow=nrow(g),ncol=nrow(g)) 
overlap.list <- array(NA,dim=c(nrow(g),nrow(g),200))
for(i in 1:nrow(g)){
for(j in 1:nrow(g)){
  count.matrix[i,j] <- length(intersect(sig.path.genes[[i]]$V1,sig.path.genes[[j]]$V1))
  if(length(intersect(sig.path.genes[[i]]$V1,sig.path.genes[[j]]$V1))>0){
  overlap.list[i,j,1:length(intersect(sig.path.genes[[i]]$V1,sig.path.genes[[j]]$V1))] <- intersect(sig.path.genes[[i]]$V1,sig.path.genes[[j]]$V1)
}
}
}

o.index<- which(count.matrix>0,arr.ind=TRUE)
o.index.ns <- o.index[which(o.index[,1]!=o.index[,2]),]
                  
