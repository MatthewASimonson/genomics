############################################################
# FOLLOWING WAS WRITTEN BY MATTHEW A. SIMONSON 2011-2012   #
############################################################


# Function 'clean_gene_list' Inputs:
# 1. name of gene list file; columns (CHR,START,STOP,NAME)
# 2. name of output file
#
# Function 'clean_gene_list' Outputs:
# 1. list of unique autosomal genes with header (NAME, CHR, START, STOP)

clean_gene_list <- function(list.file,out.filename){ # inputs must be character (in quotes)
# read in list of genes
#######################
g.list <- read.table(list.file,header=FALSE) #

# re-arrange columns
####################
g.list <- cbind(g.list[,4],g.list[,1:3])

# give proper column names:
###########################
names(g.list) <- c("NAME","CHR","START","STOP")

# Order all genes by location:
##############################

order.index <- order(as.numeric(as.character(g.list$CHR)),as.numeric(as.character(g.list$START)))
g.list.ordered <- g.list[order.index,]

# Remove possible duplicates:
#############################

g.list.un <- g.list.ordered[as.logical((duplicated(g.list.ordered$NAME)*-1)+1),] # weird syntax that should be accomplished with 'unique' but wasn't working for some reason???

# remove all genes on X and Y chromosomes and NA's:
###################################################

no.X.index <- which(g.list.un$CHR!='X')
no.X.g.list <- g.list.un[no.X.index,]

no.Y.index <- which(no.X.g.list$CHR!='Y')
no.XY.g.list <- no.X.g.list[no.Y.index,]

autosome.index <- which(is.na(no.XY.g.list$CHR)==FALSE)
autosome.list <- no.XY.g.list[autosome.index,]

# write out final list of genes without duplicates:
#####################################################

final.list <- autosome.list

write.table(final.list,file=out.filename,quote=FALSE,row.names=FALSE,col.names=TRUE)#

} # END FUNCTION clean_gene_list

# Function 'gene_snp_sets()' Inputs:
# 1. name of dataset in quotes 
# 2. up/down stream additional base positions beyond gene to include as genic in KB
# 3. name of file that lists all unique autosomal genes in quotes; header format (NAME, START, STOP, CHR) ; 
# 4. name of file that lists all genes to include in .set file in quotes; format of file should be one gene name per row, format OF NAME matching file glist-hg18 
# Function 'gene_snp_sets()' Outputs:
# 1. gene set file (.set extension)

gene_snp_sets <- function(dataset,KB,clean.list,mygenes.list){ # dataset and clean.list must be in quotes

# Read in list of cleaned autosomal genes
#########################################
mygenes.list <- read.table(mygenes.list,header=FALSE)
gene.list <- read.table(clean.list,header=TRUE) # read in generated gene list
gene.count <- nrow(gene.list)

# generate subdirectory to store all binary plink files
#######################################################
curdir <- getwd()

system(paste("mkdir ",dataset,"_plink_bin",sep=""))# generate subdirectory to store all binary plink files
system(paste("cp ",dataset,".bim ",curdir,"/",dataset,"_plink_bin",sep="")) # copy .bim file to subdirectory
system(paste("cp ",dataset,".bed ",curdir,"/",dataset,"_plink_bin",sep="")) # copy .bed file to subdirectory
system(paste("cp ",dataset,".fam ",curdir,"/",dataset,"_plink_bin",sep="")) # copy .fam file to subdirectory

setwd(paste("",curdir,"/",dataset,"_plink_bin",sep="")) # change working directory to subdirectory so individual gene files are created in proper location

rearrange.list <- gene.list[,c(2:4,1)] # swap columns so in proper order to generate 'set' file
write.table(rearrange.list,file='temp.list',quote=FALSE,row.names=FALSE,col.names=FALSE)#
  
# Generate gene set SNP file with KB window surrounding gene list using PLINK:
################################################################################

system(paste('plink --bfile ',dataset,' --make-set temp.list --make-set-border ',KB,' --write-set --out ',dataset,'',sep=''))
system('rm temp.list') # remove gene file that has no header

} # end function

###############################
# CALL AND EXECUTE FUNCTIONS: #
###############################


# INPUT EXPECTED ON COMMAND LINE:
# R CMD BATCH --no-save --no-restore '--args list.file="glist-hg18" cleanlist.filename="clean.genelist" dataset="ARIC.clean" KB=20 mygenes.list="example.list.txt"' generate_clean_set.R &
args=(commandArgs(TRUE)) # read in command line arguments

if(length(args)==0){
    print("No arguments supplied.")
}else{
    for(i in 1:length(args)){
         eval(parse(text=args[[i]]))
    }
}

clean_gene_list(list.file,cleanlist.filename)
gene_snp_sets(dataset,KB,cleanlist.filename,mygenes.list) # NOTE: this function should read in the 'out.filename' file from first function               

