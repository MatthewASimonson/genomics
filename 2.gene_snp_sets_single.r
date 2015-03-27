# Function 'gene_snp_sets()' Inputs:
# 1. name of dataset in quotes 
# 2. up/down stream additional base positions beyond gene to include as genic in KB
# 3. name of file that lists all unique autosomal genes in quotes; header format (NAME, START, STOP, CHR) ; 
# 4. name of file that lists all genes to include in .set file in quotes; format of file should be one gene name per row, gene names all capitalized, format OF NAME matching file glist-hg18 
# Function 'gene_snp_sets()' Outputs:
# 1. gene set file (.set extension)

gene_snp_sets <- function(dataset,KB,clean.list,mygenes.list){ # dataset and clean.list must be in quotes

# Read in list of cleaned autosomal genes
#########################################
genes <- read.table("KEGG.circadian.entrez",header=FALSE)
entrez <- read.table("entrezgenes.txt",header=FALSE) # read in generated gene list
pre <- entrez[,c(4,2,3,1)] # swap columns so in proper order to generate 'set' file
names(pre) <- 'V1'
set <- merge(genes,pre,by="V1")
rbx1 <- c('9978','41347351','41369313','22')
ts <- rbind(set,rbx1)
ts.f <- ts[,c(4,2,3,1)]
write.table(ts.f,file='temp.list',quote=FALSE,row.names=FALSE,col.names=FALSE)#
  
# Generate gene set SNP file with KB window surrounding gene list using PLINK:
################################################################################

system(paste('plink --bfile ',dataset,' --make-set temp.list --make-set-border ',KB,' --write-set --out ',dataset,'',sep=''))
system('rm temp.list') # remove gene file that has no header

} # end function


