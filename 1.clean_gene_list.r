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

} # end function


