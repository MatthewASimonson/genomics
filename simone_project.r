#By: Matthew A. Simonson, 10/25/12@4:54pm


map <- read.table("2012_09_19_mapping_file_matt.txt",header=TRUE,sep="\t",fill=NA)
dist <- read.table("2012_09_24_unweighted_unifrac_distance_matrix_for_matt.txt",header=TRUE,sep="\t",fill=NA)
total.index <- match(dist[,1],map[,1])
map <- map[total.index,]# only same subjects as in distance
otu <- read.table("otu_table_with_lineage.txt",header=TRUE,sep="\t",fill=NA)

# create data frame:

pairwise.dat <- matrix(NA,nrow=nrow(dist)^2,ncol=43) # this matrix has 43 columns; col1 is ID1, col2 is ID2, col 3 is distance, col 4 through 44 is phenotype difference or categorical level (factor level)

categorical.index <- c(3:13,16:18,21) # specify columns with categorical data here
numeric.index <- c(14:15,19:20,22:42) # specify columns with numeric data here

# fill top half and diagonal of dist with NA's to remove redundancy in next step
dist.fill <- dist[1:nrow(dist),2:ncol(dist)]
for(i in 1:nrow(dist.fill)){
  for(j in i:ncol(dist.fill)){
    dist.fill[i,j] <- NA
  }
   print(i) # counts to 226
}

dist[1:nrow(dist),2:ncol(dist)] <- dist.fill

# fill data matrix:

for(i in 1:nrow(dist)){
  for(j in 1:(ncol(dist)-1)){
    from.row <- i # take data from this row in distance
    from.col <- j+1 # take data from this col in distance
    to.index <- j+((i-1)*nrow(dist)) # which row in pairwise.dat to place data in
    pairwise.dat[to.index,3] <- dist[from.row,from.col]
    pairwise.dat[to.index,1] <- as.character(dist$SampleID[i])
    pairwise.dat[to.index,2] <- as.character(dist$SampleID[j])

    # categorical columns:
    sorted.cat <- apply(rbind(as.matrix(map[i,categorical.index]),as.matrix(map[j,categorical.index])),2,sort) # alphabetical order merging to prevent non-unique combinations of factors
    pairwise.dat[to.index,categorical.index+1] <- apply(sorted.cat,2,paste,collapse=".") # now paste together 
 
    # numeric columns:
    pairwise.dat[to.index,numeric.index+1] <- abs(as.numeric(as.character(map[i,numeric.index]))-as.numeric(as.character(map[j,numeric.index])))
    
      }# end loop j
  print(i) # counts to 226
    } # end loop i
#
pairwise.frame <- as.data.frame(pairwise.dat) # convert to data frame from matrix
names(pairwise.frame) <- c('ID1','ID2','DIST',names(map)[5:(ncol(map)-1)]) # give columns appropriate names
pairwise.final <- pairwise.frame[is.na(pairwise.frame$DIST)==FALSE,] # remove rows comparing individual to themself and repeat combinations
#







# non zero index:

x <- -5:5
x.noz <- x[x!=0]

# apply wilcox.test to each row:
wilc.result <- as.data.frame(matrix(NA,nrow=nrow(otu),ncol=2)) # store 
names(wilc.result) <- c('W','p-val')

for(i in 1:nrow(otu)){
wilc.result[i,1] <- wilcox.test(as.numeric(otu[i,]),as.numeric(otu[i,]))$statistic
wilc.result[i,2] <- wilcox.test(as.numeric(otu[i,]),as.numeric(otu[i,]))$p.value
}
