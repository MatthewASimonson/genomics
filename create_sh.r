msig <- as.matrix(read.table("misgdb.ranges",header=FALSE))

for(i in 1:nrow(msig)){
system(paste("cp x.make.sets ",strsplit(msig[i,1],'\\.')[[1]][1],".sh",sep=""))
system(paste("sed -i 's/x/",strsplit(msig[i,1],'\\.')[[1]][1],"/g' ",strsplit(msig[i,1],'\\.')[[1]][1],".sh",sep=""))
print(i)
}

for(i in 1:nrow(msig)){
system(paste("qsub ",strsplit(msig[i,1],'\\.')[[1]][1],".sh -A UCB00000149",sep=""))
print(i)
}
