##############################################
#                                            #
#         Check clump make-up                #
#                                            #
##############################################

library(islandNeighbours)

#sort multiPhylo by size
#ntips<-sapply(tree,function(x) length(x$tip))
#sorted.trees<-tree[order(ntips,decreasing=F)]


#make list of files to open
temp<-list.files(pattern = "clump_*")

#get smallest trees in each cluster
smallTrees <- data.frame(matrix(ncol=3,nrow=0))
for (i in 1:length(temp)) {
  t <- as.multiPhylo(read.tree(temp[[i]]))
  #print(i)
  
  #equivalent to sapply(t,function(x) length(x$tip.label))
  ntips<-list() 
  for (j in 1:length(t)) {
    nleaves <- length(t[[j]]$tip)
    ntips <- c(ntips, as.vector(nleaves))
    #print(ntips)
  }
  
  ntipsVect <- as.vector(unlist(ntips))
  #print(ntipsVect)
  info <- c(temp[i], min(as.vector(ntipsVect)), length(which(ntipsVect == min(ntipsVect))))
  smallTrees <- rbind(smallTrees, info, make.row.names = F)
}
#rbind was overwritting colnames for empty df
colnames(smallTrees) <- c('clump', 'smallestTree', 'minTrees')

#save txt
write.table(smallTrees, file = 'smallTreesPerClump.txt', sep = '\t', row.names = F, col.names = T)
