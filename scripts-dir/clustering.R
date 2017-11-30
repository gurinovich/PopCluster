args = commandArgs(trailingOnly=TRUE)

work.dir <- args[1]
system(paste0("mkdir ", work.dir,"results"),wait=T)
system(paste0("mkdir ", work.dir,"results/clusters"),wait=T)
clusters.dir <- paste0(work.dir, "results/clusters/")

##### FUNCTIONS START ##########

#returns elements (repeating elements are considered separate elements) that are in x but not i y
vectdiff <- function(x,y) {
  x <- as.vector(x)
  y <- as.vector(y)
  z <- vector()
  posx <- 1
  for (i in x) {
    posy <- 1
    for (j in y) {
      if (i == j) {
        z <- c(z,posx)
        y <- y[-posy]
        break    
      }
      posy <- posy + 1
    }
    posx <- posx + 1
  }
  if (length(z)==0) {
    x
  } else {x[-z]}
}

##### FUNCTIONS END ##########

#load data
mega.data <- read.csv(paste0(work.dir, "mega-data.csv"))

numberPCs <- 6
pcs <- vector()
for (i in 1:numberPCs) {
	pcs <- c(pcs,paste0("PC",i))
}
max.clust.size <- 100

#save original parent cluster:
write.table(mega.data[,c("Family","ID")],file=paste0(clusters.dir,"cluster.",nrow(mega.data),".txt"), quote=F, row.names=F, col.names=F)

mega.data1 <- mega.data[,pcs]
rownames(mega.data1) <- mega.data$ID
mega.data1 <- scale(mega.data1)
d <- dist(mega.data1)
fit <- hclust(d)

curr.clusters <- c(nrow(mega.data))
all.clusters <- curr.clusters
name.add <- 1
ncl <- 2

while(!(length(which(curr.clusters>max.clust.size))==0))  {
  new.cut <- as.vector(table(cutree(fit,ncl)))
  ncl <- ncl + 1  
  new.clusters <- vectdiff(new.cut,curr.clusters)
  prev.cut <- curr.clusters
  curr.clusters <- new.cut  
  if (length(which(new.clusters>max.clust.size))==0) next
  groups <- cutree(fit,k=length(new.cut))
  groups.df <- as.data.frame(groups)
  groups.df$ID <- rownames(groups.df)
  mega.data.clust <- merge(groups.df,mega.data,by="ID")
  
  
  if (new.clusters[1]==new.clusters[2]) {
    pos <- which(new.cut == new.clusters[1])
    write.table(mega.data.clust[,c("Family","ID")][which(mega.data.clust$groups==pos[1]),],file=paste0(clusters.dir,"cluster.",new.clusters[1],".",name.add,".txt"),quote=F, row.names=F,col.names=F) 
    name.add <- name.add + 1
    write.table(mega.data.clust[,c("Family","ID")][which(mega.data.clust$groups==pos[2]),],file=paste0(clusters.dir,"cluster.",new.clusters[1],".",name.add,".txt"), quote=F, row.names=F,col.names=F)
    name.add <- name.add + 1
  } else {
  
  
  for (i in new.clusters) {
    if (i > max.clust.size) {
      pos <- which(new.cut == i)
      
      if (length(pos) == 1) {
        if(i %in% all.clusters) {
          write.table(mega.data.clust[,c("Family","ID")][which(mega.data.clust$groups==pos),],file=paste0(clusters.dir,"cluster.",i,".",name.add,".txt"), quote=F, row.names=F,col.names=F)
          name.add <- name.add + 1
        } else {
          
          write.table(mega.data.clust[,c("Family","ID")][which(mega.data.clust$groups==pos),],file=paste0(clusters.dir,"cluster.",i,".txt"), quote=F, row.names=F,col.names=F) 
        }
      } else {
        files <- dir(clusters.dir,pattern=paste0("^cluster.",i,"."))
        if (length(files) == 0) {
          write.table(mega.data.clust[,c("Family","ID")][which(mega.data.clust$groups==pos[1]),],file=paste0(clusters.dir,"cluster.",i,".",name.add,".txt"), quote=F, row.names=F,col.names=F) 
          name.add <- name.add + 1
          write.table(mega.data.clust[,c("Family","ID")][which(mega.data.clust$groups==pos[2]),],file=paste0(clusters.dir,"cluster.",i,".",name.add,".txt"), quote=F, row.names=F,col.names=F)
          name.add <- name.add + 1
        } else {
          for (j in pos) {
            ctclust <- 0
            temp2 <- mega.data.clust[,c("Family","ID")][which(mega.data.clust$groups==j),]
            for (k in files) {
              temp <- read.table(paste0(clusters.dir,k))
              colnames(temp) <- c("Family","ID")
              compar <- setdiff(levels(droplevels(temp[,1])),levels(droplevels(temp2[,1])))
              if (length(compar)==0) {
                next
                next
              }
              ctclust <- ctclust + 1
            }
            if (ctclust == length(files)) {
              posfin <- j
              break
            }
          }
          
          write.table(mega.data.clust[,c("Family","ID")][which(mega.data.clust$groups==posfin),],file=paste0(clusters.dir,"cluster.",i,".",name.add,".txt"),quote=F, row.names=F,col.names=F)
          name.add <- name.add + 1 
          
        }
        
        
      }
      
    }
  }
}
  all.clusters <- c(all.clusters,new.clusters)
  prev.clusters <- new.clusters
}
