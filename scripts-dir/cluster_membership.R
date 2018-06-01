args = commandArgs(trailingOnly=TRUE)

work.dir <- args[1]
res.dir <- paste0(work.dir, "results/")
clust.dir <- paste0(res.dir, "clusters/")

clust.files <- dir(clust.dir)
clust.names <- gsub(".txt","",clust.files)

#For each cluster create a vector with cluster names that are included in the cluster:
for (i in 1:length(clust.files)) {
  assign(clust.names[i],vector())
  clust.curr <- read.csv(paste0(clust.dir,clust.files[i]),sep=" ",header=F)
  names(clust.curr) <- c("Family","ID")
  
  for (j in 1:length(clust.files)) {
    if (i == j) next
    clust.check <- read.csv(paste0(clust.dir,clust.files[j]),sep=" ",header=F)
    names(clust.check) <- c("Family","ID")
    check.var <- setdiff(clust.check$ID,clust.curr$ID)
    if (length(check.var) == 0) assign(clust.names[i],c(get(clust.names[i]),clust.names[j]))
  }
}

output <- data.frame(matrix(NA,nrow=length(clust.files),ncol=3))
names(output) <- c("Parent","Cluster1","Cluster2")
output$Parent <- clust.names

for (i in 1:length(clust.files)) {
  t1 <- sort(gsub("cluster.","",get(clust.names[i])),decreasing=T)
  if (length(t1) == 0) {
    next
  } else if (length(t1) == 1) {
    output$Cluster1[which(output$Parent == clust.names[i])] = t1[1]
  } else {
    clust.curr <- read.csv(paste0(clust.dir,clust.files[i]),sep=" ",header=F)
    names(clust.curr) <- c("Family","ID")
    for (j in 1:(length(t1)-1)) {
      check.var <- T
      clust.check <- read.csv(paste0(clust.dir,paste0("cluster.",t1,".txt")[j]),sep=" ",header=F)
      names(clust.check) <- c("Family","ID")
      for (k in (j+1):length(t1)) {
        clust.check2 <- read.csv(paste0(clust.dir,paste0("cluster.",t1,".txt")[k]),sep=" ",header=F)
        names(clust.check2) <- c("Family","ID")
        if (identical(sort(union(clust.check$ID,clust.check2$ID)),sort(as.vector(clust.curr$ID)))) {
          output$Cluster1[which(output$Parent == clust.names[i])] = t1[j]
          output$Cluster2[which(output$Parent == clust.names[i])] = t1[k]
          check.var <- F
          break
        }
      }
      if (!check.var) break
    }
    
    if (check.var) {
      check.var2 <- "CHECK."
      for (j in 1:length(t1)) {
        check.var2 <- paste0(check.var2,t1[j],"-")
      }
      check.var2 <- substr(check.var2,1,nchar(check.var2)-1)
      output$Cluster1[which(output$Parent == clust.names[i])] = check.var2
    }
    
    
  }
}    

children <- c(output$Cluster1,output$Cluster2)
children <- children[!is.na(children)]
children.used <- children[-grep("CHECK",children)]
all.clusters <- gsub("cluster.","",clust.names)
children.not.used <- setdiff(all.clusters,children.used)

index <- row.names(output[grep("CHECK",output$Cluster1),])

for (i in index) {
  j <- as.numeric(i)
  t1 <- unlist(strsplit(gsub("CHECK.","",output$Cluster1[j]),"-"))
  
  if (length(intersect(t1,children.not.used)) == 1) {
    output$Cluster1[j] <- intersect(t1,children.not.used)[1]
    children.not.used <- children.not.used[children.not.used != intersect(t1,children.not.used)[1]]
  } else {
    output$Cluster1[j] <- paste(intersect(t1,children.not.used),collapse="-")
  } 
}

index <- row.names(output[grep("-",output$Cluster1),])

for (i in index) {
  j <- as.numeric(i)
  t1 <- unlist(strsplit(output$Cluster1[j],"-"))
  
  if (length(intersect(t1,children.not.used)) == 1) {
    output$Cluster1[j] <- intersect(t1,children.not.used)[1]
    children.not.used <- children.not.used[children.not.used != intersect(t1,children.not.used)[1]]
  } else {
    output$Cluster1[j] <- paste(intersect(t1,children.not.used),collapse="-")
  } 
}

output$Parent <- gsub("cluster.","",clust.names)

###############
## identify clusters that don't have siblings for removal:
to.remove <- vector()
to.remove <- output[which(!is.na(output$Cluster1) & is.na(output$Cluster2)),"Cluster1"]
to.remove.keep <- to.remove

if (!length(to.remove) == 0) {
  for (i in to.remove) {
    temp.clust <- get(paste0("cluster.",i))
    to.remove.keep <- c(to.remove.keep,temp.clust)
  }
}

to.remove.keep <- sub("cluster.","",to.remove.keep)
to.remove.keep <- unique(to.remove.keep)

# remove clusters files from the clusters folder that are leaves without siblings:
files.to.remove <- sub("^","cluster.",to.remove.keep)
files.to.remove <- sub("$",".txt",files.to.remove)

for (i in files.to.remove) {
  system(paste0("rm ",clust.dir,i))
}

output.temp <- output
if (length(to.remove.keep) !=0 ) {
  output.temp <- output.temp[-which(output.temp$Parent %in% to.remove.keep),]
  if (nrow(output.temp) != 1) {
    output.temp <- output.temp[-which(output.temp$Cluster1 %in% to.remove.keep),]
  } else {
    output.temp$Cluster1 <- NA
    output.temp$Cluster2 <- NA
  }
}

# add leave nodes that are not present yet: cluster, . , .
all.clusters <- c(output.temp$Parent,output.temp$Cluster1,output.temp$Cluster2)
all.clusters <- all.clusters[-which(is.na(all.clusters))]
duplicated <- all.clusters[which(duplicated(all.clusters))]
addit.leaves <- all.clusters[which(!all.clusters %in% duplicated)]
#remove top parent cluster:
top.parent <- max(as.numeric(all.clusters[which(!all.clusters %in% duplicated)]))
addit.leaves <- addit.leaves[-which(addit.leaves == top.parent)]

for (i in addit.leaves) {
  row <- c(i,".",".")
  output.temp <- rbind(output.temp,row)
}

if (nrow(output.temp) == 1 & is.na(output.temp$Cluster1) & is.na(output.temp$Cluster2)) {
  output.temp <- output.temp[-1,]
}

write.table(output.temp,paste0(res.dir,"cluster_membership.txt"), quote=F, sep=",",na=".",row.names=F)
