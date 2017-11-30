args = commandArgs(trailingOnly=TRUE)

work.dir <- args[1]
files.dir <- paste0(work.dir,"full-summary/")
clust.dir <- paste0(work.dir, "clusters/")
cluster.member <- read.csv(paste0(work.dir,"cluster_membership.txt"),colClasses="character")
signif.level <- 0.05

# ----------------------- START FUNCTIONS ----------------------------

find.sibling <- function(clust,cluster.member) {
  parent <- who.is.parent(clust,cluster.member)
  if (is.na(parent)) return(NA)
  siblings <- find.children(parent,cluster.member)
  sibling <- siblings[-which(siblings %in% clust)]  
  return(sibling)
}

find.children <- function(clust,cluster.member) {
  if (cluster.member[which(cluster.member$Parent == clust),c("Cluster1")] == ".") {
    return(NA)
  }
  child1 <- cluster.member[which(cluster.member$Parent == clust),c("Cluster1")]
  child2 <- cluster.member[which(cluster.member$Parent == clust),c("Cluster2")]
  return(c(child1, child2))
}

#calculate z-score:
z.score <- function(betta1, se1, betta2, se2) {
  z.score <- abs(betta1-betta2) / sqrt(se1^2 + se2^2)
  return(z.score)
}

#TRUE if bettas are the same, FALSE if bettas are different:
same.bettas <- function(betta1, se1, betta2, se2) {
    z.score <- z.score(betta1, se1, betta2, se2)
    return( ifelse(z.score < 2, TRUE, FALSE))
}

who.is.parent <- function(sib1, cluster.member) {
  ifelse(sib1 %in% cluster.member$Cluster1,parent <- cluster.member$Parent[which(cluster.member$Cluster1 == sib1)],parent <- cluster.member$Parent[which(cluster.member$Cluster2 == sib1)])
  if(length(parent) == 0) return(NA)
  return(parent)
}


find.sibling.leaves.with.parents <- function(cluster.member) {
  empty.vec <- vector()
  leaves <- cluster.member[which(cluster.member$Cluster1 == "." & cluster.member$Cluster2 == "."),c("Parent")]
  if (length(leaves) == 0) return(empty.vec)
  parents <- vector()
  for (i in 1:length(leaves)) {
    parents[i] <- who.is.parent(leaves[i], cluster.member)
 }
  parents <- parents[!is.na(parents)]
 if (length(parents) == 0) return(empty.vec) 
  duplicated.parents <- parents[duplicated(parents)]
 if (length(duplicated.parents) == 0) return(empty.vec) 
  leaves.with.parents <- vector()
  for (i in 1:length(duplicated.parents)) {
   leaves.with.parents[2*i-1] <- find.children(duplicated.parents[i],cluster.member)[1]
   leaves.with.parents[2*i] <- find.children(duplicated.parents[i],cluster.member)[2]   
 }
 return(leaves.with.parents)
}

merge.into.parent <- function(sib1,sib2,cluster.member) {
  parent <- who.is.parent(sib1, cluster.member)
  cluster.member <- cluster.member[-c(which(cluster.member$Parent == sib1),which(cluster.member$Parent == sib2),which((cluster.member$Cluster1 == sib1 & cluster.member$Cluster2 == sib2) | (cluster.member$Cluster2 == sib1 & cluster.member$Cluster1 == sib2))),]  
  cluster.member <- rbind(cluster.member,c(parent,".","."))
  return(cluster.member)
}


report.children <- function(sib1,sib2,cluster.member) {
  to.remove <- vector()
#  top.cluster <- max(cluster.member$Parent)
#  while(who.is.parent(sib1,cluster.member) != top.cluster & length(who.is.parent(sib1,cluster.member)) != 0) {
  while(length(who.is.parent(sib1,cluster.member)) != 0 & !is.na(who.is.parent(sib1,cluster.member))) {
    to.remove <- c(to.remove,who.is.parent(sib1,cluster.member))
    sib1 <- who.is.parent(sib1,cluster.member)
  }  
  if(length(to.remove) != 0) cluster.member <- cluster.member[-which(cluster.member$Parent %in% to.remove | cluster.member$Cluster1 %in% to.remove | cluster.member$Cluster2 %in% to.remove),]  
  return(cluster.member)
}

bottom.up <- function(cluster.member,input,signif.level) {

#1) find pairs of leaves that have siblings who are also leaves and they have parent
  sibling.leaves.with.parents <- find.sibling.leaves.with.parents(cluster.member)
  
  if(length(sibling.leaves.with.parents) == 0) {return(cluster.member)}

  for (i in 1:(length(sibling.leaves.with.parents)/2)) {
    sib1 <- sibling.leaves.with.parents[2*i-1]
    sib2 <- sibling.leaves.with.parents[2*i]
    sib1.p <- input[which(input$Clusters == sib1),c("GLM.p")]
    sib2.p <- input[which(input$Clusters == sib2),c("GLM.p")]
    sib1.b <- input[which(input$Clusters == sib1),c("GLM.b")]
    sib2.b <- input[which(input$Clusters == sib2),c("GLM.b")]
    sib1.se <- input[which(input$Clusters == sib1),c("GLM.SE")]
    sib2.se <- input[which(input$Clusters == sib2),c("GLM.SE")]
    
 
      if (sib1.p < signif.level | sib2.p < signif.level) {      
        if (same.bettas(sib1.b,sib1.se,sib2.b,sib2.se)) {
          cluster.member <- merge.into.parent(sib1,sib2,cluster.member)
        } else {
          cluster.member <- report.children(sib1,sib2,cluster.member)
        }
      } else {
        cluster.member <- merge.into.parent(sib1,sib2,cluster.member)
      }
    }  
#  }

  cluster.member <- bottom.up(cluster.member,input,signif.level)  
  return(cluster.member)
}


# ----------------------- END FUNCTIONS ----------------------------

file.names <- dir(files.dir)

#output <- data.frame(matrix(NA,nrow=1,ncol=18))
#names(output) <- c("Clusters","PC1.b","PC1.p","PC1.SE","PC2.b","PC2.p","PC2.SE","PC3.b","PC3.p","PC3.SE", "PC4.b","PC4.p","PC4.SE", "GLM.b","GLM.p", "GLM.SE","Allele","Clust.sib")

output <- data.frame(matrix(NA,nrow=1,ncol=6))
names(output) <- c("Clusters","GLM.b","GLM.p", "GLM.SE","Allele","Clust.sib")


for (i in 1:length(file.names)) {
  cluster.member.t <- cluster.member
  input <- read.csv(paste0(files.dir,file.names[i]),colClasses=c("Clusters"="character"))
  
# check if there are NA entries and adjust input & cluster.member accordingly
  if (any(apply(input,2,function(x) any(is.na(x))))) {
    clusts.na <- input[which(is.na(input$GLM.b)),c("Clusters")]
    siblings <- vector()
    for (j in 1:length(clusts.na)) {
      siblings <- c(siblings,find.sibling(clusts.na[j],cluster.member.t))
    }
    clusts.na <- unique(c(clusts.na,siblings))
    leaves <- cluster.member.t[which(cluster.member.t$Cluster1 == "." & cluster.member.t$Cluster2 == "."),c("Parent")]
    #if all the NA entries are leaves:
    if (length(intersect(clusts.na,leaves)) == length(clusts.na)) {
      for (j in 1:length(clusts.na)) {
        if (any(cluster.member.t$Parent %in% clusts.na[j])) {
        cluster.member.t <- merge.into.parent(find.children(who.is.parent(clusts.na[j],cluster.member.t),cluster.member.t)[1],find.children(who.is.parent(clusts.na[j],cluster.member.t),cluster.member.t)[2],cluster.member.t)
        } else next
      }
    } else next
  }
  
  
  t1 <- bottom.up(cluster.member.t,input,signif.level)
  colnames(t1) <- c("Clusters","Cluster1","Cluster2")
#  sum(as.numeric(t1$Clusters))
#  sum(as.numeric(as.character(t1$Clusters)))

  if (nrow(t1) == 1) {
    clust <- as.character(t1$Clusters)
    clust.p <- input[which(input$Clusters == clust),c("GLM.p")]
    if (clust.p >= signif.level) next
  }


  t2 <- merge(t1,input)
  t2 <- t2[,-which(names(t2) %in% c("Cluster1","Cluster2"))]
  t2$Allele <- rep(sub(".txt","",file.names[i]),nrow(t2))

  clust.temp <- t2$Clusters 
  sib.temp <- vector()
  for (j in 1:length(clust.temp)) {
    sib.temp <- c(sib.temp, find.sibling(clust.temp[j],cluster.member))
    if (!(sib.temp[j] %in% clust.temp)) sib.temp[j] <- NA
  }

  t2$Clust.sib <- sib.temp

  output <- rbind(output,t2)
  
}

output <- output[-which(is.na(output$Allele)),]

output$norm.betta <- output$GLM.b/output$GLM.SE
#output <- output[,c("Allele","Clusters","Clust.sib","norm.betta", "GLM.p","GLM.b","GLM.SE","PC1.b",      "PC1.p","PC1.SE","PC2.b","PC2.p","PC2.SE","PC3.b","PC3.p","PC3.SE","PC4.b","PC4.p","PC4.SE")]
output <- output[,c("Allele","Clusters","Clust.sib","norm.betta", "GLM.p","GLM.b","GLM.SE")]

write.table(output,paste0(work.dir,"results-0.05.csv"),quote=F,sep=",",row.names=F)
