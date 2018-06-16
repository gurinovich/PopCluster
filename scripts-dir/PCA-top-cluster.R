args = commandArgs(trailingOnly=TRUE)

work.dir <- args[1]
system(paste0("mkdir ", work.dir,"PCA-topcluster"),wait=T)
pca.dir <- paste0(work.dir,"PCA-topcluster/")
eigensoft.dir <- args[2]
mega.data <- args[3]
genotype.data <- args[4]
IBD.file <- args[5]

#QC:
system(paste0("mkdir ", pca.dir,"QC"),wait=T)
setwd(paste0(pca.dir,"QC"))
mega.data <- read.csv(mega.data)
mega.data.temp <- mega.data[,c("Family","ID")]
cluster <- nrow(mega.data.temp)
write.table(mega.data.temp,file=paste0("cluster.",cluster,".txt"),quote=FALSE,row.names=FALSE,col.names=FALSE)
system(paste0("plink --bfile ", genotype.data, " --keep cluster.", cluster, ".txt --make-bed --out all_QC_pruned"),wait=T)
genotype.data <- tail(strsplit(genotype.data, split = "/")[[1]], n = 1)
fam <- read.table(paste0(genotype.data, ".fam"))
fam[,6]  <- 1
write.table(fam, paste0(genotype.data, "_edit.fam"), quote = F, col.names = F, row.names = F)

#Convert:
system(paste0("mkdir ", pca.dir, "Convert"),wait=T)
setwd(paste0(pca.dir, "Convert"))
file.create("par.PED.EIGENSTRAT")
fileConn <- file("par.PED.EIGENSTRAT", "w")
writeLines(c(paste0("genotypename: ", pca.dir, "QC/", genotype.data, ".bed #.ped, .bed "),
             paste0("snpname: ", pca.dir, "QC/", genotype.data, ".bim  # .map or .bim, either works "),
             paste0("indivname: ", pca.dir, "QC/", genotype.data, "_edit.fam # .ped or .fam, either works"),
             "outputformat: EIGENSTRAT",
             "genotypeoutname: all.eigenstratgeno",
             "snpoutname: all.snp",
             "indivoutname: all.ind",
             "familynames: NO"), fileConn)
close(fileConn)
system(paste0(eigensoft.dir, "convertf -p par.PED.EIGENSTRAT > convert.log"), wait=T)
rm(fileConn)

#Eigen:
system(paste0("mkdir ", pca.dir, "Eigen"),wait=T)
setwd(paste0(pca.dir, "Eigen"))
ind <- read.table(paste0(pca.dir, "Convert/all.ind"))
ibd <- read.table(IBD.file, header = T)
which.ibd.keep <- which(as.character(ibd[,2]) %in% ind[,1] & 
                          as.character(ibd[,4]) %in% ind[,1] & ibd$PI_HAT > 0.2)
length(which.ibd.keep)
if (length(which.ibd.keep) == 0) {
  ind[,3] <- "Phenotype"  
} else {
  ibd.keep <- ibd[which.ibd.keep,]
  table(table(c( as.character(ibd.keep[,2]), as.character(ibd.keep[,4]))))       
  related.ids <- c()
  for (i in 1:dim(ibd.keep)[1]){
    if ( !(as.character(ibd.keep[i,2]) %in% related.ids) & 
           !(as.character(ibd.keep[i,4]) %in% related.ids)){
      ran.index <- sample(c(2,4),1)
      related.ids <- c(related.ids, as.character(ibd.keep[i,ran.index]))
    }
  }
  write.csv(related.ids, "Related.IDS.csv", quote = F, row.names = F)
  which.rel <- which(ind[,1] %in% related.ids)
  length(which.rel)
  ind[,3] <- "Phenotype"
  ind[which.rel,3] <- "Related_InferredPC"
}
write.table(ind , "phenotype_PCA.ind", quote = F, col.names = F, row.names = F, sep = "    ")
## making a list of populations to be analyzed (don't put related in there b/c you want their eigenvectors inferred)
pop <- names(table(ind[,3]))
pop <- pop[!(pop %in% c("Ignore", "Related_InferredPC","BorderlineCR_InferredPC" ))]
write.table(pop, "pop_list.txt", quote = F, col.names = F, row.names = F)

file.create("smartpca.par")
fileConn <- file("smartpca.par", "w")
writeLines(c(paste0("genotypename: ", pca.dir, "Convert/all.eigenstratgeno"),
             paste0("snpname: " ,pca.dir, "Convert/all.snp"),
             paste0("indivname: ", pca.dir, "Eigen/phenotype_PCA.ind"),           
             "evecoutname:     phenotype.evec",
             "evaloutname:     phenotype.eval",
             "numoutevec:      30",
             "numoutlierevec:  5",
             "outliersigmathresh:  9",
             paste0("poplistname: ", pca.dir, "Eigen/pop_list.txt"),
             "snpweightoutname: phenotype_snp_wt.txt"), fileConn)
close(fileConn)

system(paste0(eigensoft.dir, "smartpca -p smartpca.par > pca.log"), wait=T)
system(paste0("sed '1d' phenotype.evec | tr -s ' ' | sed -e 's/^[ \t]*//' | cut -f 1-31 -d \" \" > top-cluster.txt"), wait=T)

# add recalculated PCs to mega.data (some subjects that are too outliers might get removed)
setwd(work.dir)
pcs <- read.table(paste0(pca.dir,"/Eigen/top-cluster.txt"))
colnames(pcs) <- c("ID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8")
pcs <- pcs[,c(1:9)]
temp <- merge(mega.data,pcs)
temp <- temp[,c(colnames(mega.data),"PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8")]

write.table(temp,file=paste0(work.dir,"mega-data.csv"),quote=F, sep=",",row.names=FALSE)
