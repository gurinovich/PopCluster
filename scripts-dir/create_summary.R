args = commandArgs(trailingOnly=TRUE)

work.dir <- args[1]
system(paste0("mkdir ", work.dir, "full-summary"))
res.dir <- paste0(work.dir,"full-summary/")
analysis.dir <- paste0(work.dir,"analysis/")
clusters.dir <- paste0(work.dir,"clusters/")
files.dir <- args[2]
PCsNumber <- 4
mega.data <- read.csv(paste0(files.dir,"mega-data.csv"))
mega.data <- subset(mega.data,select=-c(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8))
snps <- read.csv(args[3],header=F,as.is=T)
alleles <- snps$V1
covariates <- read.csv(args[4],header=F,as.is=T)
covariates <- covariates$V1
cov.string <- paste0(paste(covariates, collapse="+"),"+")

file.names <- dir(analysis.dir)

for (i in 1:length(alleles)) {
  output <- data.frame(matrix(NA,nrow=length(file.names),ncol=16))
#  names(output) <- c("Clusters","PC1.b","PC1.p","PC1.SE","PC2.b","PC2.p","PC2.SE","PC3.b","PC3.p","PC3.SE","PC4.b","PC4.p","PC4.SE","GLM.b","GLM.p","GLM.SE")
  names(output) <- c("Clusters","GLM.b","GLM.p","GLM.SE")
  output$Clusters <- gsub(".txt","",file.names)
  
  for (j in 1:length(file.names)) {
    
    cluster <- read.table(paste0(analysis.dir,file.names[j]))
    pca.names <- c("ID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30") 
    colnames(cluster) <- pca.names
    mega.data.cluster <- merge(cluster,mega.data,by="ID")  
    # 
    # for (k in 1:PCsNumber) {
    #   mod <- lm(paste0("PC",k,"~",cov.string,alleles[i]), data=mega.data.cluster)
    #   output[j,c(paste0("PC",k,".b"))] <-  as.numeric(coef(summary(mod))[,1][nrow(coef(summary(mod)))])  #betta
    #   output[j,c(paste0("PC",k,".p"))] <-  as.numeric(coef(summary(mod))[,4][nrow(coef(summary(mod)))])  #p-value
    #   output[j,c(paste0("PC",k,".SE"))] <- as.numeric(coef(summary(mod))[,2][nrow(coef(summary(mod)))]) # SE for the betta
    # }
    # 
    mod <- glm(paste0("phenotype~",cov.string,"PC1+PC2+PC3+PC4+",alleles[i]), family=binomial, data=mega.data.cluster)
    output[j,c("GLM.b")] <-  as.numeric(coef(summary(mod))[,1][nrow(coef(summary(mod)))])  #betta
    output[j,c("GLM.p")] <-  as.numeric(coef(summary(mod))[,4][nrow(coef(summary(mod)))])  #p-value
    output[j,c("GLM.SE")] <-  as.numeric(coef(summary(mod))[,2][nrow(coef(summary(mod)))])   #SE for the betta
    
  }
  write.table(output,paste0(res.dir,alleles[i],".txt"),quote=F,sep=",",row.names=F) 
}
