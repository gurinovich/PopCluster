library(dplyr)
library(testit)
library(geepack)

args = commandArgs(trailingOnly=TRUE)

work.dir <- args[1]
system(paste0("mkdir ", work.dir, "full-summary"))
res.dir <- paste0(work.dir,"full-summary/")
analysis.dir <- paste0(work.dir,"analysis/")
clusters.dir <- paste0(work.dir,"clusters/")
files.dir <- args[2]
mega.data <- read.csv(paste0(files.dir,"mega-data.csv"))
mega.data <- subset(mega.data,select=-c(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8))
snps <- read.csv(args[3],header=F,as.is=T)
alleles <- snps$V1
covariates <- read.csv(args[4],header=F,as.is=T)
covariates <- covariates$V1
cov.string <- paste0(paste(covariates, collapse="+"),"+")
PCsNumber <- as.numeric(args[5])
pcs <- vector()
for (i in 1:PCsNumber) {
    pcs <- c(pcs,paste0("PC",i))
}
PCs.string <- paste0(paste(pcs, collapse="+"),"+")
family.type <- args[6]
Nsubjects <- as.numeric(args[7])
MAF <- as.numeric(args[8])

file.names <- dir(analysis.dir)

for (i in 1:length(alleles)) {
    output <- data.frame(matrix(NA,nrow=length(file.names),ncol=4))
    names(output) <- c("Clusters","GLM.b","GLM.p","GLM.SE")
    output$Clusters <- gsub(".txt","",file.names) 
    for (j in 1:length(file.names)) {     
        cluster <- read.table(paste0(analysis.dir,file.names[j]))
        pca.names <- c("ID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30") 
        colnames(cluster) <- pca.names
        mega.data.cluster <- merge(cluster,mega.data,by="ID")  
        mega.data.cluster <- as.tbl(mega.data.cluster)
        
        table.test <- table(mega.data.cluster$phenotype)
        
        alleles.vec <- select(mega.data.cluster, alleles[i])
        AF <- sum(alleles.vec, na.rm = T)/(2*sum(complete.cases(alleles.vec)))
        MAF.t <- ifelse(AF < 0.5, AF, 1-AF)
        
        mega.data.cluster <- arrange(mega.data.cluster, family.id)
        
        if (table.test[1] >= Nsubjects & table.test[2] >= Nsubjects & MAF.t > MAF) {        
            if (!has_warning(mod <- geeglm(formula(paste0("phenotype~",cov.string, PCs.string, alleles[i])), family=family.type, id = family.id, data=mega.data.cluster, corstr = "exchangeable")))
              { 
          output[j,c("GLM.b")] <-  as.numeric(coef(summary(mod))[,1][nrow(coef(summary(mod)))])  #betta
          output[j,c("GLM.p")] <-  as.numeric(coef(summary(mod))[,4][nrow(coef(summary(mod)))])  #p-value
          output[j,c("GLM.SE")] <-  as.numeric(coef(summary(mod))[,2][nrow(coef(summary(mod)))])   #SE for the betta   
            } 
        }
        
    }
    write.table(output,paste0(res.dir,alleles[i],".txt"),quote=F,sep=",",row.names=F) 
}
