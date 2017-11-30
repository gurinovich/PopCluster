1. all_QC_pruned.bed, all_QC_pruned.bim, all_QC_pruned.fam PLINK files with genome-wide genotype data

To get all_QC_pruned.bed/bim/fam the following QC should be performed on the original PLINK file:

* standard QC: HW equilibrium, remove very rare variants, etc. 
* Exclude the SNPs in the following two regions (can be done using PLINK): 
  * MHC region (chr6)
  * region of inversion polymorphism (chr8)
* Remove AT and CG SNPs.  This can be done in R.  For example,

#### Go to R: ##############################################################
	data <- read.table("all_qc.bim", header=FALSE)
	comb <- paste(data$V5, data$V6, sep="")
	print(table(comb))
	atgc <- (comb %in% c("AT", "TA", "GC", "CG"))
	print(table(atgc))
	write.table( data[ atgc==FALSE, 2], "mylist.txt", quote=F, row.names=F, col.names=F)

	Then, use plink to remove these SNPs from your genotype data.
  
* Prune the genotype data (i.e. only keep independent SNPs).  This is done in two steps. 
First, we get a list of SNPs that we want to keep in the following command: 

	plink --bfile all_qc --indep-pairwise 50 5 0.3 --out all_qc_pruned

Then, we extract that list of SNPs in the following command:

	plink --bfile all_qc --extract all_qc_pruned.prune.in --make-bed --out all_QC_pruned
