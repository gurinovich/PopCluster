The files described below are required as input files to PopCluster. You will find example simulated files to check that the structure of input files is correct.

1. **all_QC_pruned.bed**, **all_QC_pruned.bim**, **all_QC_pruned.fam** PLINK files with genome-wide genotype data

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
  
* Prune the genotype data (i.e. only keep independent SNPs).  This is done in two steps. all_qc.bed/bim/fam is intermediate PLINK file. 

First, we get a list of SNPs that we want to keep in the following command: 

	plink --bfile all_qc --indep-pairwise 50 5 0.3 --out all_qc_pruned

Then, we extract that list of SNPs in the following command:

	plink --bfile all_qc --extract all_qc_pruned.prune.in --make-bed --out all_QC_pruned

2. **IBD.genome**: identity by descent file.

To create IBD file, we use the following PLINK command:

	plink  --file all_QC_pruned --genome --min 0.10 --out IBD-file
	
3. **mega-data-noPCs.csv**

Create **mega-data-noPCs.csv** file with the following columns:

* ID: ids for the subjects (should match the ids in provided PLINK file)
* Family: family ids for the subjects (should match the family ids in provided PLINK file)
* phenotype: 0 for controls and 1 for cases
* columns for covariates (the same names will be used in **covariates.txt** file): covariates to be included in the model (the covariates included in the **covariates.txt** file only will be tested; here we can have more)
* columns for SNPs to be tested (could be more than to be tested, but the names of the SNPs to be tested must match the names in the file **SNPs.txt**): dosages for the SNPs/alleles to be tested by the algorithm

4. **covariates.txt** file contains one column with covariates' names to be included in the association model

5. **SNPs.txt** file contains one column with variants names to be tested
 

  
