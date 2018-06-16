The files described below are required as input files to PopCluster. You will find example simulated files to check that the structure of input files is correct.

1. **all_QC_pruned.bed**, **all_QC_pruned.bim**, **all_QC_pruned.fam** PLINK files with genome-wide genotype data

To get all_QC_pruned.bed/bim/fam the following QC should be performed on the original PLINK files:

* standard QC: HW equilibrium, remove very rare variants, etc. 
* Exclude the SNPs in the following two regions (can be done using PLINK): 
  * MHC region (chr6)
  * region of inversion polymorphism (chr8)
* Remove AT and CG SNPs. This can be done in R.  For example,

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

To create IBD file from your PLINK bed/bim/fam file, use the following PLINK command:

	plink  --bfile all_QC_pruned --genome --min 0.10 --out IBD
	
3. **mega-data-noPCs.csv**

Create **mega-data-noPCs.csv** file with the following columns:

* ID: ids for the subjects (should match the ids in provided PLINK file)
* Family: family ids for the subjects (should match the family ids in provided PLINK file)
* phenotype: 0 for controls and 1 for cases
* columns for covariates (the same names will be used in **covariates.txt** file): covariates to be included in the model (the covariates included in the **covariates.txt** file only will be tested; here we can have more)
* columns for SNPs to be tested (could be more than to be tested, but the names of the SNPs to be tested must match the names in the file **SNPs.txt**): dosages for the SNPs/alleles to be tested by the algorithm
* family.id: in case of related individuals, please provide different numbers for different families. In case of no related individuals, or no knowledge about the relationships, please provide sequential numbers from 1 to the number of subjects. If there are related individuals, the algorithm will account for family structure through generalized estimating equation (GEE).

4. **covariates.txt** file contains one column with covariates' names to be included in the association model

5. **SNPs.txt** file contains one column with variants names to be tested
 
## How we simulated the example files

1. Downloaded **hapmap1.map/ped** files from http://zzz.bwh.harvard.edu/plink/tutorial.shtml on November 30, 2017
2. Converted PLINK map/ped format to bed/bim fam format using the following PLINK command:

```plink --file hapmap1 --make-bed --out all_QC_pruned```

3. Modified **all_QC_pruned.fam** file for it to have unique sample IDs:

#### Go to R: ##############################################################
	library(dplyr)
	library(data.table)
	fam <- fread("all_QC_pruned.fam")
	fam <- tbl_df(fam)
	fam$V2 <- seq(1, nrow(fam))
	write.table(fam, file = "all_QC_pruned.fam", quote = FALSE, row.names = FALSE, col.names = FALSE)

4. Created **IBG.genome** file using the following PLINK command:

```plink --bfile all_QC_pruned --genome --min 0.10 --out IBD```

5. Created **mega-data-noPCs.csv** file:

#### Go to R: ##############################################################
	library(dplyr)
	library(data.table)
	fam <- fread("all_QC_pruned.fam")
	fam <- tbl_df(fam)
	mega.data <- data_frame(ID = fam$V2,
                        Family = fam$V1,
                        phenotype = sample(0:1, nrow(fam), replace = TRUE),
                        covariate1 = sample(1:2, nrow(fam), replace = TRUE),
                        covariate2 = sample(1:2, nrow(fam), replace = TRUE),
                        SNP1 = sample(0:2, nrow(fam), replace = TRUE),
                        SNP2 = sample(0:2, nrow(fam), replace = TRUE),
                        SNP3 = sample(0:2, nrow(fam), replace = TRUE),
			family.id = seq(1, nrow(fam)))
	write.table(mega.data, file = "mega-data-noPCs.csv", quote = FALSE, sep = ",", row.names = FALSE)
