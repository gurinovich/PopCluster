# PopCluster
PopCluster: a new algorithm to identify genetic variants with ethnicity-dependent effects

## Dependencies
* Before running the program, make sure that R, PLINK, and EIGENSOFT (version 6.1.4) are installed.
* Load respective modules:

*module load EIGENSOFT*

*module load R*

*module load plink*

* Install R dplyr package if not installed already:

in R:

*install.packages("dplyr")*

* For files in eigensoft-dir folder, make sure the priviliges are set up appropriately. In eigensoft-dir:

*chmod a=rwx convertf*

*chmod a=rwx smartpca*

## Execution
To submit the job, run the following shell script:

*sh PopCluster.sh arguments.txt*

**arguments.txt** file contains arguments that include paths that need to be changed and default values for parameters that can be changed as well:

* WORKDIR: where the PopCluster.sh is located and where you want results to be saved.
* SCRIPTSDIR: scripts-dir folder here (contains R scripts)
* FILESDIR: files-dir folder here (contains input files (see README.md in ./files-dir for details))
* EIGENSOFTDIR: eigensoft-dir folder here (contains EIGENSOFT files for calculating PCA from genome-wide genotype data (see README.md in ./eigensoft-dir for details))
* ClusteringPCs: number of principal components used for clustering (default: 6)
* maxclustsize: minimum cluster sizet to be included in the analysis (default: 100)
* GLMPCs: number of principal components included in the logistic regression (default: 4)
* GLMfamily: link function to be used in the generalized linear model (default: binomial - for logistic regression)

All the directories should be complete paths and should end in "/", for example: "/restricted/projectnb/necs/Nastia-Analysis/clustering-procedure/".

* **results** folder will be created after running the PopCluster.sh and will store the results of the PopCluster
