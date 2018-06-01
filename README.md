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
To submit the job, run the following shell script that links all the files and scripts:

*sh PopCluster.sh workdir scriptsdir filesdir eigensoftdir*

* workdir: where the PopCluster.sh is located and where you want results to be saved.
* scriptsdir: scripts-dir folder here
* filesdir: files-dir folder here
* eigensoftdir: eigensoft-dir folder here

All the directories should be complete paths and should end in "/", for example: "/restricted/projectnb/necs/Nastia-Analysis/clustering-procedure/".

Example of running:

*sh PopCluster.sh /restricted/projectnb/necs/Nastia-Analysis/clustering-procedure/PopCluster/PopCluster-test2/ /restricted/projectnb/necs/Nastia-Analysis/clustering-procedure/PopCluster/PopCluster-test2/scripts-dir/ /restricted/projectnb/necs/Nastia-Analysis/clustering-procedure/PopCluster/PopCluster-test2/files-dir/ /restricted/projectnb/necs/Nastia-Analysis/clustering-procedure/PopCluster/PopCluster-test2/eigensoft-dir/*

* **scripts-dir** folder contains R scripts
* **files-dir** folder contains input files (see README.md in ./files-dir for details)
* **eigensoft-dir** folder contains EIGENSOFT files for calculating PCA from genome-wide genotype data (see README.md in ./eigensoft-dir for details)
* **results** folder will be created after running the PopCluster.sh and will store the results of the PopCluster
