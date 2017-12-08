# PopCluster
PopCluster: a new algorithm to identify genetic variants with ethnicity-dependent effects

Before running the program, make sure that R and PLINK are installed.

To submit the job, run the following shell script that links all the files and scripts:

*sh PopCluster.sh workdir scriptsdir filesdir eigensoftdir*

* workdir: where the PopCluster.sh is located and where you want results to be saved.
* scriptsdir: scripts-dir folder here
* filesdir: files-dir folder here
* eigensoftdir: eigensoft-dir folder here

All the directories should end in "/", for example: "/restricted/projectnb/necs/Nastia-Analysis/clustering-procedure/"

* **scripts-dir** folder contains R scripts
* **files-dir** folder contains input files (see README.md in ./files-dir for details)
* **eigensoft-dir** folder contains EIGENSOFT files for calculating PCA from genome-wide genotype data (see README.md in ./eigensoft-dir for details)
* **results** folder will be created after running the PopCluster.sh and will store the results of the PopCluster
