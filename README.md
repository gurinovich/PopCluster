# PopCluster
PopCluster: a new algorithm to identify genetic variants with ethnicity-dependent effects

To submit the job, run the following shell script that links all the files and scripts:

*sh PopCluster.sh workdir scriptsdir filesdir eigensoftdir*

* workdir: where the PopCluster.sh is located and where you want results to be saved.
* scriptsdir: scripts-dir folder here
* filesdir: files-dir folder here
* eigensoftdir: eigensoft-dir folder here

All the directories should end in "/", for example: "/restricted/projectnb/necs/Nastia-Analysis/clustering-procedure/"

* **scripts-dir** folder contains R scripts
* **files-dir** folder contains input files (see below for details)
* **eigensoft-dir** folder contains EIGENSOFT files for calculating PCA from genome-wide genotype data
* **results** folder will contain the results of the PopCluster
