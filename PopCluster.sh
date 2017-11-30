#!/bin/bash

#To run type qsub -l h_rt=240:00:00 -P necs -b y -cwd -j y -l mem_free=64g -m ae -M agurinov@bu.edu sh PopCluster.sh

module load EIGENSOFT
module load R
module load plink

# Things to change - INPUT ARGUMENTS:

#DIRECTORIES:
#working directory (where PopClsuter.sh is located & all other folders (scripts-dir and files-dir) & files)
WORKDIR=/restricted/projectnb/necs/Nastia-Analysis/clustering-procedure/PopCluster/
SCRIPTSDIR=$WORKDIR"scripts-dir/"
FILESDIR=$WORKDIR"files-dir/"
#EIGENSOFTDIR=/share/apps/6.0/EIG/EIG5.0.1_gnu447_x86_64/bin/
#EIGENSOFTDIR: where convertf.c and smartpca.c are located
EIGENSOFTDIR=/restricted/projectnb/necs/Nastia-Analysis/clustering-procedure/EIG-7.2.1/src/

#FILES:
#bed/bim/fam file that consists all the subjects to analyze (+ there could be more)(Plink format - QCd and imputed)
genotypeData=$FILESDIR"all_QC_pruned"
IBD=$FILESDIR"IBD.7.18.2015.genome"
#mega-file without PCs with all the subjects and SNPs to analyze (subjects must be present in $genotypeData file) and with the covariates and binary phenotype info
megaData=$FILESDIR"mega-data-noPCs.csv"
covariates=$FILESDIR"covariates.txt"
SNPs=$FILESDIR"SNPs.txt"

# The code (not to change):
Rscript $SCRIPTSDIR"PCA-top-cluster.R" $WORKDIR $EIGENSOFTDIR $megaData $genotypeData $IBD --save
Rscript $SCRIPTSDIR"clustering.R" $WORKDIR --save
Rscript $SCRIPTSDIR"cluster_membership.R" $WORKDIR --save

# Calculate PCA:
RESULTSDIR=$WORKDIR"results/"
mkdir $RESULTSDIR"PCA"
mkdir $RESULTSDIR"analysis"
CLUSTERSDIR=$RESULTSDIR"clusters/"
PCADIR=$RESULTSDIR"PCA/"

prefix="$CLUSTERSDIR"
prefix+="cluster."
suffix=".txt"

for file in $CLUSTERSDIR*
do
  clust=${file%$suffix}
  clust=${clust#$prefix}
  echo "${clust}"
  DIRECTORY="$PCADIR"
  DIRECTORY+="$clust"
  if [ ! -d "$DIRECTORY" ]
  then
     Rscript $SCRIPTSDIR"calculate-PCA.R" "$clust" $RESULTSDIR $genotypeData $EIGENSOFTDIR $IBD
  fi
done

#PopCluster
Rscript $SCRIPTSDIR"create_summary.R" $RESULTSDIR $WORKDIR $SNPs $covariates --save
Rscript $SCRIPTSDIR"combine-results.R" $RESULTSDIR --save