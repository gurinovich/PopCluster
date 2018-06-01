#!/bin/bash

source $1

#FILES:
#bed/bim/fam file that consists all the subjects to analyze (+ there could be more)(Plink format - QCd and imputed)
genotypeData=$FILESDIR"all_QC_pruned"
IBD=$FILESDIR"IBD.genome"
#mega-file without PCs with all the subjects and SNPs to analyze (subjects must be present in $genotypeData file) and with the covariates and binary phenotype info
megaData=$FILESDIR"mega-data-noPCs.csv"
covariates=$FILESDIR"covariates.txt"
SNPs=$FILESDIR"SNPs.txt"

# The code (not to change):
Rscript $SCRIPTSDIR"PCA-top-cluster.R" $WORKDIR $EIGENSOFTDIR $megaData $genotypeData $IBD --save
Rscript $SCRIPTSDIR"clustering.R" $WORKDIR $ClusteringPCs $maxclustsize --save
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
Rscript $SCRIPTSDIR"create_summary.R" $RESULTSDIR $WORKDIR $SNPs $covariates $GLMPCs $GLMfamily --save
Rscript $SCRIPTSDIR"combine-results.R" $RESULTSDIR --save
