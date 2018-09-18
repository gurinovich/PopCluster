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
wait
Rscript $SCRIPTSDIR"clustering.R" $WORKDIR $ClusteringPCs $maxclustsize --save
wait
Rscript $SCRIPTSDIR"cluster_membership.R" $WORKDIR --save
wait

# Calculate PCA:
RESULTSDIR=$WORKDIR"results/"
mkdir $RESULTSDIR"PCA"
mkdir $RESULTSDIR"analysis"
CLUSTERSDIR=$RESULTSDIR"clusters/"
PCADIR=$RESULTSDIR"PCA/"

prefix="$CLUSTERSDIR"
prefix+="cluster."
suffix=".txt"

wait

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

wait

#PopCluster
Rscript $SCRIPTSDIR"create_summary.R" $RESULTSDIR $WORKDIR $SNPs $covariates $GLMPCs $GLMfamily $Nsubjects $MAF --save
wait
Rscript $SCRIPTSDIR"combine-results.R" $RESULTSDIR $Pvalue --save

{
  echo 'results-0.05.csv file contains the results of PopCluster.'
  echo 'Columns: Allele: name of the SNP/allele; Clusters: name of the cluster (represents the size of the cluster); Clust.sib: sibling cluster (if returned); OR: odds ratio of associaiton; GLM.p: p-value of associaiton; GLM.b: estimate from regression model; GLM.SE: standard error of the estimate from the regression model.'
  echo ''
  echo 'cluster_membership.txt file has a tree structure of the final results.'
  echo ''
  echo 'PCA folder contains PCA results of the returned clusters.'
  echo ''
  echo 'full-summary folder contains complete regression results of the returned associations.'
  echo ''
  echo 'clusters folder contains family and id information for each subject in respective clusters.'
  echo ''
  echo 'analysis folder contains extracted PCA results for each of the returned clusters.'
} > $RESULTSDIR/README.txt
