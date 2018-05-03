#!/bin/bash
chrom=$1
cancer=$2

start=(`grep -w ^$chrom impute.region.txt|cut -f 2`)
end=(`grep -w ^$chrom impute.region.txt|cut -f 3`)

for((i=0;i < ${#start[@]};i++));do
python GenoImpute.py impute -g /share/data4/TCGA/SNP6_Genotype/2_shapeit_phased/$cancer/$cancer\_CEU_chr$chrom\_phased.haps.gz \
-c chr$chrom -s ${start[$i]} -e ${end[$i]} \
-o /home/jintao/Projects/impute2/TCGA/$cancer/chr$chrom
done
