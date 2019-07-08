#!/bin/bash
# @Date    : 2015-09-28 13:39:14
# @Author  : 郭金涛
# @Email   : guojt-4451@163.com
# @Version : $Id$

for i in *vcf;do /share/apps/annovar/convert2annovar.pl --format vcf4 -allsample -withfreq $i > ../../annovar/CHB_somatic_hc/${i%.*}.annova;done

for i in *annova
./annotate_variation.pl --buildver hg19 --geneanno --outfile /home/jintao/output/ESCC/oncotator/12.anno /home/jintao/output/ESCC/oncotator/12.annovar humandb


annotate_variation.pl -downdb -webfrom annovar -buildver hg19 gnomad211_genome .

/share/apps/annovar_2018Apr16/table_annovar.pl /share/data4/TCGA/Germline/annovar/test.vcf \
    -vcfinput /share/data0/reference/annovar_humandb/hg19 \
    -build hg19 -out test -otherinfo -remove -nastring . \
    -protocol refGene,avsnp150,1000g2015aug_eur,gnomad211_genome,exac03 -operation g,f,f,f,f

/share/apps/annovar_2018Apr16/prepare_annovar_user.pl -dbtype cosmic /share/data0/reference/COSMIC/grch38/CosmicMutantExport.tsv.gz \
    -vcf /share/data0/reference/COSMIC/grch38/CosmicCodingMuts.vcf.gz > hg38_cosmic88.txt 


for i in *.vcf;
do
{perl /share/apps/annovar_2018Apr16/table_annovar.pl $i -vcfinput /share/data0/reference/annovar_humandb \
    -buildver hg38 -out ${i%.*} -otherinfo -remove \
    -protocol refGene,knownGene,ensGene,avsnp150,ljb26_all,dbnsfp35a,intervar_20180118,exac03,clinvar_20190305,regsnpintron \
    -operation g,g,g,f,f,f,f,f,f,f -nastring . &}
done

SIFT_score
SIFT_pred
MetaSVM_score
MetaSVM_pred
CADD_raw
CADD_phred
ExAC_EAS
CLNALLELEID
CLNDN
CLNDISDB
CLNREVSTAT
CLNSIG