#!/bin/bash
# @Date    : 2015-09-28 13:39:14
# @Author  : 郭金涛
# @Email   : guojt-4451@163.com
# @Version : $Id$

for i in *vcf;do /share/apps/annovar/convert2annovar.pl --format vcf4 -allsample -withfreq $i > ../../annovar/CHB_somatic_hc/${i%.*}.annova;done

for i in *annova
./annotate_variation.pl --buildver hg19 --geneanno --outfile /home/jintao/output/ESCC/oncotator/12.anno /home/jintao/output/ESCC/oncotator/12.annovar humandb


annotate_variation.pl -downdb -webfrom annovar -buildver hg19 dbnsfp30a humandb/
table_annovar.pl ex1.avinput humandb/ -protocol dbnsfp30a -operation f -build hg19 -nastring .