#!/bin/bash
# @Date    : 2015-12-11 13:06:32
# @Author  : 郭金涛
# @Email   : guojt-4451@163.com
# @Version : maf2bed

awk -F '\t' '$10~/SNP/{OFS="\t";print "chr"$5,$6,$11,$13,$16,$1}' CHB_somatic_hc.all.p0.05.maf | sed 1i"CHROM\tPOS\tREF\tALT\tSAMPLE\tGENE" > CHB_somatic_hc.all.p0.05.bed

awk -F '\t' '$10~/SNP/{OFS="\t";print "chr"$5,$6,".",$11,$13}' CHB_somatic_hc.all.p0.05.maf | sed 1i"#CHROM\tPOS\tID\tREF\tALT" > CHB_somatic_hc.all.p0.05.bed.vcf
