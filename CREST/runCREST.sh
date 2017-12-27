#!/bin/bash
# @Date    : 2017-07-26 10:48:41
# @Author  : Jintao Guo
# @Email   : guojt-4451@163.com



gfServer start localhost 50000 -stepSize=5 -canStop -log=blatServer.log $ref.2bit &

# get soft-clipping positions
perl extractSClip.pl -i $tumor.bam --ref_genome $ref.fa -o $output_dir
perl extractSClip.pl -i $normal.bam --ref_genome $ref.fa -o $output_dir

# remove germline events
perl countDiff.pl -d $tumor.bam.cover -g $normal.bam.cover

# detect SV
perl CREST.pl -f $somatic.cover -d $tumor.bam -g $normal.bam --ref_genome $ref.fa \
-t $ref.2bit -o $output_dir --blatserver localhost --blatport 50000

# visualization
perl bam2html.pl -d $tumor.bam -g $normal.bam -i $somatic.predSV --ref_genome $ref.fa -o $tumor.predSV.html


gfServer stop localhost 50000