#!/bin/bash
# @Date    : 2016-03-30 12:32:26
# @Author  : 郭金涛
# @Email   : guojt-4451@163.com
# @Version : $Id$

for i in *.predSV.txt;
do
    awk '{split(FILENAME,a, "."); OFS="\t";print a[1],$0}' $i >>../escc.predSV.txt
done

