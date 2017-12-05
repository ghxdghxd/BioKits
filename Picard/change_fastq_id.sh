#!/bin/bash

set -o nounset                              # Treat unset variables as an error

for i in *.fastq;
do
	awk '{if($0~/^[+@]SRR/){split($0,a,".");print a[1]"."a[2]"/"a[3]}else{print $0}}' $i >./new/$i
done