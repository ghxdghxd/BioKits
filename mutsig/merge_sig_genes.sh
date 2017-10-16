#!/bin/bash

set -o nounset                              # Treat unset variables as an error

paste -d '\t' *.mutsig.sig_genes.txt | sed 1d |awk '{
    a=b=c=0;
    for(i=0;i<113;i++)
    {
        a+=$(8+15*i);
        b+=$(9+15*i);
        c+=$(10+15*i);
    }
    OFS="\t";
    print $1,$2,$3,$4,a,b,c
}'
