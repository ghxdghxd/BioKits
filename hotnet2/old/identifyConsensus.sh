#!/bin/bash

set -o nounset                              # Treat unset variables as an error

python identifyConsensus.py --results_files /home/jintao/hotnet2/p0.05/somatic_hc/hint+hi2012/delta_0.000149211019373/results.json /home/jintao/hotnet2/p0.05/somatic_hc/hprd/delta_0.000605680846217/results.json /home/jintao/hotnet2/p0.05/somatic_hc/iref/delta_0.000440583520307/results.json /home/jintao/hotnet2/p0.05/somatic_hc/multinet/delta_0.000431922864868/results.json --networks hint+hi2012 hprd iref multinet -o /home/jintao/hotnet2/identifyConsensus.json




python identifyConsensus.py --results_files /home/jintao/hotnet2/new/mutation/p0.05_NaN/hint+hi2012/delta_0.000141676004445/results.json /home/jintao/hotnet2/new/mutation/p0.05_NaN/hprd/delta_0.000114736125609/results.json /home/jintao/hotnet2/new/mutation/p0.05_NaN/iref/delta_0.000215263625089/results.json /home/jintao/hotnet2/new/mutation/p0.05_NaN/multinet/delta_0.000137816111905/results.json --networks hint+hi2012 hprd iref multinet -o /home/jintao/hotnet2/new/mutation/p0.05_NaN/identifyConsensus.json



python /share/apps/hotnet2-1.0.0/identifyConsensus.py --results_files ./hint+hi2012/delta_0.00014906608532/results.json ./hprd/delta_8.13760587238e-05/results.json ./iref/delta_8.41201582657e-05/results.json ./multinet/delta_0.00011346305738/results.json --networks hint+hi2012 hprd iref multinet -o identifyConsensus.json


python /share/apps/hotnet2-1.0.0/identifyConsensus.py --results_files ./hint+hi2012/delta_0.00014906608532/results.json -ms 3 --networks hint+hi2012 -o test.json


grep -v -e ":" -e "}" -e "{" -e "]" identifyConsensus.json|sed 's/"//g;s/,//g;s/ //g' > identifyConsensus.gene

