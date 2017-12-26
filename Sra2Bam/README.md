# SraAnalysis
## sra2fastq
```
usage: sra2fq
optional arguments:
  -h, --help            show this help message and exit
  -i File               input file
  -I Files              list of input files
  -r Regular_expression input files
  -p INT                gzip multiple samples simultaneously [1]

Dependencies:
  --fastq_dump DIR      PATH to picard

Output arguments:
  -o DIR                output dir or output file [/home/g]

Clustered:
  --qsub                run crest in cluster [False]
  --nodes STR           name of nodes (e.g: n1,n2,...)
  -n INT                number of nodes [1]
```
## sra2bam
```
usage: sra2bam.py
optional arguments:
  -h, --help            show this help message and exit
  -i File               input file
  -I Files              list of input files
  -f Regular_expression input files
  -p INT                gzip multiple samples simultaneously [1]
  -t INT                number of threads to allocate to each sample [1]

Required arguments:
  --fastq_dump DIR      PATH to fastq_dump
  -r FILE               faidx indexed reference sequence file

Output arguments:
  -o DIR                output dir or output file [/home/g]

Clustered:
  --qsub                run crest in cluster [False]
  --nodes STR           name of nodes (e.g: n1,n2,...)
  -n INT                number of nodes [1]
```
