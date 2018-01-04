# Picard

## java.io.IOException: No space left on device
```java
java -jar picard.jar TMP_DIR=PATH_TO_TMP....
```
## mark duplicates
```
mark_duplicates.py

optional arguments:
  -h, --help            show this help message and exit
  --path2picard DIR     PATH to picard
  -i FILE               input bamfile
  -I FILE               list of input bamfile
  -o DIR                output dir or output file [/home/g]
  -p INT                analyze multiple samples simultaneously [1]
  -t INT                number of threads to allocate to each sample [1]
  -m STR                memory [2g]
  --qsub                run crest in cluster
  --remove_duplicates   remove PCR duplicates
  --validation_stringency {STRICT, LENIENT, SILENT}
                        Validation stringency for all SAM files read by this
                        program [STRICT]
  --node STR            name of nodes
  -n INT                number of nodes [1]
```
