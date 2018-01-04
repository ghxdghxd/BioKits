# oncotator
## Install
### Anaconda3 python3
```python
conda create -n oncotator python=2.7
source activate oncotator
pip install numpy
pip install pyvcf
git clone https://github.com/broadinstitute/oncotator.git
cd oncotator
python setup.py install
pip install biopython or pandas .....
```

## run oncotator
```
usage: oncotator.py [-h] [-v] -i file --input-format Format
                    [--output-format Format] [--qsub node]
```
1 输入为VCF格式

2 输入为MAFLITE格式

NCBI_Build	sample	chr	start	end	ref_allele	alt_allele
hg19	7  55259515  55259515  T  G
hg19	7  140453136  140453136 A T
hg19	1  120612003  120612004  GG  -
hg19	8  145138175  145138176  -  G
hg19	12  42538391  42538401  GGAGCGAGCAG  -
