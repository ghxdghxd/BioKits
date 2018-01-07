# oncotator
## Install oncotator
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
```shell
oncotator.py [-h] [-v] -i file --input-format Format
            [--output-format Format] [--qsub node]
```
1 输入为VCF格式

2 输入为MAFLITE格式

|NCBI_Build|sample|chr|start|end|ref_allele|alt_allele|
|:--------:|:----:|:-:|:---:|:-:|:--------:|:--------:|
| hg19 | 7 | 55259515 | 55259515 | T | G
| hg19 | 7 | 140453136 | 140453136 | A | T
| hg19 | 1 | 120612003 | 120612004 | GG | -
| hg19 | 8 | 145138175 | 145138176 | - | G
| hg19 | 12 | 42538391 | 42538401 | GGAGCGAGCAG | -


## Oncotator Error解决方法记录
### 1. 注释 varscan VCF出错
#### **由于VCF格式不正确(oncotator只能识别标准VCF格式)**
|   |标准VCF|varScan VCF|
|:-:|:-:|:--------:|
|##头部注释|无RD|有RD|
|##头部注释|##FORMAT=<ID=AD,Number=**R**|##FORMAT=<ID=AD,Number=**1**|
|FORMAT: AD|复合形式 **:RD,AD:**|分别显示 **RD**:**AD**|

#### **解决方法**
1. 修改头部注释
\#\#FORMAT=<ID=AD,Number=1
改为
\#\#FORMAT=<ID=AD,Number=R,
否则出现 ValueError: invalid literal for float(): 932,1
2. 修改FORMAT AD的部分
    改为与标准VCF格式一致，AD包括ref_depth,alt_depth
    代码如下
```shell
awk -F '\t' '$0 !~ /^##FORMAT=<ID=RD/ {OFS="\t";
if($0 ~ /^##FORMAT=<ID=AD/){gsub("Number=1", "Number=R",$0); print $0; next};
if($0 ~ /^#/){print $0; next};
split($9, a, ":"); split($10, b, ":");
for(i in a){if(a[i] == "RD"){ind = i}}; gsub("RD:", "", $9);
b[ind + 1] = b[ind]","b[ind + 1]; delete b[ind]; c = b[1];
for(i = 2; i <= length(b); i++){if(i != ind){c = c":"b[i]}};
$10 = c; print $0}' varScan.vcf > NEW.varScan.vcf
```
