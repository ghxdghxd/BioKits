export PATH=/share/apps/breakdancer/lib/breakdancer-max1.4.5-unstable-66-4e44b43:$PATH
tumorBam=/r730_iscsi/analysis/lch/data_release/alignment/HCC2218C_S3/HCC2218C_S3.bam.dedup.bam
normalBam=/r730_iscsi/analysis/lch/data_release/alignment/HCC2218BL_S4/HCC2218BL_S4.bam.dedup.bam

bam2cfg.pl -q 35 -g -h $tumorBam $normalBam > test.cfg
breakdancer-max -t -h -q 35 -y 40 -d test.ctx test.cfg > test.ctx

bam2cfg.pl -q 35 -g -h NA19238_chr21_del_inv.bam NA19240_chr21_del_inv.bam > NA19.cfg
breakdancer-max -t -h -q 35 -y 40 -d NA19 NA19.cfg > NA19


for i in `ls /home/jintao/Projects/genokon/sjzp/Illumina_B170[1-9]*.bam`;
do 
sample=`basename $i`
sample=${sample%%.*}
cat>job<< EOF
#!/bin/bash
#PBS -N $sample
#PBS -l nodes=1:compute-0-4
#PBS -j oe
source /etc/profile.d/set.sh
cd /home/jintao/Projects/genokon/sjzp/breakdancer
/share/apps/breakdancer/lib/breakdancer-maxunstable/bam2cfg.pl -g -h \
$i /home/jintao/Projects/genokon/sjzp/Illumina_B17NC.bam.dedup.bam > $sample.cfg
# breakdancer-max -h -d $sample $sample.cfg > $sample.out
EOF
qsub job
done

for i in *.cfg;
do
awk -F '\t' '{OFS="\t";name=$1;gsub("readgroup:","", name);$5="lib:"name;print $0}' $i > $i.1
mv $i.1 $i
done

for i in `ls /home/jintao/Projects/genokon/sjzp/Illumina_B170[1-9]*.bam`;
do 
sample=`basename $i`
sample=${sample%%.*}
cat>job<< EOF
#!/bin/bash
#PBS -N $sample
#PBS -l nodes=1:compute-0-4
#PBS -j oe
source /etc/profile.d/set.sh
cd /home/jintao/Projects/genokon/sjzp/breakdancer
breakdancer-max -h -d $sample $sample.cfg > $sample.out
EOF
qsub job
done







for i in *segments.p_value;
do
    awk '$2!~/[GTd]/{print}' $i|grep -v NA|tee|perl ../mergeSegments.pl - --ref-arm-sizes ../arm_sizes.nochr.txt --amp-threshold 0.5 --del-threshold -0.5 --output-basename ${i%%.*}
    sed 1d ${i%%.*}.events.tsv|grep -wv neutral|tee| intersectBed -a - -b ../hg19.RefSeq.nochr.bed -wa -wb | sed 1i"`head -1 ${i%%.*}.events.tsv`\tchrom\ttxStart\ttxEnd\tgene" - | grep -v chrY >${i%%.*}.events.anno
done