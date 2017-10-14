#!/bin/bash
# @Date    : 2015-09-21 10:54:26
# @Author  : 郭金涛
# @Email   : guojt-4451@163.com
# @Version : 1.0

for i in *.maf;
do
    name=${i%%.*}T;
    grep -Ev "#|Hugo_Symbol|^Unknown" $i |awk -F '\t' '$14!~/rs/{
        if($202=="True"){
            OFS="\t";
            gsub(/__UNKNOWN__/,"NA",$0);
            $3="XMU";
            $4=$203;
            $8="+";
            $16=$33=$17="'$name'";
            print $0}}' >> ../tumor.hc.maf
done

#sed 1,4d $i | grep -Ev "Hugo_Symbol|^Unknown" |awk -F '\t' '{

for i in *.maf;
do
    name=${i%%.*}T;
    sed 1,4d $i | grep -Ev "Hugo_Symbol" |awk -F '\t' '{
        if($202=="True"){
            OFS="\t";
            gsub(/__UNKNOWN__/,"NA",$0);
            $3="XMU";
            $4=$203;
            $8="+";
            $16=$33=$17="'$name'";
            print $0}}' >> ../tumor.hc.all.maf
done

for i in *.vcf.maf;
do
    name=${i%%.*};
    sed 1,4d $i | grep -Ev "Hugo_Symbol|^Unknown" |awk -F '\t' '{
        OFS="\t";
        gsub(/__UNKNOWN__/,"NA",$0);
        $16=$33=$17=$34="'$name'";
        print $0}' >> all.maf
done


## including IGR
for i in *.maf;
do
    name=${i%%.*}T;
    grep -Ev "#|Hugo_Symbol" $i |awk -F '\t' '{
        if($202=="True"){
            OFS="\t";
            gsub(/__UNKNOWN__/,"NA",$0);
            $3="XMU";
            $4=$203;
            $8="+";
            $16=$33=$17="'$name'";
            print $0}}' >> ../tumor.hc.all.maf
done


awk -F '\t' '{OFS="\t";print $1,$16,$9,$5,$6,$7,$11,$12,$13}'

python hotnet2.2.py --hotnet2 /share/data4/hotnet2-1.0.0/ \
--snv_file /home/jintao/paper_mutation/japan/hotnet2/snv.p0.05.txt \
-o /home/jintao/paper_mutation/japan/hotnet2/ -p 40


289 gencode_transcript_status

## maf Variant_Classification 9
## maf 10 SNP/INS/DEL
## maf 289 KNOWN/NOVEL/PUTATIVE


python change_variant.py  >>>>  CHB_somatic_hc.all.maf.1

## SNV
awk -F '\t' '{if ($307=="Silent" && $10=="SNP"){print $0}}' CHB_somatic_hc.all.maf.1 |cut -f 5-7|sort -u|wc



#INDEL


awk -F '\t' '{if ($10!="SNP"){print $0}}' CHB_somatic_hc.all.maf.1 |cut -f 5-7|sort -u|wc


awk -F '\t' '$35~/T>G|A>C/{if ($10=="SNP"){print $0}}' CHB_somatic_hc.all.maf |cut -f 66|awk -F '' '{OFS="";print $10,$11,$12}'| tr '[a-z]' '[A-Z]'|sort |uniq -c |sort -nr|less




##ESCC Germline, maf2
##DP4, 123
##FREQ, 160
#GATK
for i in *.maf;
do
    grep -v -e "Hugo_Symbol" -e "#" -e "^Unknown" $i |awk -F '\t' '{OFS="\t";print $1,$5,$6,$7,$11,$13}' >> gatk.maflite.txt
done




awk -F '\t' '{OFS="\t";if($10~/SNP/){split($66,a,"");b=a[9]a[10]a[11]a[12]a[13];if(b!="GGTGG"&& $13!="G"){print $0}}else{print $0}}' CHB_somatic_hc.all.2016815.IGR.maf > CHB_somatic_hc.all.2016815.IGR.rmTG.maf