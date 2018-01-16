########### call CNV ##############################
for i in `seq 1 22` X Y;
do
    python varscan.copynumber.2.2.py -r /share/data4/Genome/hg19/hg19.chr.fa -l chr$i -i /home/jintao/varsan_CNV/ZH.txt --min-coverage 20 -o /home/jintao/varsan_CNV/CNV2/ZH --mpileup_option B,q1 /share/apps/VarScan/VarScan.v2.3.9.jar --qsub
done

python varscan.copynumber.2.2.py -r /share/data4/Genome/hg19/hg19.chr.fa -l chr17 -i /home/jintao/varsan_CNV/ESCC-D1.txt --min-coverage 20 -o /home/jintao/varsan_CNV/CNV2 --mpileup_option B,q1 /share/apps/VarScan/VarScan.v2.3.9.jar --qsub

/home/jintao/samtools-0.1.18/bin/samtools mpileup -B -q1 -r chr17 -f /share/data4/Genome/hg19/hg19.chr.fa /share/data2/ESCC/SRP_PTMD/ESCC-D1N.dedup.bam /share/data2/ESCC/SRP_PTMD/ESCC-D1T.dedup.bam | java -jar /share/apps/VarScan/VarScan.v2.3.9.jar copynumber - /home/jintao/varsan_CNV/CNV2/ESCC-D1.chr17 --mpileup 1 --min-coverage 20 --min-segment-size 100

# python ~/qsub.py
#### merge copynumber#################
for i in `l *chr1.copynumber`;
do
    head -1 $i > ../merge/${i%%.*}.copynumber
    for j in `seq 1 22` X;
    do
        sed 1d ${i%%.*}.chr$j.copynumber >> ../merge/${i%%.*}.copynumber
    done
done

for i in *.copynumber;
do
    java -jar /share/apps/VarScan/VarScan.v2.3.9.jar copyCaller $i --output-file ${i%%.*}.copynumber.called
done

######## 运行copyCaller 重新定位##############
for i in *.called;
do
    perl ~/recenter_varscan_cn.pl $i
    mv $i.[rc]* recentered
    mv ${i%.*}* raw
done

for i in *.called;
do
cent=`cat $i.centerinfo`
    java -jar /share/apps/VarScan/VarScan.v2.3.9.jar copyCaller $i --output-file ${i%%.*}.copynumber.called.recentered --recenter-up $cent
done

###### R DNAcopy 分析CNV ##################
for i in *.copynumber.called.recentered;
do
    Rscript varscan2_copynumber.R $i ${i%%*}
done

####### 合并 segments.p_value, 注释 CNV
for i in *segments.p_value;
do
    perl ~/varsan_CNV/mergeSegments.pl $i --ref-arm-sizes ~/varsan_CNV/arm_sizes.txt --amp-threshold 0.9 --del-threshold -0.9 --output-basename ${i%%.*}
    sed 1d ${i%%.*}.events.tsv|grep -wv neutral|tee| intersectBed -a - -b ~/varsan_CNV/hg19.RefSeq.bed -wa -wb | sed 1i"`head -1 ${i%%.*}.events.tsv`\tchrom\ttxStart\ttxEnd\tgene" - | grep -v chrY >${i%%.*}.events.anno
done

#####合并所有样本
for i in *.events.anno;
do
    awk 'NR>1{gsub(".events.anno","",FILENAME);OFS="\t";print FILENAME,$0}' $i >> allsample.events.anno.txt
done


# for i in *.events.tsv;do echo $i `grep -v neutral $i|sed 1d|cut -f 8|sort|uniq -c|awk '{print $1}'|paste -s` >> cnv.number.txt;done
# grep -v neutral *.events.tsv|grep focal |sed 1i"sample\t`head -1 12.events.tsv`"|sed 's/.events.tsv:/\t/g' > all.focal.tsv
# awk '$10<1000{print}' all.focal.tsv|cut -f 1,9|sort |uniq -c >cnv.number.1000.txt



####### filter gene to hotnet2

for i in *.events.anno;
do
    awk 'NR>1{gsub(".events.anno","",FILENAME);OFS="\t";print FILENAME,$0}' >> ../allsample.events.anno
done

sed 1i"sample\t`head -1 merge/12.events.anno`" allsample.events.anno

# Rscript
all_event <- read.table("allsample.CNVs.events.anno", header=T, stringsAsFactors=F, sep="\t")

event <- unique(all_event[, c("sample", "event_type", "gene")])

event$gene1 <- apply(event, 1, function(x){
    if(x[2]=="amplification"){
        return(paste0(x[3],"(A)"))
        }else{
            return(paste0(x[3], "(D)"))
        }
    })

mat <- as.data.frame.matrix(table(event[, c("sample", "gene")]))

count <- apply(mat, 2, function(x){
    count=0
    for (i in 1:length(x)) {
        if(x[i]==0){
            count = count + 1
        }
    }
    return(count)
    })

count <- 302-count
hotnet2 <- sapply(unique(event$sample), function(x){
    m <- event[which(event$sample==x & event$gene %in% names(which(count>100))), ]
    return(paste(sort(m$gene1), collapse="\t"))
    })

write.table(hotnet2, "cna.100.txt", col.names=F, quote=F, sep="\t")
################

library(parallel)
all_event <- read.table("allsample.CNVs.events.anno", header=T, stringsAsFactors=F, sep="\t")

a <- paste(all_event$sample, all_event$chrom, all_event$chr_start, all_event$chr_stop, sep=":")
u <- unique(paste(all_event$sample, all_event$chrom, all_event$chr_start, all_event$chr_stop, sep=":"))

gene <- mclapply(u, function(x){
    return(paste(sort(unique(all_event[which(a==x),"gene"])), collapse=","))
   }, mc.cores=48)

uu <- cbind(u, as.data.frame(unlist(gene)))
all_event$name <- a
all_event_u <- unique(all_event[, c(1:14,19)])
all_event_u$gene <- uu[match(all_event_u$name, uu$u),2]
write.table(all_event_u,"all.txt.1",col.names=T, row.names=F, quote=F, sep="\t")



####### GISTIC ################
## type 1
for i in *.markers;
do
    sed 1d $i >> markers
    sed 1d ${i%%.*}.segmentation.value >> segmentation
done

for i in *.segmentation.value;
do
    sed 1d $i >> segmentation
done

sort -u markers|sort -k1V -k2V -k3n |sed 's/chr//g'> markers.nochr.txt
awk '{$6=$6-1;OFS="\t";gsub("chr", "", $2);print $0}' segmentation> segmentation.nochr.txt


###

for i in *.events.tsv;
do
    awk 'NR>1{OFS="\t";print $1":"$2,$1,$2"\n"$1":"$3,$1,$3}' $i >>markers.event
    awk 'NR>1{OFS="\t";split(FILENAME, a, ".");print a[1],$1,$2,$3,$6,$4}' $i >>segmentationfile.event
done


###########################  add ####################
for i in *.segmentation.value;do Rscript ~/varsan_CNV/addSeg.R $i ${i%%.*};done

for i in *.events.tsv;
do
    awk 'NR>1{OFS="\t";split(FILENAME, a, ".");print a[1],$1,$2,$3,$6,$4}' $i >>segmentationfile.event
    awk '{gsub("X", "", $1);gsub("\\.", "-",$1);OFS="\t";print $0}' ${i%%.*}.left.seg >>segmentationfile.event
done
sort -k1V -k2V -k3n -k4n segmentationfile.event >  ../GISTIC/segmentationfile.302.event.add.txt
awk '{OFS="\t";print $2":"$3,$2,$3"\n"$2":"$4,$2,$4}' segmentationfile.event | sort -u |sort -k1V -k2V -k3n >  ../GISTIC/markers.302.event.add.txt

#########################################################

# type 2
cut -f 2-3 segmentationfile.txt > markers.txt
cut -f 2,4 segmentationfile.txt >> markers.txt
sort -k1V markers.txt >markers.txt.1
mv markers.txt.1 markers.txt
nl markers.txt > markers.num.txt





#########

for i in *.events.tsv;
do
    awk 'NR>1{OFS="\t";split(FILENAME, a, ".");print a[1],$1,$2,$3,$6,$4}' $i | sed 's/chr//g' >> segmentedFile
done

sed 1i"barcode\tchromosome\tstart\tstop\tnum.mark\tseg.mean" segmentedFile >../JISTIC/segmentationfile.event.txt

awk 'NR>1{print $2":"$3"\t"$2"\t"$3"\n"$2":"$4"\t"$2"\t"$4}' segmentationfile.event.txt | sort -u | sort -k1V -k2V -k3n | sed 1i"Probe\tChrom\tBasePair" >markers.event.txt



java -Xmx200G -classpath ~/JISTIC/JISTIC.jar JISTIC.convertSEG segmentationfile.event.txt markers.event.txt excludedregions=/home/jintao/JISTIC/glioexample/CNV.XY.txt IncludeRemovedMarkers verbose > MatrixFile

java -Xmx1500m -classpath ~/JISTIC/JISTIC.jar JISTIC.filterMarkers MatrixFile
 > FilteredMatrixFile

ffkitnjava -Xmx250G -classpath ~/JISTIC/JISTIC.jar JISTIC.Distribution spec=focal/GISTICFocal.spec copynumber=matrixFile locations=hg19_Gene_Info.txt bands=cytoBand.txt

java -Xmx250G -classpath ~/JISTIC/JISTIC.jar JISTIC.Distribution spec=limited/GISTICLimited.spec copynumber=MatrixFile locations=hg19_Gene_Info.txt bands=cytoBand.txt

java -Xmx1500m -classpath ~/JISTIC/JISTIC.jar JISTIC.Convert2IGV output






python hotnet2.3.py --hotnet2 /share/apps/hotnet2-1.0.0 --snv_file /home/jintao/hotnet2/new/mutation/p0.05_302_snv_indel_cna_100_noY_0.5/snv.p0.05.txt -o /home/jintao/hotnet2/new/mutation/p0.05_302_snv_indel_cna_100_noY_0.5 -p 48 --cna_file /home/jintao/hotnet2/new/mutation/p0.05_302_snv_indel_cna_100_noY_0.5/cna.100.txt

