#!/usr/bin/env Rscript
# @Date    : 2017-03-6 11:40:41
# @Author  : Jintao Guo
# @Email   : guojt-4451@163.com
# @Version : $Id$


## LUAD

maf <- read.csv("LUAD-TP.final_analysis_set.maf", sep="\t", header=T, stringsAsFactors = F)
mutsig <- read.csv("sig_genes.txt", sep="\t", header=T, stringsAsFactors=F)
snv <- maf[which(maf$Hugo_Symbol %in% mutsig$gene[which(mutsig$p <0.01)]), c("Hugo_Symbol", "Tumor_Sample_Barcode")]
snv <- as.data.frame.matrix(t(table(snv)))
rownames(snv) <- sapply(rownames(snv), function(x){return(substr(x,1,12))})

snv <- apply(snv, 1, function(x){paste(colnames(snv)[which(x>0)],collapse="\t")})
write.table(snv, "snv.txt",col.names=F, row.names=T,quote=F,sep="\t")

