geneToRegion <- read.table("~/Projects/ESCC/varscan_CNV/CNV2/merge/0.5/JISTIC/limited/geneToRegion.txt", 
                           header=F, sep="\t",stringsAsFactors = F)
colnames(geneToRegion) <- c("gene", "pos")

peak <- read.table("~/Projects/ESCC/varscan_CNV/CNV2/merge/0.5/JISTIC/limited/peaks.txt", 
                   header=T, sep="\t",stringsAsFactors = F)
amp_matrix <- read.table("~/Projects/ESCC/varscan_CNV/CNV2/merge/0.5/JISTIC/limited/AMP.tumors.matrix", 
                         header=T, sep="\t",stringsAsFactors = F)
amp_continuous_matrix <- read.table("~/Projects/ESCC/varscan_CNV/CNV2/merge/0.5/JISTIC/limited/AMP.tumors.continuous.matrix", 
                         header=T, sep="\t",stringsAsFactors = F)

del_matrix <- read.table("~/Projects/ESCC/varscan_CNV/CNV2/merge/0.5/JISTIC/limited/DEL.tumors.matrix", 
                         header=T, sep="\t",stringsAsFactors = F)
del_continuous_matrix <- read.table("~/Projects/ESCC/varscan_CNV/CNV2/merge/0.5/JISTIC/limited/DEL.tumors.continuous.matrix", 
                                    header=T, sep="\t",stringsAsFactors = F)
# all_lesions.conf_99.txt
## amp
peak_amp <- as.data.frame(unique(geneToRegion[grep("PEAK:AMP", geneToRegion$pos), 2]))
colnames(peak_amp) <- "pos"
peak_amp$pos <- as.character(peak_amp$pos)
peak_amp$`Unique Name` <- sapply(1:nrow(peak_amp), function(x){
  return(paste("Amplification Peak", x))
})
peak_amp$Descriptor <- sapply(peak_amp$pos, function(x){
  a=strsplit(x, "\\(")
  b=strsplit(unlist(a)[2], "\\)")
  return(unlist(b)[1])
})
peak_amp$`Wide Peak Limits` <- sapply(peak_amp$pos, function(x){
  a=strsplit(x, "chr")
  b=strsplit(unlist(a)[2], "\\(")
  return(paste0("chr",unlist(b)[1]))
})
peak_amp$`Peak Limits` = peak_amp$`Wide Peak Limits`
peak_amp$`Region Limits` = peak_amp$`Wide Peak Limits`

peak_amp$`q values` = peak$qvalue[match(peak_amp$pos, peak$peak)]
peak_amp$`Residual q values after removing segments shared with higher peaks` = peak_amp$`q values`
peak_amp$`Broad or Focal` = ""
peak_amp$`Amplitude Threshold` = "0: t<0.25; 1: 0.25<t< 0.5; 2: t>0.5"
peak_amp <- cbind(peak_amp, amp_matrix[match(peak_amp$pos, amp_matrix$Region), -1])

peak_continuous_amp <- peak_amp[,1:10]
peak_continuous_amp$`Unique Name` <- paste(peak_continuous_amp$`Unique Name`, "- CN values")
peak_continuous_amp$`Amplitude Threshold` <- "Actual Copy Change Given"
peak_continuous_amp <- cbind(peak_continuous_amp, amp_continuous_matrix[match(peak_amp$pos, amp_continuous_matrix$Region), -1])

## del
peak_del <- as.data.frame(unique(geneToRegion[grep("PEAK:DEL", geneToRegion$pos), 2]))
colnames(peak_del) <- "pos"
peak_del$pos <- as.character(peak_del$pos)
peak_del$`Unique Name` <- sapply(1:nrow(peak_del), function(x){
  return(paste("Deletion Peak", x))
})
peak_del$Descriptor <- sapply(peak_del$pos, function(x){
  a=strsplit(x, "\\(");b=strsplit(unlist(a)[2], "\\)");
  return(unlist(b)[1])
})
peak_del$`Wide Peak Limits` <- sapply(peak_del$pos, function(x){
  a=strsplit(x, "chr");b=strsplit(unlist(a)[2], "\\(");
  return(paste0("chr",unlist(b)[1]))
})
peak_del$`Peak Limits` = peak_del$`Wide Peak Limits`
peak_del$`Region Limits` = peak_del$`Wide Peak Limits`

peak_del$`q values` = peak$qvalue[match(peak_del$pos, peak$peak)]
peak_del$`Residual q values after removing segments shared with higher peaks` = peak_del$`q values`
peak_del$`Broad or Focal` = ""
peak_del$`Amplitude Threshold` = "0: t>-0.25; 1: -0.25>t> -0.5; 2: t<-0.5"
peak_del <- cbind(peak_del, del_matrix[match(peak_del$pos, del_matrix$Region), -1])

peak_continuous_del <- peak_del[,1:10]
peak_continuous_del$`Unique Name` <- paste(peak_continuous_del$`Unique Name`, "- CN values")
peak_continuous_del$`Amplitude Threshold` <- "Actual Copy Change Given"
peak_continuous_del <- cbind(peak_continuous_del, del_continuous_matrix[match(peak_del$pos, del_continuous_matrix$Region), -1])

write.table(rbind(peak_amp, peak_del, peak_continuous_amp, peak_continuous_del)[,-1], 
  "all_lesions.conf_99.txt", col.names = T, row.names = F, sep="\t", quote = F)

# amp_genes.conf_99.txt
## cytoband
## q value
## residual q value
## wide peak boundaries
## genes in wide peak

amp_gene <- peak_amp[,c(1,3,7,8,5)]
colnames(amp_gene) <- c("pos","cytoband", "q value", "residual q value", "wide peak boundaries")
head(amp_gene)
amp_gene$`genes in wide peak` <- sapply(1:nrow(amp_gene), function(i){
  chrom = unlist(strsplit(amp_gene$cytoband[i],"[pq]"))[1]
  a <- read.table(paste0("~/service/hpc/Projects/ESCC/varscan_CNV/CNV2/merge/0.5/JISTIC/limited/AMP.genes.chr",chrom,".matrix"),
                  header=T, sep="\t", stringsAsFactors = F)
  b <- a[,match(gsub("-", ".", gsub("[\\:\\-\\(\\)]", ".", amp_gene$pos[i])), colnames(a))]
  return(paste(a$X[which(b>0)], collapse = "\t"))
})

write.table(rbind(t(amp_gene[,2:5]), t(str_split(amp_gene$`genes in wide peak`, pattern="\t", simplify = T))),
 "amp_genes.conf_99.txt", sep="\t", col.names = F, row.names = T, quote = F)

del_gene <- peak_del[,c(1,3,7,8,5)]
colnames(del_gene) <- c("pos","cytoband", "q value", "residual q value", "wide peak boundaries")
del_gene$`genes in wide peak` <- sapply(1:nrow(del_gene), function(i){
  chrom = unlist(strsplit(del_gene$cytoband[i],"[pq]"))[1]
  a <- read.table(paste0("~/service/hpc/Projects/ESCC/varscan_CNV/CNV2/merge/0.5/JISTIC/limited/DEL.genes.chr",chrom,".matrix"),
                  header=T, sep="\t", stringsAsFactors = F)
  b <- a[,match(gsub("-", ".", gsub("[\\:\\-\\(\\)]", ".", del_gene$pos[i])), colnames(a))]
  return(paste(a$X[which(b>0)], collapse = "\t"))
})
write.table(rbind(t(del_gene[,2:5]), t(str_split(del_gene$`genes in wide peak`, pattern="\t", simplify = T))),
  "del_genes.conf_99.txt", sep="\t", col.names = F, row.names = T, quote = F)



### for hotnet2
# gene <- read.table("gene.matrix", header=T, sep="\t", stringsAsFactors=F)
# geneToregion <- read.table('geneToRegion.txt', header=F, sep="\t", stringsAsFactors=F)
# gene <- gene[which(gene$Gene %in% geneToregion$V1[grep("PEAK", geneToregion$V2)]), ]
# a <- apply(gene[,2:303], c(1,2), function(x){if(x!="N"){return(paste0("(",x,")"))}else(return(""))})
# for(i in 1:302){a[,i] <- paste0(gene[,1], a[,i])}
# b <- apply(a, c(1,2), function(x){if(length(grep("\\(", x))>0){return(x)}else{return("")}})
# b<-t(b)
# c <- apply(b, 1, function(x){x=unique(unlist(x));paste(x, collapse="\t")})
# names(c) <- gsub("\\.", "-",gsub("X","",names(c)))
# write.table(c, "gene.matrix.2", col.names=F,row.names=T, quote=F, sep="\t")


laml.gistic = readGistic(gisticAllLesionsFile = "all_lesions.conf_99.txt", 
  gisticAmpGenesFile = "amp_genes.conf_99.txt", 
  gisticDelGenesFile = "del_genes.conf_99.txt",
  gisticScoresFile = "gscore.gp_gistic.gistic.txt")

gcp = gisticChromPlot(gistic = laml.gistic, markBands = "all", markBandsCol = "black",
  width=15, height=5, file = "scnv")
gbp = gisticBubblePlot(gistic = laml.gistic)


laml.gistic.p = readGistic(gisticAllLesionsFile = "all_lesions.conf_99.txt", 
  gisticAmpGenesFile = "amp_genes.conf_99.txt", 
  gisticDelGenesFile = "del_genes.conf_99.txt",
  gisticScoresFile = "pval.gp_gistic.gistic.txt")

gcp = gisticChromPlot(gistic = laml.gistic.p, markBands = "all", markBandsCol = "black",
  fdrCutOff = 0.001, width=15, height=5, file = "scnv.p")


laml.gistic.peak = readGistic(gisticAllLesionsFile = "all_lesions.conf_99.txt", 
  gisticAmpGenesFile = "amp_genes.conf_99.txt", 
  gisticDelGenesFile = "del_genes.conf_99.txt",
  gisticScoresFile = "peak.gp_gistic.gistic.txt")

gcp = gisticChromPlot(gistic = laml.gistic.peak, markBands = "all", markBandsCol = "black",
  width=15, height=5, file = "scnv.peak")