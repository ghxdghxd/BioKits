args<-commandArgs(TRUE)
file_called <- args[1]
file_name <- args[2]


library(DNAcopy)

# files <- list.files(pattern="*.copynumber.called")
# cn <- read.table(files[1], header = T)
# cn <- cn[, c("chrom", "chr_start", "adjusted_log_ratio")]
# colnames(cn)[3] <- as.character(strsplit(files[1], ".copynumber.called"))

# for(i in files[2:length(files)]){
# 	cn_tmp <- read.table(i, header = T)
# 	cn_tmp <- cn_tmp[, c("chrom", "chr_start", "adjusted_log_ratio")]
# 	colnames(cn_tmp)[3] <- as.character(strsplit(i, ".copynumber.called"))
# 	cn <- merge(cn, cn_tmp, all=T)
# }


cn <- read.table(file_called, header = T)
CNA.object <- CNA( genomdat = cn$adjusted_log_ratio, chrom = cn$chrom, maploc = cn$chr_start, data.type= 'logratio', sampleid = file_name)
CNA.smoothed <- smooth.CNA(CNA.object)
segment <- segment(CNA.smoothed, verbose=1, min.width=2, undo.SD=3)
p.segment <- segments.p(segment)
pdf(paste0(file_name, ".pdf"))
plot(segment, plot.type="w")
plot(segment, plot.type="s")
plot(segment, plot.type="p")
plot(segment, plot.type="c")
dev.off()
write.table(p.segment, file=paste0(file_name, ".segments.p_value"), quote = F, col.names=T, row.names=F, sep="\t")
