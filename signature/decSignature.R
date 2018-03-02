#!/usr/bin/env Rscript
# @Date    : 2018-03-01 12:49:13
# @Author  : Jintao Guo
# @Email   : guojt-4451@163.com

Args <- commandArgs(T)

input = Args[1]  # with header: sampel chr pos ref alt

library(deconstructSigs)

sigs.input <- mut.to.sigs.input(mut.ref = input, 
                                sample.id = "sample", 
                                chr = "chr", 
                                pos = "pos", 
                                ref = "ref", 
                                alt = "alt")

sigs = whichSignatures(tumor.ref = sigs.input, 
                        signatures.ref = signatures.cosmic, 
                        sample.id = rownames(sigs.input), 
                        contexts.needed = TRUE,
                        tri.counts.method = 'default')

plotSig <- function(sigOut) {
    colnames(sigOut$weights) <- gsub("Signature.", "Sig", colnames(sigOut$weights))
    plotSignatures(sigOut)
}

pdf(paste0(rownames(sigs.input), ".decSig.pdf"))
plotSig(sigs)
makePie(sigs)
dev.off()

