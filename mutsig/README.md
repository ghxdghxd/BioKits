# MutSigCV
## 嗅觉基因，OR, 显著的原因
肿瘤异质性的突变过程， 一些基因，如ORs，虽然对肿瘤没有功能影响，
但因为能够迅速积累突变，而被鉴定出来。

A tab-delimited report of significant mutations,
listed in descending order from most significant to least significant.
The "nnei","x", and "X" values in the MutSig output analysis give insight
into how the background mutation rate is calculated for a given gene.
nnei gives the number of neighboring genes that are pooled together to
compute the background mutation rate for that gene;
these genes are not necessarily adjacent on the genome, but rather they have nearby covariate values.
x gives the number of mutated bases in these neighboring genes that are either silent or non-coding,
while X gives the total number of bases related to these neighboring genes.
