# [hotnet2](https://github.com/raphael-group/hotnet2)
hotnet2是一种新型计算机算法，能够筛选庞大的遗传数据，发现相互作用基因，而这些基因一旦突变就会导致多种癌症的发生发展。
## 背景
预先定义一个合理的基因集合或组合的统计需求
标准化分析突变相关通路及蛋白复合物的方法
这项研究没有选择癌症遗传学研究，而是采用了与众不同的方式，寻找癌症样品中频繁出现的单个基因的突变。基因并不会常常独自作战，主要还是与其它基因形成网络和途径，调控细胞功能。在某些情况下，途径中多个基因中出现一个突变，就会引发故障发生，导致癌症。因为有害突变可以分布在多个这样的基因网络，因此难以通过统计检验发现它们。而Hotnet2运算方法则能在网络水平上分析基因，帮助科学家们识别罕见，但又重要的癌症突变。
HotNet2 算法是通过将病患突变数据投射到一张基因相互作用图谱上，然后寻找比偶发突变更常见的突变之间的相互作用网络，这一程序能作为heat
sources寻找经常突变的基因。通过分析图谱上分布和聚集的方式，HotNet2 能找到与癌症相关的的“热”网络。
Hotnet的最初版本已经被用于识别急性髓细胞白血病、卵巢癌和几个其它类型癌症中的重要网络，目前这一版本也经过修改，可以用于处理更大和更复杂的泛癌症数据集。

## subnetworks.json结构
+ stats : {Minimum edge weight δ: { Minimum subnetwork size1 : { expected :XXX,  observed : XXX,  pval : XXX} Minimum subnetwork size2 : { expected : XXX,  observed : XXX,  pval : XXX}}}

+ deltas : [ deltas1, deltas2, deltas3, deltas4]        Minimum edge weight δ

+ mutation_matrices : {deltas 1:{}, deltas2 : {}, deltas3 : {}, deltas4 : {}}    基因不同位点的类型（snv indel 或 cnv)

+ typeToSamples

+ subnetworks：{deltas :[ network0, network1] } network: {edges :{

+ ks Minimum subnetwork size

+ sampleToTypes： 样本类型

## runHotnet2.py:
## with or without mutsig
```shell
runhotnet2.py
--hotnet2 $PATH2hotnet2 \
--maf_file  $maf_file \
[ --mutsig_file $mutsig_file ] \
[ --threshold q0.01 ] \
[ -o $output_dir ]\
[ -p $thread_num ] \
[ --cna_file $Scna_file ]
```

## visnetwork.R:
```shell
visnetwork.R mat.txt heat.txt
```
