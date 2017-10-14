#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2015-11-08 21:56:31
# @Author  : 郭金涛
# @Email   : guojt-4451@163.com
# @Version : $Id$

import sys


def main():
    with open("/home/g/apps/Projects/00ESCC/result/somatic_all/CHB_somatic_hc.all.maf") as f:
        maf = f.readlines()
    a = ["variant1"]
    b = ["variant2"]
    for line in maf[1:]:
        line = line.split("\t")
        if line[8] in ["Intron", "5'UTR", "3'UTR", "RNA", "5'Flank",
                       "lincRNA"]:
            a.append("NA")
        elif line[8] == "Silent":
            a.append("Silent")
        else:
            a.append("non_silent")
        if line[158] == "frameshift_variant":
            b.append("Frame")
        elif line[158] in ["inframe_deletion", "inframe_insertion"]:
            b.append("Inframe")
        elif line[158] == "missense":
            b.append("Missense")
        elif line[158] == "stop_lost":
            b.append("Lost")
        elif line[158] == "stop_gained":
            b.append("Gain")
        elif line[158] == "synonymous_variant":
            b.append("Synonymous")
        elif line[158] == "splice_region_variant":
            b.append("Splice")
        else:
            b.append("NA")

    with open("variant1.txt", "w") as f:
        f.write("\n".join(a))
    with open("variant2.txt", "w") as f:
        f.write("\n".join(b))

if __name__ == '__main__':
    sys.exit(main())
