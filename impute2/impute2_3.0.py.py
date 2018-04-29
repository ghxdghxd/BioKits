#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# impute2_3.0.py.py
# @Author : JT Guo
# @Email  : guojt-4451@163.com
# @Date   : 2018-4-28 09:09:45


import os
import tempfile

version = 3.0


def makedir(new_dir, exist_dir=None):
    """Make a directory. If it doesn't exist, handling concurrent race conditions.
    """
    if exist_dir:
        new_dir = os.path.join(exist_dir, new_dir)
    if os.path.exists(new_dir):
        print("The " + new_dir + " is already exist")
    else:
        print("Make " + new_dir)
        os.makedirs(new_dir)


def impute(l):
    chrom = l[0]
    int_start = l[1]
    int_end = l[2]
    ref_map = "/share/data4/Genome/ALL_1000G_phase1integrated_v3_impute/1000GP_Phase3/genetic_map_chr" + \
        chrom + "_combined_b37.txt"
    ref_hap = "/share/data4/Genome/ALL_1000G_phase1integrated_v3_impute/1000GP_Phase3/1000GP_Phase3_chr" + \
        chrom + ".hap.gz"
    ref_legend = "/share/data4/Genome/ALL_1000G_phase1integrated_v3_impute/1000GP_Phase3/1000GP_Phase3_chr" + \
        chrom + ".legend.gz"
    genotype_g = "/home/jintao/impute2/PRAD/input/" + chrom + ".sorted.txt"
    strand_g = "/home/jintao/impute2/PRAD/input/" + chrom + ".strand"
    sample_g = "/home/jintao/impute2/PRAD/input/sample_g.txt"
    outfile1 = "/home/jintao/impute2/PRAD/output/" + chrom + \
        "." + int_start + "-" + int_end + ".prephasing.impute2"
    outfile2 = "/home/jintao/impute2/PRAD/output/" + chrom + \
        "." + int_start + "-" + int_end + ".prephased.impute2"
    ftmp = tempfile.NamedTemporaryFile()
    ftmp.write("#!/bin/bash\n")
    ftmp.write(
        "#PBS -o " + os.path.split(os.path.realpath(__file__))[0] + "/log\n")
    ftmp.write("#PBS -l nodes=1:ppn=2,walltime=10:00:00\n")
    ftmp.write("#PBS -j oe\ncd $PBS_O_WORKDIR\n")
    if chrom != "X":
        cmd1 = impute2 +\
            " -prephase_g " +\
            " -m " + ref_map +\
            " -g " + genotype_g +\
            " -int " + str(int_start) + " " + str(int_end) +\
            " -Ne 20000" +\
            " -o " + str(outfile1)
        cmd2 = impute2 +\
            " -use_prephased_g" +\
            " -m " + ref_map +\
            " -h " + ref_hap +\
            " -l " + ref_legend +\
            " -known_haps_g " + outfile1 + "_haps" +\
            " -strand_g " + strand_g +\
            " -int " + str(int_start) + " " + str(int_end) +\
            " -Ne 20000" +\
            " -o " + str(outfile2) +\
            "-phase"
        ftmp.write(cmd1 + "\n")
        ftmp.write(cmd2)
    else:
        cmd = impute2 +\
            " -prephase_g" +\
            " -chrX" +\
            " -m /share/data4/Genome/ALL_1000G_phase1integrated_v3_impute/genetic_map_chrX_nonPAR_combined_b37.txt" +\
            " -h /share/data4/Genome/ALL_1000G_phase1integrated_v3_impute/ALL_1000G_phase1integrated_v3_chrX_nonPAR_impute.hap.gz" +\
            " -l /share/data4/Genome/ALL_1000G_phase1integrated_v3_impute/ALL_1000G_phase1integrated_v3_chrX_nonPAR_impute.legend.gz" +\
            " -g " + genotype_g +\
            " -sample_g " + sample_g +\
            " -int " + str(int_start) + " " + str(int_end) +\
            " -Ne 20000" +\
            " -o " + str(outfile2)
        ftmp.write(cmd)
        # cmd1 = impute2 +\
        #     " -prephase_g" +\
        #     " -chrX" +\
        #     " -m /share/data4/Genome/ALL_1000G_phase1integrated_v3_impute/genetic_map_chrX_nonPAR_combined_b37.txt" +\
        #     " -g " + genotype_g +\
        #     " -sample_g " + sample_g +\
        #     " -int " + str(int_start) + " " + str(int_end) +\
        #     " -Ne 20000" +\
        #     " -o " + str(outfile1)
        # cmd2 = impute2 +\
        #     " -use_prephased_g" +\
        #     " -chrX" +\
        #     " -m /share/data4/Genome/ALL_1000G_phase1integrated_v3_impute/genetic_map_chrX_nonPAR_combined_b37.txt" +\
        #     " -h /share/data4/Genome/ALL_1000G_phase1integrated_v3_impute/ALL_1000G_phase1integrated_v3_chrX_nonPAR_impute.hap.gz" +\
        #     " -l /share/data4/Genome/ALL_1000G_phase1integrated_v3_impute/ALL_1000G_phase1integrated_v3_chrX_nonPAR_impute.legend.gz" +\
        #     " -known_haps_g " + outfile1 + "_haps" +\
        #     " -int " + str(int_start) + " " + str(int_end) +\
        #     " -Ne 20000" +\
        #     " -o " + str(outfile2) +\
        #     " -phase"
        # ftmp.write(cmd1 + "\n")
        # ftmp.write(cmd2)
    ftmp.seek(0)
    # print(ftmp.read())
    os.system("qsub " + ftmp.name)
    ftmp.close()


if __name__ == '__main__':
    impute2 = "/share/apps/impute_v2.3.2_x86_64_static/impute2"
    with open("/home/jintao/impute2/PRAD/input/ng_3094_S2.hg19.bed") as f:
        loci_list = map(lambda x: x.strip().split("\t"), f.readlines())
    for l in loci_list:
        impute(l)
