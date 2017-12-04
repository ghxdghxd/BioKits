#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2017-12-04 14:47:05
# @Author  : 郭金涛
# @Email   : guojt-4451@163.com

import glob
import os
import sys
import tempfile
import multiprocessing
import subprocess

version = "2.1"


def get_args():
    """Get arguments from commond line"""
    try:
        import argparse
    except ImportError as imerr:
        print("\033[1;31m" + str(imerr) + " \033[0m")
        sys.exit()

    parser = argparse.ArgumentParser(usage="%(prog)s",
                                     fromfile_prefix_chars='@',
                                     description=__doc__)

    required_group = parser.add_argument_group("Required arguments")
    required_group.add_argument("--fastq_dump",
                                metavar="DIR",
                                dest="fastq_dump",
                                help="PATH to fastq_dump")
    required_group.add_argument("-r",
                                metavar="FILE",
                                dest="ref_fa",
                                required=True,
                                help="faidx indexed reference sequence file")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-i",
                       metavar="File",
                       dest="input_file",
                       help="input file")
    group.add_argument("-I",
                       metavar="Files",
                       dest="input_list",
                       help="list of input files")
    group.add_argument("-f",
                       metavar="Regular_expression",
                       type=str,
                       dest="input_re",
                       help="input files")

    group = parser.add_argument_group("Output arguments")
    group.add_argument("-o",
                       metavar="DIR",
                       dest="output_dir",
                       default=os.environ['HOME'],
                       help="output dir or output file \
                       [" + os.environ['HOME'] + "]")
    parser.add_argument("-p",
                        metavar="INT",
                        default=1,
                        dest="processes_number",
                        type=str,
                        help="gzip multiple samples simultaneously [1]")
    parser.add_argument("-t",
                        metavar="INT",
                        default=1,
                        dest="threads_number",
                        type=str,
                        help="number of threads to allocate to \
                        each sample [1]")
    clustered_group = parser.add_argument_group("Clustered")
    clustered_group.add_argument("--qsub",
                                 action="store_true",
                                 default=False,
                                 dest="qsub",
                                 help="run crest in cluster [False]")
    clustered_group.add_argument("--nodes",
                                 metavar="STR",
                                 dest="node_name",
                                 type=str,
                                 help="name of nodes (e.g: n1,n2,...)")
    clustered_group.add_argument("-n",
                                 metavar="INT",
                                 default=1,
                                 dest="nodes_number",
                                 type=int,
                                 help="number of nodes [1]")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    else:
        check_dependencies(["bwa", "samtools"])
        parser.parse_args().output_dir = os.path.abspath(
            parser.parse_args().output_dir)
        return parser.parse_args()


def check_dependencies(tools):
    """Ensure required tools are present.
    """
    print("Checking required dependencies......")
    if not isinstance(tools, list):
        tools = [tools]
    for t in tools:
        try:
            subprocess.check_call("which " + t, stderr=subprocess.STDOUT)
        except OSError:
            raise OSError("\033[1;31m" + __file__ + " requires " +
                          t + "\033[0m")
    sys.exit()


def sra2bam(f):
    sra2fq = args.fastq_dump + " --split-3 --gzip " + f + \
        " -O " + args.output_dir
    print("Running: \n", sra2fq)
    os.system(sra2fq)
    fastq_name = os.path.basename(f).split(".sra")[0]
    fastq_r1 = os.path.join(args.output_dir, fastq_name + "_1.fastq.gz")
    fastq_r2 = os.path.join(args.output_dir, fastq_name + "_2.fastq.gz")
    bwa_mem = "bwa mem -t " + args.threads_number + \
        " -M -R \"@RG\\tID:" + fastq_name + "\\tLB:Hg19\\tPL:Illumina\\tPU:Barcode\\tSM:" + \
        fastq_name + "\\tCREATE_INDEX:True\" " + \
        args.ref_fa + " " + fastq_r1 + " " + fastq_r2 + " | samtools view -u -t " + \
        args.ref_fa + ".fai -S - -b | samtools sort - " + \
        os.path.join(args.output_dir, fastq_name + ".sorted")
    bam_index = "samtools index " + \
        os.path.join(args.output_dir, fastq_name + ".sorted.bam")
    print("Running: \n", bwa_mem)
    os.system(bwa_mem)
    print("Running: \n", bam_index)
    os.system(bam_index)
    os.remove(fastq_r1)
    os.remove(fastq_r2)


def qsub_sra2bam(f):
    fastq_name = os.path.basename(f).split(".sra")[0]
    fastq_r1 = os.path.join(args.output_dir, fastq_name + "_1.fastq.gz")
    fastq_r2 = os.path.join(args.output_dir, fastq_name + "_2.fastq.gz")
    ftmp = tempfile.NamedTemporaryFile()
    ftmp.write("#!/bin/bash\n")
    ftmp.write("#PBS -N " + fastq_name[-7:] + "-sra2bam\n")
    ftmp.write(
        "#PBS -o " + os.path.split(os.path.realpath(__file__))[0] + "/log\n")

    if args.node_name:
        ftmp.write(
            "#PBS -l nodes=1:" + args.node_name + ":ppn=1,walltime=100:00:00\n")
    else:
        ftmp.write(
            "#PBS -l nodes=1:ppn=1,walltime=100:00:00\n")

    ftmp.write("#PBS -j oe\ncd $PBS_O_WORKDIR\n")
    sra2fq = args.fastq_dump + " --split-3 --gzip " + \
        f + " -O " + args.output_dir
    bwa_mem = "bwa mem -t " + args.threads_number + \
        " -M -R \"@RG\\tID:" + fastq_name + "\\tLB:Hg19\\tPL:Illumina\\tPU:Barcode\\tSM:" + \
        fastq_name + "\\tCREATE_INDEX:True\" " + \
        args.ref_fa + " " + fastq_r1 + " " + fastq_r2 + " | samtools view -u -t " + \
        args.ref_fa + ".fai -S - -b | samtools sort - " + \
        os.path.join(args.output_dir, fastq_name + ".sorted")
    bam_index = "samtools index " + \
        os.path.join(args.output_dir, fastq_name + ".sorted.bam")
    ftmp.write(sra2fq + "\n")
    ftmp.write(bwa_mem + "\n")
    ftmp.write(bam_index)
    ftmp.seek(0)
    os.system("qsub " + ftmp.name)
    ftmp.close()


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


def main():
    global args
    args = get_args()
    pool = multiprocessing.Pool(processes=int(args.processes_number))
    makedir(os.path.split(os.path.realpath(__file__))[0] + "/log")
    if args.input_file:
        if args.qsub:
            qsub_sra2bam(args.input_file)
        else:
            sra2bam(args.input_file)
    else:
        if args.input_list:
            with open(args.input_list) as f:
                f_list = map(lambda x: x.strip(), f.readlines())
        elif args.input_re:
            f_list = glob.glob(args.input_re)

        if args.qsub:
            pool.map(qsub_sra2bam, f_list)
        else:
            pool.map(sra2bam, f_list)
        pool.close()
        pool.join()


if __name__ == '__main__':
    sys.exit(main())
