#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2017-12-4 14:47:05
# @Author  : 郭金涛
# @Email   : guojt-4451@163.com
# @Version : 2.0

import glob
import os
import sys
import tempfile
import multiprocessing


def get_args():
    """get args"""
    try:
        import argparse
    except ImportError:
        print("\033[1;31m" + "ImportError: \033[0m" + __file__ + " requires `argparse`\n")
        sys.exit()

    parser = argparse.ArgumentParser(prog="sra2fq", usage="%(prog)s",
                                     fromfile_prefix_chars='@', description=__doc__)

    group = parser.add_argument_group("Dependencies")
    group.add_argument("--fastq_dump", metavar="DIR", dest="fastq_dump",
                       help="PATH to picard")

    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "-i", metavar="File", dest="input_file", help="input file")
    group.add_argument(
        "-I", metavar="Files", dest="input_list", help="list of input files")
    group.add_argument(
        "-r", metavar="Regular_expression", type=str, dest="input_re", help="input files")

    group = parser.add_argument_group("Output arguments")
    group.add_argument("-o", metavar="DIR", dest="output_dir",
                       default=os.environ['HOME'],
                       help="output dir or output file [" + os.environ['HOME'] + "]")
    parser.add_argument("-p", metavar="INT", default=1,
                        dest="processes_number", type=str, help="gzip multiple samples simultaneously  [1]")

    clustered_group = parser.add_argument_group("Clustered")
    clustered_group.add_argument("--qsub", action="store_true",
                                 default=False, dest="qsub", help="run crest in cluster [False]")
    clustered_group.add_argument("--nodes", metavar="STR", dest="node_name",
                                 type=str, help="name of nodes (e.g: n1,n2,...)")
    clustered_group.add_argument("-n", metavar="INT", default=1, dest="nodes_number",
                                 type=int, help="number of nodes [1]")

    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    else:
        return args


def sra2fq(f):
    cmd = args.fastq_dump + u" --split-3 --gzip " + \
        f + u" -O " + args.output_dir + u" --defline-qual '+' --defline-seq '@$ac-$si/$ri'"
    print(cmd)
    os.system(cmd)


def qsub_sra2fq(f):
    """Run sra2fq in nodes."""
    ftmp = tempfile.NamedTemporaryFile()
    ftmp.write(u"#!/bin/bash\n")
    # ftmp.write("#PBS -N " + jobname[-8:] + "-sra2fq\n")
    ftmp.write(
        u"#PBS -o " + os.path.split(os.path.realpath(__file__))[0] + u"/log\n")

    if args.node_name:
        ftmp.write(
            u"#PBS -l nodes=1:" + args.node_name + u":ppn=1,walltime=100:00:00\n")
    else:
        ftmp.write(
            u"#PBS -l nodes=1:ppn=1,walltime=100:00:00\n")

    ftmp.write(u"#PBS -j oe\ncd $PBS_O_WORKDIR\n")
    cmd = args.fastq_dump + u" --split-3 --gzip " + \
        f + u" -O " + args.output_dir
    ftmp.write(cmd)
    ftmp.seek(0)
    print(ftmp.read())
    # os.system(u"qsub " + ftmp.name)
    ftmp.close()


def makedir(new_dir, exist_dir=None):
    """Make a directory. If it doesn't exist, handling concurrent race conditions.
    """
    if exist_dir:
        new_dir = os.path.join(exist_dir, new_dir)
    if os.path.exists(new_dir):
        print(u"The " + new_dir + u" is already exist")
    else:
        print(u"Make " + new_dir)
        os.makedirs(new_dir)


def main():
    global args
    args = get_args()
    pool = multiprocessing.Pool(processes=int(args.processes_number))
    makedir(os.path.split(os.path.realpath(__file__))[0] + "/log")
    if args.input_file:
        if args.qsub:
            qsub_sra2fq(args.input_file)
        else:
            sra2fq(args.input_file)
    else:
        if args.input_list:
            with open(args.input_list) as f:
                f_list = map(lambda x: x.strip(), f.readlines())
        elif args.input_re:
            f_list = glob.glob(args.input_re)

        if args.qsub:
            pool.map(qsub_sra2fq, f_list)
        else:
            pool.map(sra2fq, f_list)
        pool.close()
        pool.join()


if __name__ == '__main__':
    sys.exit(main())
