#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2017-12-5 09:42:59
# @Author  : Jintao Guo
# @Email   : guojt-4451@163.com

import os
import sys
import tempfile
import multiprocessing

version = "2.0"


def main():
    global args
    args = get_args()
    pool = multiprocessing.Pool(processes=int(args.processes_number))

    if args.bam_list:
        with open(args.bam_list) as f:
            bam_list = map(lambda x: x.strip(), f.readlines())
            if args.qsub:
                pool.map(qsub_mark_duplicates, bam_list)
            else:
                pool.map(mark_duplicates, bam_list)
            pool.close()
            pool.join()
    elif args.bam_file:
        if args.qsub:
            qsub_mark_duplicates(args.bam_file)
        else:
            mark_duplicates(args.bam_file)


def get_args():
    """Get arguments from commond line"""
    try:
        import argparse
    except ImportError as imerr:
        print("\033[1;31m" + str(imerr) + " \033[0m")
        sys.exit()

    parser = argparse.ArgumentParser(usage="%(prog)s", fromfile_prefix_chars='@', description=__doc__)

    parser.add_argument("--path2picard",
                        metavar="DIR",
                        dest="picard",
                        default=".",
                        help="PATH to picard")

    group = parser.add_mutually_exclusive_group()
    group.add_argument("-i",
                       metavar="FILE",
                       dest="bam_file",
                       help="input bamfile")
    group.add_argument("-I",
                       metavar="FILE",
                       dest="bam_list",
                       help="list of input bamfile")

    parser.add_argument("-o",
                        metavar="DIR",
                        dest="output_dir",
                        default=os.environ['HOME'],
                        help="output dir or output file \
                        [" + os.environ['HOME'] + "]")
    parser.add_argument("-p",
                        metavar="INT",
                        default=1,
                        dest="processes_number",
                        type=int,
                        help="analyze multiple samples simultaneously [1]")
    parser.add_argument("-t",
                        metavar="INT",
                        default=1,
                        dest="threads_number",
                        type=int,
                        help="number of threads to allocate \
                        to each sample [1]")
    parser.add_argument("-m",
                        metavar="STR",
                        default="2g",
                        dest="Xmx",
                        help="memory [2g]")

    parser.add_argument("--qsub",
                        action="store_true",
                        default=False,
                        dest="qsub",
                        help="run crest in cluster")
    parser.add_argument("--remove_duplicates",
                        action="store_true",
                        default=False,
                        dest="rmdup",
                        help="remove PCR duplicates")
    parser.add_argument("--validation_stringency",
                        type=str,
                        default="STRICT",
                        metavar="{STRICT, LENIENT, SILENT}",
                        dest="validation_stringency",
                        help="Validation stringency for all SAM files \
                        read by this program [STRICT]")
    parser.add_argument("--node",
                        metavar="STR",
                        dest="node",
                        help="name of nodes")
    parser.add_argument("-n",
                        metavar="INT",
                        default=1,
                        dest="nodes_number",
                        type=int,
                        help="number of nodes [1]")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    else:
        args = parser.parse_args()
        args.picard = os.path.realpath(args.picard)
        args.output_dir = os.path.realpath(args.output_dir)
        args.threads_number = str(args.threads_number)
        args.processes_number = str(args.processes_number)
        return args


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


def mark_duplicates(bam):
    bam_name = os.path.basename(bam).split(".sorted.bam")[0]
    cmd = "java -XX:ParallelGCThreads=" + args.threads_number + \
        " -Xmx" + args.Xmx + \
        " -jar " + args.picard + "/MarkDuplicates.jar" +\
        " REMOVE_DUPLICATES=" + str(args.rmdup) +\
        " MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000" +\
        " VALIDATION_STRINGENCY=" + args.validation_stringency +\
        " INPUT=" + bam + \
        " OUTPUT=" + args.output_dir + "/" + bam_name + ".dedup.bam" +\
        " METRICS_FILE=" + args.output_dir + "/" + bam_name + ".metrics" +\
        " CREATE_INDEX=true ASSUME_SORTED=true"
    print(cmd)
    os.system(cmd)


def qsub_mark_duplicates(bam):
    bam_name = os.path.basename(bam).split(".sorted.bam")[0]
    ftmp = tempfile.NamedTemporaryFile()
    ftmp.write("#!/bin/bash\n")
    ftmp.write("#PBS -N " + bam_name[-10:] + "-PTMD\n")
    ftmp.write("#PBS -o " + os.path.split(os.path.realpath(__file__))[0] + "/log\n")
    if args.node:
        ftmp.write("#PBS -l nodes=1:" + args.node +
                   ":ppn=" + args.threads_number +
                   ",mem=" + args.Xmx + ",walltime=100:00:00\n")
    else:
        ftmp.write("#PBS -l nodes=1:ppn=" + args.threads_number +
                   ",mem=" + args.Xmx + ",walltime=100:00:00\n")
    ftmp.write("#PBS -j oe\ncd $PBS_O_WORKDIR\n")

    cmd = "java -XX:ParallelGCThreads=" + args.threads_number + \
        " -Xmx" + args.Xmx + \
        " -jar " + args.picard + "/MarkDuplicates.jar" +\
        " REMOVE_DUPLICATES=" + str(args.rmdup) +\
        " MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000" +\
        " VALIDATION_STRINGENCY=" + args.validation_stringency +\
        " INPUT=" + bam + \
        " OUTPUT=" + args.output_dir + "/" + bam_name + ".dedup.bam" +\
        " METRICS_FILE=" + args.output_dir + "/" + bam_name + ".metrics" +\
        " CREATE_INDEX=true ASSUME_SORTED=true"
    ftmp.write(cmd)
    ftmp.seek(0)
    # print ftmp.read()
    os.system("qsub " + ftmp.name)
    ftmp.close()


if __name__ == '__main__':
    sys.exit(main())
