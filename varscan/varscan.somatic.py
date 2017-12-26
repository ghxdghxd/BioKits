#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2017-12-7 17:41:23
# @Author  : Jintao Guo
# @Email   : guojt-4451@163.com

import os
import sys
import multiprocessing
import tempfile
import subprocess

usage = "%(prog)s [-h | --help] [-v | --version] \n\
                (PATH_to_varsacn) (ref_fa) (bam_paired) \n\
                [-l bed_file] [-o output_dir] [-p|--max-processors] [-c|--min-coverage] \n\
                [-m|--mpileup_options] [-m|--mpileup_options] [-m|--mpileup_options ...] \n\
                [--qsub] [--node-name]"

__version__ = " \n %(prog)s     V2.3"


def get_args3():
    """Get arguments from commond line."""
    try:
        import argparse as ap
    except ImportError as imerr:
        print("\n\033[1;31m" + str(imerr) + " \033[0m")
        sys.exit()

    parser = ap.ArgumentParser(usage=usage,
                               fromfile_prefix_chars='@',
                               description=__doc__,
                               formatter_class=ap.RawTextHelpFormatter)
    parser.add_argument('-v', '--version',
                        action='version',
                        version=__version__)
    parser.add_argument("varscan",
                        metavar="PATH_to_varsacn",
                        default=".",
                        help="PATH to varsacn.jar")
    parser.add_argument("ref_fa",
                        metavar="ref_fa",
                        help="faidx indexed reference sequence file")
    parser.add_argument("bam_paired",
                        metavar="bam_paired",
                        help="list of input bam_files (e.g: normal_bam tumor_bam)")
    parser.add_argument("-l",
                        metavar="bed_file",
                        dest="bed",
                        help="BED or position list file containing a list of regions or sites where pileup or "
                        "BCF should be generated [null]")
    parser.add_argument("-o",
                        metavar="output_dir",
                        dest="output_dir",
                        default=os.environ['HOME'],
                        help="output dir or output file [" +
                        os.environ['HOME'] + "]")
    parser.add_argument("-p", "--max-processors",
                        metavar="",
                        default=1,
                        dest="processes_number",
                        type=str,
                        help="analysis multiple samples simultaneously [1]")
    parser.add_argument("-c", "--min-coverage",
                        metavar="",
                        default=8,
                        dest="min_coverage",
                        type=str,
                        help="Minimum coverage in normal and tumor to call variant [8]")
    parser.add_argument("-m", "--mpileup_options",
                        metavar="",
                        dest="mpileup_option",
                        action="append",
                        help="samtools mpileup options (e.g: -m q30 -m C50)")

    clustered_group = parser.add_argument_group("Clustered arguments")
    clustered_group.add_argument("--qsub",
                                 action="store_true",
                                 default=False,
                                 dest="qsub",
                                 help="run crest in cluster [False]")
    clustered_group.add_argument("--node-name",
                                 metavar="",
                                 dest="node_name",
                                 type=str,
                                 help="name of nodes (e.g: compute-0-1,...)")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    else:
        args = parser.parse_args()
        check_dependencies(["samtools"])
        args.output_dir = os.path.realpath(args.output_dir)
        args.processes_number = str(args.processes_number)
        if args.mpileup_option:
            args.mpileup_option = "-" + " -".join(args.mpileup_option)
        return args


def check_dependencies(tools):
    """Ensure required tools are present."""
    print("Checking required dependencies......")
    if isinstance(tools, list):
        pass
    else:
        tools = [tools]
    for t in tools:
        subp = subprocess.Popen(["which", t], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        if subp.stderr.read():
            print("\033[1;31m" + "OSError: \033[1;33m" + __file__ + " requires " + t + "\033[0m")


def varscan(bams, run=False):
    """Running samtools mpileup + varscan2 mpileup2snp."""
    normal_bam = bams.split()[0]
    tumor_bam = bams.split()[1]
    output_name = os.path.basename(get_same_name(normal_bam, tumor_bam))
    cmd = ["samtools", "mpileup"]
    if args.mpileup_option:
        cmd.append(args.mpileup_option)
    if args.bed:
        cmd.append("-l " + args.bed)
    cmd.append("-f " + args.ref_fa)
    cmd.append(normal_bam)
    cmd.append(tumor_bam)
    cmd.append("| java -jar " + args.varscan)
    cmd.append("somatic - " + args.output_dir + "/" + output_name)
    cmd.append("--mpileup 1 --output-vcf 1 --min-coverage " + args.min_coverage)
    if(run):
        print(cmd)
        os.system(cmd)
    else:
        return(output_name, cmd)


def qsub_varscan(bams):
    """Qsub varscan2."""
    output_name, cmd = varscan(bams, run=False)
    ftmp = tempfile.NamedTemporaryFile()
    ftmp.write("#!/bin/bash\n")
    ftmp.write("#PBS -N " + output_name[-12:] + "-VS\n")
    ftmp.write("#PBS -o " + os.path.split(os.path.realpath(__file__))[0] + "/log\n")

    if args.node_name:
        ftmp.write("#PBS -l nodes=1:" + args.node_name + ":ppn=2,walltime=100:00:00\n")
    else:
        ftmp.write("#PBS -l nodes=1:ppn=2,walltime=100:00:00\n")

    ftmp.write("#PBS -j oe\ncd $PBS_O_WORKDIR\n")
    ftmp.write(cmd)
    ftmp.seek(0)
    # print cmd
    # print ftmp.read()
    os.system("qsub " + ftmp.name)
    ftmp.close()


def makedir(new_dir, exist_dir=None):
    """Make a directory.

    If it doesn't exist, handling concurrent race conditions.
    """
    if exist_dir:
        new_dir = os.path.join(exist_dir, new_dir)
    if os.path.exists(new_dir):
        print("The " + new_dir + " is already exist")
    else:
        print("Make " + new_dir)
        os.makedirs(new_dir)


def get_same_name(str1, str2):
    """Get same name."""
    for i in range(1, len(str1)):
        if str1[:-i] == str2[:-i] and str1[:-i][-1].isalnum():
            return os.path.basename(str1[:-i])
            break


def main():
    """Main."""
    global args
    args = get_args3()

    print(args.varscan)
    sys.exit()

    pool = multiprocessing.Pool(processes=int(args.processes_number))
    with open(args.bam_paired) as f:
        bam_list = map(lambda x: x.strip(), f.readlines())
        if args.qsub:
            makedir(os.path.split(os.path.realpath(__file__))[0] + "/log")
            pool.map(qsub_varscan, bam_list)
        else:
            pool.map(varscan, bam_list)
        pool.close()
        pool.join()


if __name__ == '__main__':
    sys.exit(main())
