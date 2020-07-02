#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2016-10-14 17:41:23
# @Author  : Jintao Guo
# @Email   : guojt-4451@163.com

"""Extract genotype from bam using samtools(mpileup) + varscan(mpileup2cns).
"""

import os
import sys
import multiprocessing
import tempfile
import subprocess

version = "1.0"
__version__ = "%(prog)s "

def get_args3():
    """Get arguments from commond line."""
    try:
        import argparse as ap
    except ImportError as imerr:
        print("\033[1;31m" + str(imerr) + " \033[0m")
        sys.exit()

    parser = ap.ArgumentParser(prog="bam2GT",
                               fromfile_prefix_chars='@',
                               description=__doc__,
                               formatter_class=ap.RawTextHelpFormatter)
    parser.add_argument('-v', '--version',
                        action='version',
                        version=__version__)

    required_group = parser.add_argument_group("Required arguments")
    required_group.add_argument("varscan",
                                metavar="varsacn",
                                default=".",
                                help="PATH to varsacn.jar")
    required_group.add_argument("-r",
                                metavar="FILE",
                                dest="ref_fa",
                                required=True,
                                help="faidx indexed reference sequence file")
    required_group.add_argument("-i",
                                metavar="FILE",
                                dest="bam_list",
                                required=True,
                                help="list of input bam_files")

    parser.add_argument("-l",
                        metavar="FILE",
                        dest="bed",
                        help="BED or position list file containing "
                        "a list of regions or sites where pileup or "
                        "BCF should be generated [null]")
    parser.add_argument("-o",
                        metavar="DIR",
                        dest="output_dir",
                        default=os.environ['HOME'],
                        help="output dir or output file [" +
                        os.environ['HOME'] + "]")
    parser.add_argument("-p",
                        metavar="INT",
                        default=1,
                        dest="processes_number",
                        type=str,
                        help="analysis multiple samples simultaneously [1]")
    parser.add_argument("--min-coverage",
                        metavar="INT",
                        default=10,
                        dest="min_coverage",
                        type=str,
                        help="Minimum coverage in normal and tumor to call "
                        "variant [10]")
    parser.add_argument("--min-reads2",
                        metavar="INT",
                        default=10,
                        dest="min_reads2",
                        type=str,
                        help="Minimum coverage in normal and tumor to call "
                        "variant [2]")
    parser.add_argument("--mpileup_option",
                        metavar="STR",
                        dest="mpileup_option",
                        help="samtools mpileup options (e.g:q30,C50)")

    clustered_group = parser.add_argument_group("Clustered arguments")
    clustered_group.add_argument("--qsub",
                                 action="store_true",
                                 default=False,
                                 dest="qsub",
                                 help="run in cluster [False]")
    clustered_group.add_argument("-n",
                                 metavar="INT",
                                 dest="node_number",
                                 type=str,
                                 help="number of nodes [1]")
    clustered_group.add_argument("--nodes",
                                 metavar="STR",
                                 dest="node",
                                 type=str,
                                 help="name of nodes (e.g: n1,n2,...)")

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()
    else:
        check_dependencies(["samtools"])
        args = parser.parse_args()
        args.output_dir = os.path.realpath(args.output_dir)
        args.processes_number = str(args.processes_number)
        args.mpileup_option = "-" + " -".join(str(args.mpileup_option).split(","))
        return args



def check_dependencies(tools):
    """Ensure required tools are present."""
    if not isinstance(tools, list):
        tools = [tools]
    try:
        for tool in tools:
            try:
                subprocess.check_output(["which", tool]).strip()
            except subprocess.CalledProcessError:
                print("\033[1;31m" + __file__ + " requires " +
                      tool + "\033[0m")
    except:
        sys.exit()

def bam2gt(bam):
    """Running samtools mpileup + varscan2 mpileup2snp."""
    #output_name = os.path.basename(bam).split(".")[0][:-1]
    output_name = os.path.basename(bam).split(".")[0]
    dir_name = os.path.dirname(bam)
    cmd = "samtools mpileup "
    if args.mpileup_option:
        cmd = cmd + args.mpileup_option
    if args.bed:
        cmd = cmd + " -l " + args.bed
    cmd = cmd + " -f " + args.ref_fa + " " + bam + \
        " | java -jar " + args.varscan + " mpileup2cns" + \
        " - --output-vcf 1 --min-coverage " + args.min_coverage + \
        " --min-reads2 " + args.min_reads2 + \
        " | bgzip > " + dir_name + "/" + output_name + ".cns.vcf.bz"
    tabix_cmd = "tabix" + dir_name + "/" + output_name + ".cns.vcf.bz"
    #print([output_name, cmd])
    # return([output_name, cmd])
    os.system(cmd)
    os.system(tabix_cmd)


def qsub_bam2gt(bam):
    """Qsub varscan2."""
    output_name = os.path.basename(bam).split(".")[0]
    dir_name = os.path.dirname(bam)
    cmd = "samtools mpileup "
    if args.mpileup_option:
        cmd = cmd + args.mpileup_option
    if args.bed:
        cmd = cmd + " -l " + args.bed
    cmd = cmd + " -f " + args.ref_fa + " " + bam + \
        " | java -jar " + args.varscan + " mpileup2cns" + \
        " - --output-vcf 1 --min-coverage " + args.min_coverage + \
        " --min-reads2 " + args.min_reads2 + \
        " | bgzip > " + dir_name + "/" + output_name + ".cns.vcf.bz"
    tabix_cmd = "tabix" + dir_name + "/" + output_name + ".cns.vcf.bz"

    ftmp = tempfile.NamedTemporaryFile()
    ftmp.write("#!/bin/bash\n")
    ftmp.write("#PBS -N " + cmd[0] + "-cns-GT\n")
    ftmp.write(
        "#PBS -o " + os.path.split(os.path.realpath(__file__))[0] + "/log\n")

    if args.node_name:
        ftmp.write("#PBS -l nodes=1:" + args.node_name +
                   ":ppn=2,walltime=100:00:00\n")
    else:
        ftmp.write("#PBS -l nodes=1:ppn=2,walltime=100:00:00\n")

    ftmp.write("#PBS -j oe\ncd $PBS_O_WORKDIR\n")
    ftmp.write("source /etc/profile.d/set.sh\n")
    ftmp.write(cmd[1])
    ftmp.write(tabix_cmd)
    ftmp.seek(0)
    # print(cmd)
    # print(ftmp.read())
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


# def get_same_name(str1, str2):
#     """Get same name."""
#     for i in range(1, len(str1)):
#         if str1[:-i] == str2[:-i] and str1[:-i][-1].isalnum():
#             return os.path.basename(str1[:-i])


def main():
    """Main."""
    global args
    args = get_args3()
    # print(args)
    # sys.exit()
    pool = multiprocessing.Pool(processes=int(args.processes_number))
    with open(args.bam_list) as f:
        bam_list = map(lambda x: x.strip(), f.readlines())
        if args.qsub:
            makedir(os.path.split(os.path.realpath(__file__))[0] + "/log")
            pool.map(qsub_bam2gt, bam_list)
        else:
            # print(bam2gt(list(bam_list)[0]))
            pool.map(bam2gt, bam_list)
        pool.close()
        pool.join()

if __name__ == '__main__':
    sys.exit(main())
