#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2016-06-05 17:41:23
# @Author  : Jintao Guo
# @Email   : guojt-4451@163.com

"""varscan2 copynumber.

variant detection in massively parallel sequencing data
"""

import os
import sys
import multiprocessing
import tempfile
import subprocess

version = "2.1"


def get_args():
    """Get arguments from commond line."""
    try:
        import argparse
    except ImportError as imerr:
        print "\033[1;31m" + str(imerr) + " \033[0m"
        sys.exit()

    parser = argparse.ArgumentParser(prog="varsacn.py",
                                     version="%(prog)s " + version,
                                     fromfile_prefix_chars='@',
                                     description=__doc__)
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
                                dest="bam_paired",
                                required=True,
                                help="list of input bam_files "
                                "(e.g: normal_bam tumor_bam)")
    parser.add_argument("-l",
                        metavar="STR",
                        dest="region",
                        help="regions or sites where pileup "
                        "should be generated [null]")
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
                        default=8,
                        dest="min_coverage",
                        type=str,
                        help="Minimum coverage in normal and tumor to call "
                        "variant [8]")
    parser.add_argument("--mpileup_option",
                        metavar="STR",
                        dest="mpileup_option",
                        help="samtools mpileup options (e.g:q30,C50)")

    clustered_group = parser.add_argument_group("Clustered arguments")
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

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    else:
        check_dependencies(["samtools"])
        args = parser.parse_args()
        args.output_dir = os.path.realpath(
            args.output_dir)
        args.processes_number = str(
            args.processes_number)
        args.mpileup_option = "-" + \
            " -".join(str(args.mpileup_option).split(","))
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


def varscan(bams):
    """Running varscan2."""
    normal_bam = bams.split()[0]
    tumor_bam = bams.split()[1]
    output_name = os.path.basename(get_same_name(normal_bam, tumor_bam))
    if args.mpileup_option:
        if args.region:
            cmd = "/home/jintao/samtools-0.1.18/bin/samtools mpileup " + \
                args.mpileup_option + \
                " -r " + args.region + \
                " -f " + args.ref_fa + " " + normal_bam + " " + tumor_bam + \
                " | java -jar " + args.varscan + \
                " copynumber - " + args.output_dir + "/" + output_name + \
                "." + str(args.region) + " --mpileup 1 --min-coverage " + \
                args.min_coverage + " --min-segment-size 100"
        else:
            cmd = "/home/jintao/samtools-0.1.18/bin/samtools mpileup " + \
                args.mpileup_option + \
                " -f " + args.ref_fa + " " + normal_bam + " " + tumor_bam + \
                " | java -jar " + \
                args.varscan + " copynumber - " + args.output_dir + "/" + \
                output_name + " --mpileup 1 --min-coverage " + \
                args.min_coverage + " --min-segment-size 100"
    else:
        if args.region:
            cmd = "samtools mpileup -r " + args.region + \
                " -f " + args.ref_fa + " " + normal_bam + " " + tumor_bam + \
                " | java -jar " + \
                args.varscan + " copynumber - " + args.output_dir + "/" + \
                output_name + "." + str(args.region) + \
                " --mpileup 1 --min-coverage " + \
                args.min_coverage + " --min-segment-size 100"
        else:
            cmd = "samtools mpileup -f " + args.ref_fa + \
                " " + normal_bam + " " + tumor_bam + " | java -jar " + \
                args.varscan + " copynumber - " + args.output_dir + "/" + \
                output_name + " --mpileup 1 --min-coverage " + \
                args.min_coverage + " --min-segment-size 100"
    print cmd
    # os.system(cmd)


def qsub_varscan(bams):
    """Qsub varscan2."""
    normal_bam = bams.split()[0]
    tumor_bam = bams.split()[1]
    output_name = os.path.basename(get_same_name(normal_bam, tumor_bam))

    if args.mpileup_option:
        if args.region:
            cmd = "/home/jintao/samtools-0.1.18/bin/samtools mpileup " + \
                args.mpileup_option + \
                " -r " + args.region + \
                " -f " + args.ref_fa + " " + normal_bam + " " + tumor_bam + \
                " | java -jar " + args.varscan + \
                " copynumber - " + args.output_dir + "/" + output_name + \
                "." + str(args.region) + " --mpileup 1 --min-coverage " + \
                args.min_coverage + " --min-segment-size 100"
            cmd1 = "java -jar " + args.varscan + " copyCaller " + \
                args.output_dir + "/" + output_name + "." + \
                str(args.region) + ".copynumber --output-file " + \
                args.output_dir + "/" + output_name + "." + \
                str(args.region) + ".copynumber.called"
        else:
            cmd = "/home/jintao/samtools-0.1.18/bin/samtools mpileup " + \
                args.mpileup_option + \
                " -f " + args.ref_fa + " " + normal_bam + " " + tumor_bam + \
                " | java -jar " + \
                args.varscan + " copynumber - " + args.output_dir + "/" + \
                output_name + " --mpileup 1 --min-coverage " + \
                args.min_coverage + " --min-segment-size 100"
            cmd1 = "java -jar " + args.varscan + " copyCaller " + \
                args.output_dir + "/" + output_name + \
                ".copynumber --output-file " + \
                args.output_dir + "/" + output_name + ".copynumber.called"
    else:
        if args.region:
            cmd = "samtools mpileup -r " + args.region + \
                " -f " + args.ref_fa + " " + normal_bam + " " + tumor_bam + \
                " | java -jar " + \
                args.varscan + " copynumber - " + args.output_dir + "/" + \
                output_name + "." + str(args.region) + \
                " --mpileup 1 --min-coverage " + \
                args.min_coverage + " --min-segment-size 100"
            cmd1 = "java -jar " + args.varscan + " copyCaller " + \
                args.output_dir + "/" + output_name + "." + \
                str(args.region) + ".copynumber --output-file " + \
                args.output_dir + "/" + output_name + "." + \
                str(args.region) + ".copynumber.called"
        else:
            cmd = "samtools mpileup -f " + args.ref_fa + \
                " " + normal_bam + " " + tumor_bam + " | java -jar " + \
                args.varscan + " copynumber - " + args.output_dir + "/" + \
                output_name + " --mpileup 1 --min-coverage " + \
                args.min_coverage + " --min-segment-size 100"
            cmd1 = "java -jar " + args.varscan + " copyCaller " + \
                args.output_dir + "/" + output_name + \
                ".copynumber --output-file " + \
                args.output_dir + "/" + output_name + ".copynumber.called"

    ftmp = tempfile.NamedTemporaryFile()
    ftmp.write("#!/bin/bash\n")
    ftmp.write("#PBS -N " + output_name[-12:] + "CNV\n")
    ftmp.write(
        "#PBS -o " + os.path.split(os.path.realpath(__file__))[0] + "/log\n")
    if args.node_name:
        ftmp.write("#PBS -l nodes=1:" + args.node_name +
                   ":ppn=1,walltime=100:00:00\n")
    else:
        ftmp.write("#PBS -l nodes=1:ppn=1,walltime=100:00:00\n")

    ftmp.write("#PBS -j oe\ncd $PBS_O_WORKDIR\n")
    ftmp.write("source /etc/profile.d/set.sh\n")
    ftmp.write(cmd + "\n")
    ftmp.write(cmd1)
    ftmp.seek(0)
    # print cmd
    print ftmp.read()
    # os.system("qsub " + ftmp.name)
    ftmp.close()


def makedir(new_dir, exist_dir=None):
    """Make a directory.

    If it doesn't exist, handling concurrent race conditions.
    """
    if exist_dir:
        new_dir = os.path.join(exist_dir, new_dir)
    if os.path.exists(new_dir):
        print "The " + new_dir + " is already exist"
    else:
        print "Make " + new_dir
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
    args = get_args()
    print args.varscan
    pool = multiprocessing.Pool(processes=int(args.processes_number))
    with open(args.bam_paired) as f:
        bam_list = map(lambda x: x.strip(), f.readlines())
        # print bam_list
        if args.qsub:
            makedir(os.path.split(os.path.realpath(__file__))[0] + "/log")
            pool.map(qsub_varscan, bam_list)
        else:
            pool.map(varscan, bam_list)
        pool.close()
        pool.join()

if __name__ == '__main__':
    sys.exit(main())
