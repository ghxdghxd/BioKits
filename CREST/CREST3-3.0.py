#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2015-08-21 22:05:39
# @Author  : 郭金涛
# @Email   : guojt-4451@163.com


"""
CREST (Clipping Reveals Structure) is an algorithm for detecting genomic
structural variations at base-pair resolution using next-generation sequencing
data.
"""

import os
import sys
import subprocess
import multiprocessing
import time

version = "3.0"


def main():
    args = get_args3()
    if args.qsub:
        makedir(os.path.split(os.path.realpath(__file__))[0] + "/log")
    else:
        makedir(args.output_dir)
        run_processes_crest(args)


def get_args3():
    """get args
    """
    try:
        import argparse
    except ImportError:
        raise ImportError(__file__ + " requires `argparse`")

    parser = argparse.ArgumentParser(prog="crest.py",
                                     usage="%(prog)s",
                                     description=__doc__,
                                     fromfile_prefix_chars='@')
    parser.add_argument('-v', '--version',
                        action='version',
                        version="%(prog)s " + version)
    group = parser.add_argument_group("Dependencies")
    group.add_argument("--path2crest",
                       metavar="DIR",
                       dest="crest",
                       default=".",
                       help="PATH to CREST")
    group.add_argument("--path2cap3",
                       metavar="FILE",
                       dest="cap3",
                       default=".",
                       help="PATH to CAP3")

    group = parser.add_argument_group("Options")
    group.add_argument("-r",
                       metavar="FILE",
                       dest="ref_fa",
                       help="faidx indexed reference sequence file")
    group.add_argument("-t",
                       metavar="FILE",
                       dest="ref_2bit",
                       help="The reference genome in 2bit format")
    group.add_argument("-i",
                       metavar="FILE",
                       dest="bams_file",
                       help="list of input bamfiles [tumorBam \\t normalBam]")
    group.add_argument("-o",
                       metavar="DIR",
                       dest="output_dir",
                       default=os.environ['HOME'],
                       help="output dir or output file \
                       [" + os.environ['HOME'] + "]")
    group.add_argument("-p",
                       metavar="INT",
                       default=1,
                       dest="processes_number",
                       type=int,
                       help="number of processes [1]")

    group = parser.add_argument_group("Cluster")
    group.add_argument("--qsub",
                       action="store_true",
                       default=False,
                       dest="qsub",
                       help="run crest in cluster")
    group.add_argument("-n",
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
        check_dependencies("blat2")
        return args


def check_dependencies(tools):
    """Ensure required tools are present.
    """
    print("Checking required dependencies......")
    if not isinstance(tools, list):
        tools = [tools]
    for t in tools:
        try:
            print(subprocess.check_output(["which", t]).strip())
        except OSError:
            raise OSError("\033[1;31m" + __file__ + " requires " +
                          t + "\033[0m")
            sys.exit()


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


def run_crest(cmd_list):
    for cmd in cmd_list:
        print(cmd)
        subprocess.call(cmd, shell=True)


def run_processes_crest(args):
    """run crest in a computer
    """
    gfserver_start = "gfServer start localhost 50000 -stepSize=5 -canStop -log=blatServer.log " + \
        args.ref_2bit + " &"
    print(gfserver_start)
    os.system(gfserver_start)

    gfserver_check = "gfServer files localhost 50000"
    while 1:
        subp = subprocess.Popen(gfserver_check,
                                stderr=subprocess.PIPE,
                                stdout=subprocess.PIPE,
                                shell=True)
        if subp.stderr.read():
            print("\033[1;31m Waiting for BLAT server to start...\033[0m")
            time.sleep(10)
        else:
            print("\033[1;31m BLAT server is running!\033[0m")
            break

    extractsclip_tumor_cmd = ["perl",
                              os.path.join(os.path.abspath(args.crest),
                                           "extractSClip.pl"),
                              "-i", "tumor_bam",
                              "--ref_genome", args.ref_fa,
                              "-o", os.path.abspath(args.output_dir)]

    extractsclip_normal_cmd = ["perl",
                               os.path.join(os.path.abspath(args.crest),
                                            "extractSClip.pl"),
                               "-i", "normal_bam",
                               "--ref_genome", args.ref_fa,
                               "-o", os.path.abspath(args.output_dir)]

    countdiff_cmd = ["perl",
                     os.path.join(os.path.abspath(args.crest),
                                  "countDiff.pl"),
                     "-d", "tumor_bam_cover",
                     "-g", "normal_bam_cover"]

    crest_cmd = ["perl",
                 os.path.join(os.path.abspath(args.crest), "CREST.pl"),
                 "-f", "somatic_cover",
                 "-d", "tumor_bam",
                 "-g", "normal_bam",
                 "--ref_genome", args.ref_fa,
                 "-t", args.ref_2bit,
                 "-o", os.path.abspath(args.output_dir),
                 "--blatserver localhost --blatport 50000"]

    bam2html_cmd = ["perl",
                    os.path.join(os.path.abspath(args.crest), "bam2html.pl"),
                    "-d", "tumor_bam",
                    "-g", "normal_bam",
                    "-i", "somatic_predSV",
                    "--ref_genome", "ref_fa",
                    "-o", "tumor_predSV.html"]

    pool = multiprocessing.Pool(processes=args.processes_number)

    with open(args.bams_file) as bams_file:
        for line in bams_file.readlines():
            bams = line.strip().split()
            extractsclip_tumor_cmd[3] = crest_cmd[5] = \
                bam2html_cmd[3] = bams[0]

            bam2html_cmd[7] = os.path.join(
                os.path.abspath(args.output_dir),
                os.path.basename(bams[0]) + ".predSV.txt")
            bam2html_cmd[11] = os.path.join(
                os.path.abspath(args.output_dir),
                os.path.basename(bams[0]) + ".predSV.html")
            if len(bams) == 1:
                crest_cmd[3] = os.path.join(os.path.abspath(args.output_dir),
                                            os.path.basename(bams[0]) +
                                            ".cover")
                pool.apply_async(run_crest,
                                 [" ".join(extractsclip_tumor_cmd),
                                  " ".join(crest_cmd[:6] + crest_cmd[8:]),
                                  " ".join(bam2html_cmd[:4] +
                                           bam2html_cmd[6:])])
            elif len(bams) == 2:
                extractsclip_normal_cmd[3] = crest_cmd[7] = \
                    bam2html_cmd[5] = bams[1]
                countdiff_cmd[3] = os.path.join(os.path.abspath(
                    args.output_dir), os.path.basename(bams[0]) + ".cover")
                countdiff_cmd[5] = os.path.join(os.path.abspath(
                    args.output_dir), os.path.basename(bams[1]) + ".cover")
                crest_cmd[3] = os.path.join(os.path.abspath(args.output_dir),
                                            os.path.basename(bams[0]) +
                                            ".cover.somatic.cover")
                pool.apply_async(run_crest,
                                 ([" ".join(extractsclip_tumor_cmd),
                                   " ".join(extractsclip_normal_cmd),
                                   " ".join(countdiff_cmd),
                                   " ".join(crest_cmd),
                                   " ".join(bam2html_cmd)], ))
        pool.close()
        pool.join()
    gfserver_stop = "gfServer stop localhost 50000"
    print(gfserver_stop)
    os.system(gfserver_stop)

if __name__ == "__main__":
    sys.exit(main())
