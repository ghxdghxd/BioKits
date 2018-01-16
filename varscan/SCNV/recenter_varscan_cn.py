#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2016-09-26 20:13:43
# @Author  : Jintao Guo
# @Email   : guojt-4451@163.com
# @Version : 1.0

"""Recenter copynumber."""

import os
import sys
import pandas as pd
import tempfile
import subprocess

version = "1.0"


def get_args3():
    """Get arguments from commond line."""
    try:
        import argparse
    except ImportError as imerr:
        print("\033[1;31m" + str(imerr) + " \033[0m")
        sys.exit()

    parser = argparse.ArgumentParser(usage="%(prog)s",
                                     fromfile_prefix_chars='@',
                                     description=__doc__)
    parser.add_argument('-v', '--version',
                        action='version',
                        version="%(prog)s " + version)

    parser.add_argument("-i",
                        metavar="File",
                        dest="input_called",
                        help="copynumber called file produced by "
                        "varscan 'copyCaller'")
    parser.add_argument("-o",
                        metavar="Dir",
                        dest="output_dir",
                        help="output dir")
    parser.add_argument("--min-coverage",
                        metavar="INT",
                        default=8,
                        dest="min_coverage",
                        type=str,
                        help="Minimum coverage in normal and tumor to call "
                        "variant [8]")
    parser.add_argument("--min-segment-size",
                        metavar="INT",
                        default=20,
                        dest="min_segment_size",
                        type=str,
                        help="Minimum number of consecutive bases to report a "
                        "segment [20]")
    parser.add_argument("--max-segment-size",
                        metavar="INT",
                        default=100,
                        dest="max_segment_size",
                        type=str,
                        help="Max size before a new segment is made [100]")

    clustered_group = parser.add_argument_group("Clustered arguments")
    clustered_group.add_argument("--qsub",
                                 action="store_true",
                                 default=False,
                                 dest="qsub",
                                 help="run in cluster [False]")
    clustered_group.add_argument("--node",
                                 metavar="STR",
                                 dest="node",
                                 type=str,
                                 help="name of nodes (e.g: n1,n2,...)")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    else:
        # check_dependencies(["group_name"])
        args = parser.parse_args()
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


def cmd_copycaller(input_copynumber, output_called, recenter=None):
    """Get copyCaller command."""
    cmd = "java -jar /share/apps/VarScan/VarScan.v2.3.9.jar copyCaller " + \
        input_copynumber
    if recenter > 0:
        cmd = cmd + " --recenter-up " + recenter + " --output-file " + \
            output_called + ".copynumber.recentered.called"
    elif recenter < 0:
        recenter = 0 - recenter
        cmd = cmd + " --recenter-down " + recenter + " --output-file " + \
            output_called + ".copynumber.recentered.called"
    else:
        cmd = cmd + " --output-file " + output_called + \
            ".copynumber.recentered.called"
    cmd = cmd + " --min-coverage 20 --amp-threshold 0.9" + \
        " --del-threshold 0.9 --min-region-size 20"
    return(cmd)


# def cmd_recenter(input_called):
#     """"Get recenter command."""
#     cmd = "awk 'BEGIN{num=sum=0}{if(NR>1){num+=$4;sum+=$4*$7}}" \
#         "END{printf(\"%.16f\",sum/num)}' " + input_called
#     return(cmd)


def cmd_dnacopy(input_recentered_called):
    """DNAcopy."""
    cmd = "Rscript /home/jintao/varsan_CNV/varscan2_copynumber.R " + \
        input_recentered_called
    return(cmd)


# def qsub(cmds):
#     """Qsub."""
#     ftmp = tempfile.NamedTemporaryFile()
#     ftmp.write("#!/bin/bash\n")
#     ftmp.write("#PBS -N jobname CNV\n")
#     ftmp.write(
#         "#PBS -o " + os.path.split(os.path.realpath(__file__))[0] + "/log\n")
#     ftmp.write("#PBS -l nodes=1:ppn=1,walltime=10:00:00\n")
#     if jobid:
#         ftmp.write("#PBS -W depend=afterok:" + jobid.readlines()[0])
#     ftmp.write("#PBS -j oe\ncd $PBS_O_WORKDIR\n")
#     ftmp.write("python ")
#     ftmp.seek(0)
#     printftmp.read()
#     # os.system("qsub " + ftmp.name)
#     ftmp.close()


def run():
    """Run."""
    cmd_copycaller = cmd_copycaller()
    cmd_recenter = cmd_recenter()
    cmd_dnacopy = cmd_dnacopy()


def main():
    """Main."""
    args = get_args3()
    centerinfo = pd.DataFrame.from_csv(
        args.input_called, sep="\t", index_col=None)
    recenter = sum(centerinfo['num_positions'] *
                   centerinfo['adjusted_log_ratio']) / \
        sum(centerinfo['num_positions'])

    cmd = cmd_copycaller(args.input_called, args.output_dir, recenter)

    print(cmd)
    # if args.qsub:
    #     script = "python3 " + os.path.join(os.getcwd(),
    #                                        os.path.basename(__file__))
    #     print(script)
    #     print(args)
    #     # subprocess.call(script, shell=True)
    # else:
    #     run()

if __name__ == '__main__':
    sys.exit(main())
