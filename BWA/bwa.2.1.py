#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2016-09-17 21:47:05
# @Author  : Jintao Guo
# @Email   : guojt-4451@163.com

"""BWA (alignment via Burrows-Wheeler transformation)."""

import os
import sys
import subprocess
import multiprocessing
import tempfile

version = "2.1"

def main():
    """Main function."""
    global args
    args = get_args3()
    check_dependencies(["bwa", "samtools"])
    makedir(args.output_dir)
    print(args.output_dir)
    pool = multiprocessing.Pool(processes=int(args.processes_number))

    with open(args.fastqs) as fq:
        fastq_list = map(lambda x: x.strip(), fq.readlines())
        if args.qsub:
            pool.map(qsub_run_bwa, fastq_list)
        else:
            pool.map(run_bwa, fastq_list)
        pool.close()
        pool.join()

        # if len(tumor) % args.nodes == 0:
        #     tumor = zip(*([iter(tumor)]) * 20)
        #     normal = zip(*([iter(normal)]) * 20)
        # else:
        #     tmp = tumor[-(len(tumor) % 20):]
        #     tumor = zip(*([iter(tumor)]) * 20)
        #     tumor[-1] = list(tumor[-1]) + tmp
        #     print len(tumor[-1])
        #     tmp = tumor[-(len(tumor) % 20):]
        #     tumor = zip(*([iter(tumor)]) * 20)
        #     tumor.append(tuple(tmp))
        #     tmp = normal[-(len(normal) % 20):]
        #     normal = zip(*([iter(normal)]) * 20)
        #     normal.append(tuple(tmp))


def get_args3():
    """Get arguments from commond line."""
    try:
        import argparse
    except ImportError as imerr:
        print("\033[1;31m" + str(imerr) + " \033[0m")
        sys.exit()

    parser = argparse.ArgumentParser(prog="prog",
                                     usage="%(prog)s",
                                     fromfile_prefix_chars='@',
                                     description=__doc__)
    parser.add_argument('-v', '--version',
                        action='version',
                        version="%(prog)s " + version)

    required_group = parser.add_argument_group("Required arguments")
    required_group.add_argument("-r", "--ref_fa",
                                metavar="FILE",
                                dest="ref_fa",
                                required=True,
                                help="faidx indexed reference sequence file")
    required_group.add_argument("-i", "--input_list",
                                metavar="FILE",
                                dest="fastqs",
                                required=True,
                                help="list of input fastq_files " +
                                "(e.g: fastq_r1 fastq_r2)")

    parser.add_argument("-o", "--output_dir",
                        metavar="DIR",
                        dest="output_dir",
                        default=os.environ['HOME'],
                        help="output dir or output file " +
                        "[" + os.environ['HOME'] + "]")
    parser.add_argument("-p", "--processes_number",
                        metavar="INT",
                        default=1,
                        dest="processes_number",
                        type=int,
                        help="analyze multiple samples simultaneously [1]")
    parser.add_argument("-t", "--threads_number",
                        metavar="INT",
                        default=1,
                        dest="threads_number",
                        type=int,
                        help="number of threads to allocate to "
                        "each sample [1]")

    clustered_group = parser.add_argument_group("Clustered arguments")
    clustered_group.add_argument("-q", "--qsub",
                                 action="store_true",
                                 default=False,
                                 dest="qsub",
                                 help="run crest in cluster [False]")
    clustered_group.add_argument("--nodes_name",
                                 metavar="STR",
                                 dest="node",
                                 type=str,
                                 help="name of nodes (e.g: n1,n2,...)")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    else:
        args = parser.parse_args()
        args.output_dir = os.path.realpath(args.output_dir)
        args.threads_number = str(args.threads_number)
        args.processes_number = str(args.processes_number)
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


def run_bwa(fastqs, run = True):
    """Run bwa."""
    fastq_r1 = fastqs.split()[0]
    fastq_r2 = fastqs.split()[1]
    fastq_name = get_same_name(fastq_r1, fastq_r2)
    bwa_mem = b"bwa mem -t " + args.threads_number.encode() + \
        b" -M -R \"@RG\\tID:" + fastq_name.encode() + \
        b"\\tLB:HG\\tPL:Illumina\\tPU:Barcode\\tSM:" + fastq_name.encode() + \
        b"\\tCREATE_INDEX:True\" " + \
        args.ref_fa.encode() + b" " + fastq_r1.encode() + b" " + fastq_r2.encode() + \
        b" | samtools view -u -t " + args.ref_fa.encode() + \
        b".fai -S - -b | samtools sort -@" + args.threads_number.encode() + b" - -o " + \
        args.output_dir.encode() + b"/" + fastq_name.encode() + b".sorted.bam\n"
    bam_index = b"samtools index " + args.output_dir.encode() + \
        b"/" + fastq_name.encode() + b".sorted.bam"
    if run:
        print("Running: ", bwa_mem)
        os.system(bwa_mem)
        print("Running: ", bam_index)
        os.system(bam_index)
    else:
        return fastq_name, bwa_mem, bam_index


def qsub_run_bwa(fastqs):
    """Run bwa in clusters."""
    fastq_name, bwa_mem, bam_index = run_bwa(fastqs, run=False)
    # fastq_r1 = fastqs.split()[0]
    # fastq_r2 = fastqs.split()[1]
    # fastq_name = get_same_name(fastq_r1, fastq_r2)
    ftmp = tempfile.NamedTemporaryFile()
    ftmp.write(b"#!/bin/bash\n")
    ftmp.write(b"#PBS -N bwa-" + fastq_name.encode() + b"\n")
    ftmp.write(b"#PBS -o " + args.output_dir.encode() + b"/" + fastq_name.encode() + b"\n")
    if args.node:
        ftmp.write(b"#PBS -l nodes=1:" + args.node.encode() +
                   b":ppn=" + args.threads_number.encode() + b",walltime=100:00:00\n")
    else:
        ftmp.write(b"#PBS -l nodes=1:ppn=" + args.threads_number.encode() + b",walltime=100:00:00\n")
    ftmp.write(b"#PBS -j oe\ncd $PBS_O_WORKDIR\n")

    # bwa_mem = b"bwa mem -t " + args.threads_number.encode() + \
    #     b" -M -R \"@RG\\tID:" + fastq_name.encode() + \
    #     b"\\tLB:HG\\tPL:Illumina\\tPU:Barcode\\tSM:" + fastq_name.encode() + \
    #     b"\\tCREATE_INDEX:True\" " + \
    #     args.ref_fa.encode() + b" " + fastq_r1.encode() + b" " + fastq_r2.encode() + \
    #     b" | samtools view -u -t " + args.ref_fa.encode() + \
    #     b".fai -S - -b | samtools sort - " + \
    #     args.output_dir.encode() + b"/" + fastq_name.encode() + b".sorted\n"
    # bam_index = b"samtools index " + args.output_dir.encode() + \
    #     b"/" + fastq_name.encode() + b".sorted.bam"
    ftmp.write(b"source /etc/profile.d/set.sh\n")
    ftmp.write(bwa_mem)
    ftmp.write(bam_index)
    ftmp.seek(0)
    # print(ftmp.read())
    os.system("qsub " + ftmp.name)
    ftmp.close()


def get_same_name(str1, str2):
    """Get the same str between two names."""
    for i in range(1, len(str1)):
        if str1[:-i] == str2[:-i] and str1[:-i][-1].isalnum():
            return os.path.basename(str1[:-i])
            # break


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

if __name__ == '__main__':
    sys.exit(main())
