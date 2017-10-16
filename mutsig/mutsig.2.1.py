#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2017-04-30 16:12:49
# @Author  : Jintao Guo
# @Email   : guojt-4451@163.com

"""Annotation by Oncotator.

# Input Format
# 1. MAF file or MAFLITE
Hugo_Symbol, Tumor_Sample_Barcode, Variant_Classification,
Chromosome, Start_position, End_position, Reference_Allele,
Tumor_Seq_Allele1, Tumor_Seq_Allele2
"""

import os
import sys
import tempfile

__version__ = "2.1"


def makedir(new_dir, exist_dir=None):
    """Make a directory.

    If it doesn't exist, handling concurrent race conditions.
    """
    if exist_dir:
        new_dir = os.path.join(exist_dir, new_dir)
    if os.path.exists(new_dir):
        pass
    else:
        os.makedirs(new_dir)


def qsub_mutsig(maf, node=None):
    """run mutsig."""
    ftmp = tempfile.NamedTemporaryFile()
    ftmp.write(b"#!/bin/bash\n")
    ftmp.write(b"#PBS -N Mutsig\n")
    if node:
        ftmp.write(b"#PBS -l nodes=1:" + node.encode('utf-8') + b"\n")
    else:
        ftmp.write(b"#PBS -l nodes=1\n")
    ftmp.write(b"#PBS -j oe\n")
    cmd = b"/share/apps/MutSigCV_1.4/run_MutSigCV.sh "
    cmd = cmd + b"/share/apps/MATLAB/MATLAB_Compiler_Runtime/v81 "
    cmd = cmd + maf.encode('utf-8')
    cmd = cmd + b" /share/apps/MutSigCV_1.4/exome_full192.coverage.txt "
    cmd = cmd + b"/share/apps/MutSigCV_1.4/gene.covariates.txt "
    cmd = cmd + maf.encode('utf-8') + b".mutsig "
    cmd = cmd + b"/share/apps/MutSigCV_1.4/mutation_type_dictionary_file.txt "
    cmd = cmd + b"/share/apps/MutSigCV_1.4/chr_files_hg19"
    ftmp.write(cmd)
    ftmp.seek(0)
    # print(ftmp.read())s
    os.system("qsub " + ftmp.name)
    ftmp.close()


def get_args3():
    """Get arguments from commond line."""
    try:
        import argparse as ap
    except ImportError as imerr:
        print("\033[1;31m" + str(imerr) + " \033[0m")
        sys.exit()

    parser = ap.ArgumentParser(usage="%(prog)s",
                               fromfile_prefix_chars='@',
                               description=__doc__,
                               formatter_class=ap.RawTextHelpFormatter)
    parser.add_argument('-v', '--version',
                        action='version',
                        version="%(prog)s\t" + __version__)

    parser.add_argument("-i",
                        metavar="FILE",
                        dest="file_maf",
                        required=True,
                        help="input maf file")

    clustered_group = parser.add_argument_group("Clustered arguments")
    clustered_group.add_argument("--qsub",
                                 metavar="STR",
                                 dest="node",
                                 type=str,
                                 help="name of nodes (e.g: n1,n2,...)")
    clustered_group.add_argument("-n",
                                 metavar="INT",
                                 dest="node_number",
                                 type=str,
                                 help="number of nodes [1]")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    else:
        args = parser.parse_args()
        return args


def main():
    """Main."""
    global args
    args = get_args3()
    args.file_maf = os.path.abspath(args.file_maf)
    if args.node:
        qsub_mutsig(args.file_maf, args.node)
    else:
        qsub_mutsig(args.file_maf)


if __name__ == '__main__':
    sys.exit(main())
