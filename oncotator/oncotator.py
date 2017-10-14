#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2017-7-8 14:09:48
# @Author  : Jintao Guo
# @Email   : guojt-4451@163.com

"""Annotation by Oncotator.

Example:
        $ oncotator.py -i input.vcf
                              --input-format VCF or MAFLITE
                              --qsub n1

# Input Format
# VCF: results of VarScan2
# MAFLITE: chr start end ref_allele alt_allele [others]
    # Single nucleotide variants
      chr4 150 150 A T
    # Insertions
      chr4 150 151 - T
    # Deletions
      chr4 150 150 A -
"""

import os
import sys
import tempfile
import linecache

version = "2.2"


def makedir(new_dir, exist_dir=None):
    """Make a directory.

    If it doesn't exist, handling concurrent race conditions.
    """
    if exist_dir:
        new_dir = os.path.join(exist_dir, new_dir)
    if os.path.exists(new_dir):
        print("The " + new_dir + " is already exist")
    else:
        print("Makedir " + new_dir)
        os.makedirs(new_dir)


def oncotator(input_file, input_format, output_format, dbdir, node=None):
    """Run oncotator in HPC. --infer-onps ===> DNP"""
    input_file = os.path.realpath(input_file)
    name, e = os.path.splitext(input_file)
    ftmp = tempfile.NamedTemporaryFile()
    ftmp.write(b"#!/bin/bash\n")
    if node:
        ftmp.write(b"#PBS -l nodes=1:" + node + "\n")
    ftmp.write(b"#PBS -j oe\n")
    export = b"source /etc/profile.d/set.sh\n"
    act = b"source activate oncotator\n"
    cmd = b"Oncotator -v -i " + input_format.encode('utf-8') + \
        b" --db-dir " + dbdir.encode('utf-8') + \
        b" -o " + output_format.encode('utf-8') + b" --skip-no-alt " + \
        input_file.encode('utf-8') + b" " + input_file.encode('utf-8') + \
        b".maf hg19\n"
    deact = b"source deactivate"
    ftmp.write(export)
    ftmp.write(act)
    ftmp.write(cmd)
    ftmp.write(deact)
    ftmp.seek(0)
    # print(ftmp.read())
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
                        version="%(prog)s " + version)
    parser.add_argument("-i",
                        metavar="file",
                        dest="file",
                        required=True,
                        help="input file")
    parser.add_argument("--input-format",
                        metavar="Format",
                        dest="input_format",
                        required=True,
                        choices=("VCF", "MAFLITE", "SEG_FILE"),
                        help="input file format [VCF or MAFLITE]")
    parser.add_argument("--output-format",
                        metavar="Format",
                        dest="output_format",
                        default="TCGAMAF",
                        choices=("TCGAMAF", "VCF", "SIMPLE_TSV",
                                 "TCGAVCF", "SIMPLE_BED", "GENE_LIST"),
                        help="output file format [TCGAMAF]")
    clustered_group = parser.add_argument_group("Clustered arguments")
    clustered_group.add_argument("--qsub",
                                 metavar="node",
                                 dest="node",
                                 type=str,
                                 help="name of nodes (e.g: n1,n2,...)")

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
    dbdir = "/share/data0/reference/Genome/oncotator_v1_ds_April052016"
    if(os.path.exists(linecache.getline(args.file, 1).strip())):
        with open(args.file) as f:
            file_list = map(lambda x: x.strip(), f.readlines())
        for line in file_list:
            if args.node:
                oncotator(line.strip(), args.input_format,
                          args.output_format, dbdir, args.node)
            else:
                oncotator(line.strip(), args.input_format,
                          args.output_format, dbdir)
    else:
        if args.node:
            oncotator(args.file, args.input_format,
                      args.output_format, dbdir, args.node)
        else:
            oncotator(args.file, args.input_format, args.output_format,
                      dbdir)


if __name__ == '__main__':
    sys.exit(main())
