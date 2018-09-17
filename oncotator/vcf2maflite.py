#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2016-10-12 20:55:26
# @Author  : Jintao Guo
# @Email   : guojt-4451@163.com
# @Version : 1.0

"""Convert vcf to input format.

# Single nucleotide variants
  chr4 150 150 A T
# Insertions
  chr4 150 151 - T
# Deletions
  chr4 150 150 A -
"""


import sys
import subprocess


version = "1.0"


def get_args3():
    """Get arguments from commond line."""
    try:
        import argparse
    except ImportError as imerr:
        print("\033[1;31m" + str(imerr) + " \033[0m")
        sys.exit()

    parser = argparse.ArgumentParser(
        prog="vcf2input",
        usage="%(prog)s",
        fromfile_prefix_chars='@',
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-v', '--version',
                        action='version',
                        version="%(prog)s " + version)

    parser.add_argument('-i', '--input',
                        metavar="FILE",
                        dest="input_file",
                        type=str,
                        help="input file for converting to input format")
    parser.add_argument('-o', '--output',
                        metavar="FILE",
                        dest="output_file",
                        type=str,
                        help="output file")

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


def snp(line):
    """Convert insert."""
    return("\t".join([":".join(line), line[0], line[1], line[1], line[2],
                      line[3]]))


def insert(line):
    """Convert insert."""
    chr = line[0]
    start = line[1]
    end = str(int(line[1]) + 1)
    ref = "-"
    alt = line[3][len(line[2]):]
    return("\t".join([":".join(line), chr, start, end, ref, alt]))


def deletion(line):
    """Convert deletion."""
    chr = line[0]
    start = line[1]
    ref = line[2][len(line[3]):]
    alt = "-"
    end = str(int(line[1]) + len(ref) - 1)
    return("\t".join([":".join(line), chr, start, end, ref, alt]))


def main():
    """Main."""
    args = get_args3()
    with open(args.input_file) as f:
        bed = map(lambda x: x.strip().split(), f.readlines())
    bed = list(filter(lambda x: x[3].find(",") < 0, bed))
    bed_snp = filter(lambda x: len(x[2]) == len(x[3]) == 1, bed)
    bed_snp = map(snp, bed_snp)
    bed_ins = filter(lambda x: len(x[2]) < len(x[3]), bed)
    bed_ins = map(insert, bed_ins)
    bed_del = filter(lambda x: len(x[2]) > len(x[3]), bed)
    bed_del = map(deletion, bed_del)
    list(bed_snp).extend(bed_ins)
    list(bed_snp).extend(bed_del)
    with open(args.output_file, "w") as f:
        f.write("ID\tchr\tstart\tend\tref_allele\talt_allele\n")
        f.write("\n".join(bed_snp))

if __name__ == '__main__':
    sys.exit(main())
