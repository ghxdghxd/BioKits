#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2016-12-14 19:40:00
# @Author  : Jintao Guo
# @Email   : guojt-4451@163.com
# @Version : 1.0

"""fastqc."""

import os
import sys
import tempfile

version = "1.0"


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
                        metavar="FILE",
                        dest="input",
                        required=True,
                        help="fastq file")
    parser.add_argument("-o",
                        metavar="DIR",
                        dest="output_dir",
                        default=os.environ['HOME'],
                        help="output dir or output file [" +
                        os.environ['HOME'] + "]")
    clustered_group = parser.add_argument_group("Clustered arguments")
    clustered_group.add_argument("--qsub",
                                 action="store_true",
                                 default=False,
                                 dest="qsub",
                                 help="run in cluster [False]")
    clustered_group.add_argument("--nodes",
                                 metavar="STR",
                                 dest="node",
                                 type=str,
                                 help="name of nodes (e.g: n1,n2,...)")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    else:
        args = parser.parse_args()
        return args


def fastqc(args):
    """Qsub."""
    name = os.path.basename(args.input)
    ftmp = tempfile.NamedTemporaryFile()
    ftmp.write(b"#!/bin/bash\n")
    ftmp.write(b"#PBS -N QC" + name.split('.')[0].encode() + b"\n")
    ftmp.write(
        b"#PBS -o " + args.output_dir.encode() + b"/log\n")
    ftmp.write(b"#PBS -j oe\ncd $PBS_O_WORKDIR\nsource /etc/profile.d/set.sh\n")
    ftmp.write(b"fastqc " + args.input.encode() + b" -o " + args.output_dir.encode())
    ftmp.seek(0)
    os.system("qsub " + ftmp.name)
    print(ftmp.read())
    ftmp.close()


def main():
    """Main."""
    args = get_args3()
    if(args.qsub):
        fastqc(args)

if __name__ == '__main__':
    sys.exit(main())
