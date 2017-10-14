#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2015-08-23 19:20:27
# @Author  : 郭金涛
# @Email   : guojt-4451@163.com

import os
import sys
import glob
import json
import numpy

version = "2.0"


def main():
    args = get_args()
    if args.hotnet2.endswith("/"):
        args.hotnet2 = args.hotnet2[:-1]
    influence_matrices_dir = glob.glob(args.hotnet2 + "/influence_matrices/*/")
    influence_matrices = {}
    for d in influence_matrices_dir:
        name = os.path.basename(os.path.dirname(d))
        influence_matrices[name] = {}
        influence_matrices[name]["edge_list"] = glob.glob(d + "/*edge_list")[0]
        influence_matrices[name]["index_genes"] = glob.glob(
            d + "/*index_genes")[0]
        influence_matrices[name]["mat"] = glob.glob(d + "/*mat")[0]
        influence_matrices[name]["permuted"] = glob.glob(
            d + "permuted/")[0] + "##NUM##/" + \
            os.path.basename(glob.glob(d + "/*mat")[0])
        influence_matrices[name]["num_permutations"] = str(len(
            glob.glob(d + "permuted/*/")))
    makedir(args.output_dir)
    hotnet2(args, influence_matrices)


def get_args():
    """Get arguments from commond line"""
    try:
        import argparse
    except ImportError as imerr:
        print "\033[1;31m" + str(imerr) + " \033[0m"
        sys.exit()

    parser = argparse.ArgumentParser(prog="hotnet2.py",
                                     usage="%(prog)s",
                                     version="%(prog)s " + version,
                                     fromfile_prefix_chars='@',
                                     description="description")

    parser.add_argument("--hotnet2", metavar="DIR", dest="hotnet2",
                        default=".", help="PATH to hotnet2")
    parser.add_argument("--mutsig_score_file",
                        metavar="FILE",
                        dest="mutsig_score_file",
                        required=True,
                        help='MutSig score file (gene to q-value).')

    parser.add_argument('--threshold',
                        metavar="FLOAT",
                        dest="threshold",
                        type=float,
                        default=1.0,
                        help='Maximum q-value threshold.')
    parser.add_argument("-o",
                        metavar="DIR",
                        dest="output_dir",
                        default=os.environ['HOME'],
                        help="output dir or output file \
                        [" + os.environ['HOME'] + "]")
    parser.add_argument("-p",
                        metavar="INT",
                        default=1,
                        dest="processes_number",
                        type=int,
                        help="number of processes [1]")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    else:
        return parser.parse_args()


def makedir(new_dir, exist_dir=None):
    """Make a directory. If it doesn't exist, handling concurrent race conditions.
    """
    if exist_dir:
        new_dir = os.path.join(exist_dir, new_dir)
    if os.path.exists(new_dir):
        print "The " + new_dir + " is already exist"
    else:
        print "Make " + new_dir
        os.makedirs(new_dir)


def hotnet2(args, influence_matrices):

    generateheat_cmd = "python " + args.hotnet2 + "/generateHeat.py mutsig --mutsig_score_file " + args.mutsig_score_file + \
        " --threshold " + str(args.threshold) + \
        " --output_file " + args.output_dir + \
        "/heat." + str(args.threshold) + ".json"

    print generateheat_cmd
    os.system(generateheat_cmd)

    for net in influence_matrices:
        makedir(os.path.join(args.output_dir, net))
        findthreshold_cmd = "python " + args.hotnet2 + "/bin/findThreshold.py network --runname " + net + "_DeltaSelection" +\
            " --infmat_index_file " + influence_matrices[net]["index_genes"] +\
            " --heat_file " + args.output_dir + "/heat." + str(args.threshold) + ".json" +\
            " --num_permutations " + influence_matrices[net]["num_permutations"] +\
            " --permuted_networks_path " + influence_matrices[net]["permuted"] +\
            " --output_file " + args.output_dir + "/" + net + "/delta." + net + "." + str(args.threshold) + ".json" + \
            " --num_cores " + str(args.processes_number)

        print findthreshold_cmd
        os.system(findthreshold_cmd)

        with open(args.output_dir + "/" + net + "/delta." + net + "." +
                  str(args.threshold) + ".json") as f:
            s = json.load(f)
            num = ["5", "10", "15", "20"]
            d = []
            for i in num:
                d.append(numpy.median(s["deltas"][i]))
            delta = str(numpy.min(d))
        findcomponents_cmd = "python " + args.hotnet2 + "/bin/findComponents.py --runname " + net + \
            " --infmat_file " + influence_matrices[net]["mat"] + \
            " --infmat_index_file " + influence_matrices[net]["index_genes"] +\
            " --infmat_name PPR --delta " + delta + \
            " --heat_file " + args.output_dir + "/heat." + str(args.threshold) + ".json" + \
            " --output_directory " + args.output_dir + "/" + net +\
            " heat --num_permutations " + \
            influence_matrices[net]["num_permutations"]

        print findcomponents_cmd
        os.system(findcomponents_cmd)

        makeresultswebsite_cmd = "python " + args.hotnet2 + "/bin/makeResultsWebsite.py --results_files " + args.output_dir + "/" + net + "/delta_" + delta + "/results.json" +\
            " --edge_file " + influence_matrices[net]["edge_list"] +\
            " --network_name " + net +\
            " --output_directory " + args.output_dir + "/" + net + \
            "/delta_" + delta + "/viz"

        print makeresultswebsite_cmd
        os.system(makeresultswebsite_cmd)

if __name__ == '__main__':
    sys.exit(main())
