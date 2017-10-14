#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2015-08-23 19:20:27
# @Author  : 郭金涛
# @Email   : guojt-4451@163.com

"""Run sample hotnet2."""

import os
import sys
import glob

version = "2.2"


def main():
    """Main."""
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
    """Get arguments from commond line."""
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
    parser.add_argument("--snv_file",
                        metavar="FILE",
                        dest="snv_file",
                        required=True,
                        help='Path to a tab-separated file containing SNVs'
                        ' where the first column of each line is a sample ID'
                        ' and subsequent columns contain the names of genes '
                        'with SNVs in that sample. Lines starting with "#" '
                        'will be ignored.')
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
                        type=int,
                        help="number of processes [1]")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    else:
        return parser.parse_args()


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


def hotnet2(args, influence_matrices):
    """Run hotnet2."""
    sample_cmd = "cut -f 1 " + args.snv_file + " |sort -u > " + \
        args.output_dir + "/sample.txt"
    print sample_cmd
    os.system(sample_cmd)
    gene_cmd = "cut -f 2 " + args.snv_file + " |sort -u > " + \
        args.output_dir + "/gene.txt"
    print gene_cmd
    os.system(gene_cmd)
    generateheat_cmd = "python " + args.hotnet2 + \
        "/generateHeat.py mutation --snv_file " + args.snv_file + \
        " --sample_file " + args.output_dir + "/sample.txt" + \
        " --gene_file " + args.output_dir + "/gene.txt" + \
        " --output_file " + args.output_dir + "/heat.json"
    print generateheat_cmd
    # os.system(generateheat_cmd)

    # results = []

    for net in influence_matrices:
        makedir(os.path.join(args.output_dir, net))
        simple_cmd = "python " + args.hotnet2 + "/runHotNet2.py --runname " + \
            net + "_SimpleRuns" + \
            " --network_name " + net + \
            " --infmat_file " + influence_matrices[net]["mat"] + \
            " --infmat_index_file " + influence_matrices[net]["index_genes"] +\
            " --heat_file " + args.output_dir + "/heat.json" +\
            " --permuted_networks_path " + \
            influence_matrices[net]["permuted"] +\
            " -o " + args.output_dir + "/" + net + \
            " --num_cores " + str(args.processes_number)
        print simple_cmd
        # os.system(simple_cmd)
        results_net = glob.glob(
            args.output_dir + "/" + net + "/*/results.json")
        makeresultswebsite_cmd = "python " + args.hotnet2 + \
            "/bin/makeResultsWebsite.py --results_files " + \
            " ".join(results_net) +\
            " --edge_file " + influence_matrices[net]["edge_list"] +\
            " --network_name " + net +\
            " --output_directory " + args.output_dir + "/" + net + "/viz"
        print makeresultswebsite_cmd
        # os.system(makeresultswebsite_cmd)
        # results.extend(results_net)

if __name__ == '__main__':
    sys.exit(main())
