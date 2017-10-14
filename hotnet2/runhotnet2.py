#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date     : 2017-03-18 19:20:27
# @Author   : Jintao Guo
# @Email    : guojt-4451@163.com

"""Run hotnet2.

Return the sig genes and networks
"""

import os
import sys
import glob
try:
    import argparse
    import json
    import pandas as pd
except ImportError as imerr:
    print("\033[1;31m" + str(imerr) + " \033[0m")
    sys.exit()

__version__ = "3.1"


def main():
    """Main."""
    args = get_args3()
    args.hotnet2 = args.hotnet2.rstrip("/")
    args.output_dir = args.output_dir.rstrip("/")
    print(args)
    if args.maf_file is not None and args.mutsig_file is not None:
        snv_genes = filtergene(args.mutsig_file, args.threshold)
        samples, snv_mat = maf2snv(args.maf_file, snv_genes)
    elif args.maf_file is not None and args.mutsig_file is None:
        samples, snv_mat = maf2snv(args.maf_file)
    elif args.mutsig_file is not None and args.maf_file is None:
        pass
    makedir(args.output_dir)
    with open(args.output_dir + "/snv." + args.threshold, "w") as f:
        f.write("\n".join(["\t".join(x) for x in snv_mat]))

    snv_file = args.output_dir + "/snv." + args.threshold

    influence_matrices_dir = glob.glob(args.hotnet2 + "/influence_matrices/*/")
    influence_matrices = {}
    for directory in influence_matrices_dir:
        name = os.path.basename(os.path.dirname(directory))
        influence_matrices[name] = {}
        influence_matrices[name]["edge_list"] = glob.glob(
            directory + "/*edge_list")[0]
        influence_matrices[name]["index_genes"] = glob.glob(
            directory + "/*index_genes")[0]
        influence_matrices[name]["mat"] = glob.glob(directory + "/*mat")[0]
        influence_matrices[name]["permuted"] = glob.glob(
            directory + "permuted/")[0] + "##NUM##/" + \
            os.path.basename(glob.glob(directory + "/*mat")[0])
        influence_matrices[name]["num_permutations"] = str(len(
            glob.glob(directory + "permuted/*/")))
    hotnet2(args, influence_matrices, snv_file)
    visnetwork(args, snv_file)


def get_args3():
    """Get arguments from commond line."""
    parser = argparse.ArgumentParser(prog="hotnet2.py",
                                     usage="%(prog)s",
                                     fromfile_prefix_chars='@',
                                     description=__doc__)
    parser.add_argument('-v', '--version',
                        action='version',
                        version="%(prog)s " + __version__)
    parser.add_argument("--hotnet2", metavar="DIR", dest="hotnet2",
                        default=".", help="PATH to hotnet2")
    parser.add_argument("--maf_file",
                        metavar="FILE",
                        type=str,
                        dest="maf_file",
                        default=None,
                        help='Path to a maf file')
    parser.add_argument("--mutsig_file",
                        metavar="FILE",
                        dest="mutsig_file",
                        default=None,
                        help='MutSig score file (gene to q-value).')
    parser.add_argument('--threshold',
                        metavar="FLOAT",
                        dest="threshold",
                        type=str,
                        default="p0.05",
                        help='Maximum p-value or q-value threshold.'
                        '(e.g. p0.05 or q0.05)')
    parser.add_argument("--cna_file",
                        metavar="FILE",
                        dest="cna_file",
                        type=str,
                        default=None,
                        help='Path to a tab-separated file containing CNAs '
                        'where the first column of each line is a sample ID '
                        'and subsequent columns contain gene names followed '
                        'by "(A)" or "(D)" indicating an amplification or '
                        'deletion in that gene for the sample. Lines starting '
                        'with "#" will be ignored.')
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
        print("The " + new_dir + " is already exist")
    else:
        print("Make " + new_dir)
        os.makedirs(new_dir)


def filtergene(mutsig_file, threshold):
    """Filter genes from mutsig with threshold."""
    if threshold[0].isalnum():
        mutsig = pd.read_table(mutsig_file, sep="\t")
        # .replace("<", "")
        mutsig.loc[:, threshold[0]] = [x
                                       for x in mutsig.loc[:, threshold[0]]]
        return(list(mutsig[mutsig.loc[:, threshold[0]].astype(float) <
                           float(threshold[1:])].loc[:, "gene"]))
    else:
        print("the threshold should be p + value or q + value")
        sys.exit()


def maf2snv(maf_file, genes=None):
    """Maf to snv."""
    # file with different encoding
    # firstly, read to byte
    with open(maf_file, "rb") as f:
        maf = [x.strip().split(b"\t") for x in f.readlines()]
    # secondly, read to dataframe
    pmaf = pd.DataFrame(maf)
    # rename colnames by the first line, then del it
    pmaf.columns = [x.decode() for x in pmaf.ix[0, :]]
    pmaf.drop(0, inplace=True)
    tmaf = pmaf.loc[
        :, ['Hugo_Symbol', 'Tumor_Sample_Barcode']].drop_duplicates()
    tmaf = tmaf.applymap(lambda x: x.decode())
    if genes is not None:
        tmaf = tmaf[tmaf.loc[:, "Hugo_Symbol"].isin(genes)]
    tmaf.loc[:,
             "Tumor_Sample_Barcode"] = [x[:12] for x in
                                        tmaf.loc[:, "Tumor_Sample_Barcode"]]
    samples = set(tmaf.loc[:, "Tumor_Sample_Barcode"])
    snv_mat = [[x] for x in samples]
    for x in snv_mat:
        x.extend(list(tmaf[
            tmaf.loc[:,
                     "Tumor_Sample_Barcode"] == x[0]].loc[:, "Hugo_Symbol"]))
    return(samples, snv_mat)


def hotnet2(args, influence_matrices, snv_file):
    """Run hotnet2."""
    cmd = "cat " + snv_file
    if args.cna_file:
        cmd = cmd + " " + args.cna_file
    sample_cmd = cmd + "|cut -f 1 |sort -u > " + \
        args.output_dir + "/sample.txt"
    gene_cmd = cmd + \
        "|cut -f 2- |sed 's/\\t/\\n/g;s/(A)//g;s/(D)//g'| sort -u > " + \
        args.output_dir + "/gene.txt"

    print(sample_cmd)
    os.system(sample_cmd)
    print(gene_cmd)
    os.system(gene_cmd)
    generateheat_cmd = "python " + args.hotnet2 + \
        "/generateHeat.py mutation --snv_file " + snv_file
    if args.cna_file:
        generateheat_cmd = generateheat_cmd + " --cna_file " + args.cna_file
    generateheat_cmd = generateheat_cmd + " --sample_file " + \
        args.output_dir + "/sample.txt --gene_file " + \
        args.output_dir + "/gene.txt --output_file " + \
        args.output_dir + "/heat.json"
    print(generateheat_cmd)
    os.system(generateheat_cmd)

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
        print(simple_cmd)
        os.system(simple_cmd)
        results_net = glob.glob(
            args.output_dir + "/" + net + "/*/results.json")
        makeresultswebsite_cmd = "python " + args.hotnet2 + \
            "/bin/makeResultsWebsite.py --results_files " + \
            " ".join(results_net) +\
            " --edge_file " + influence_matrices[net]["edge_list"] +\
            " --network_name " + net +\
            " --output_directory " + args.output_dir + "/" + net + "/viz"
        print(makeresultswebsite_cmd)
        os.system(makeresultswebsite_cmd)


def visnetwork(args, snv_file):
    """Visnetwork.

    To filter the sig genes and network from hotnet2 output dir and
    creat mat.txt and heat.txt files.
    """
    file_sigs = []
    for root, _, files in os.walk(args.output_dir):
        for file in files:
            if file.endswith("heat.json"):
                heat_file = os.path.join(root, file)
            if file.endswith("subnetworks.json"):
                file_sigs.append((os.path.join(root, file)))
    print(heat_file)
    print(file_sigs)
    if len(file_sigs) == 0:
        print("the hotnet2 output_dir is NULL ")
        sys.exit()

    mat = pd.DataFrame()
    for f in file_sigs:
        a = subnetworkmat(f)
        mat = pd.concat([mat, a])
    mat.to_csv(os.path.join(args.output_dir, "mat.txt"),
               sep="\t", header=True, index=False)

    heatfile = json.load(open(heat_file))
    heatfile = pd.DataFrame.from_dict([heatfile['heat']])
    genes = set(mat.loc[:, "source"].tolist() + mat.loc[:, "target"].tolist())
    heatfile = heatfile.loc[:, genes]
    heatfile = heatfile.T
    heatfile.columns = ["heatscore"]

    if snv_file is not None and os.path.exists(snv_file):
        snv_genes = []
        with open(snv_file) as f:
            for line in f:
                word = line.strip().split()
                snv_genes.extend(word[1:])
        snv_genes = set(snv_genes)
        heatfile.loc[
            set(heatfile.index.intersection(snv_genes)), "type"] = "SNV"

    if args.cna_file is not None and os.path.exists(args.cna_file):
        cna_genes = []
        with open(args.cna_file) as f:
            for line in f:
                word = line.strip().split()
                word1 = [x.split("(")[0] for x in word[1:]]
                cna_genes.extend(word1)
        cna_genes = set(cna_genes)
        heatfile.loc[
            set(heatfile.index.intersection(cna_genes)), "type"] = "SCNA"
        if len(snv_genes) > 0 and len(cna_genes) > 0:
            heatfile.loc[snv_genes.intersection(cna_genes), "type"] = "BOTH"

    heatfile.index.name = "gene"
    heatfile.to_csv(os.path.join(args.output_dir, "heat.txt"),
                    sep="\t", header=True, index=True)


def subnetworkmat(file):
    """Read subnetwork.json.

    Return: subnetwork deltas k p from to
    """
    data = json.load(open(file))
    deltas = sorted(data['deltas'])

    for delta in deltas:
        p_mat = pd.DataFrame.from_dict(data['stats'][delta])
        k = p_mat.columns[p_mat.loc["pval"] < 0.05]
        if len(k) > 0:
            k = min([int(x) for x in k])
            pval = p_mat.loc['pval', str(k)]
            mat = [x for x in data['subnetworks'][delta]
                   if len(x['nodes']) >= int(k)]
            mat_out = pd.DataFrame()
            for m_value in mat:
                m_edges = pd.DataFrame.from_dict(m_value['edges'])
                mat_out = pd.concat([mat_out, m_edges])
            mat_out.loc[:, "networks"] = [x[0]
                                          for x in mat_out.loc[:, "networks"]]
            mat_out.loc[:, "deltas"] = delta
            mat_out.loc[:, "ksize"] = k
            mat_out.loc[:, "pval"] = pval
            return mat_out
            break


if __name__ == '__main__':
    sys.exit(main())
