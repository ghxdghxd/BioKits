#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# tableReport.py
# @Author : JT Guo (guojt-4451@163.com)
# @Date   : 2018/9/17 21:37:46

"""Anonotated genes using DAVID"""

import sys
import logging
import traceback as tb
import suds.metrics as metrics
from suds import *
from suds.client import Client
from datetime import datetime
import subprocess
import mygene

__version__ = "%(prog)s"

def setup_logging():
    if sys.version_info < (2, 5):
        fmt = '%(asctime)s [%(levelname)s] @%(filename)s:%(lineno)d\n%(message)s\n'
    else:
        fmt = '%(asctime)s [%(levelname)s] %(funcName)s() @%(filename)s:%(lineno)d\n%(message)s\n'
    logging.basicConfig(level=logging.INFO, format=fmt)

def get_args3():
    """Get arguments from commond line."""
    try:
        import argparse as ap
    except ImportError as imerr:
        print("\033[1;31m" + str(imerr) + " \033[0m")
        sys.exit()
    parser = ap.ArgumentParser(fromfile_prefix_chars='@',
                               description=__doc__,
                               formatter_class=ap.RawTextHelpFormatter)
    parser.add_argument('-v', '--version',
                        action='version',
                        version=__version__)
    parser.add_argument("--inputIds",
                        metavar="FILE",
                        dest="inputIds",
                        required=True,
                        help="exp. Gene List")
    parser.add_argument("--idType",
                        metavar="STR",
                        dest="idType",
                        required=True,
                        choices=["OFFICIAL_GENE_SYMBOL", "AFFYMETRIX_3PRIME_IVT_ID", "AFFYMETRIX_EXON_ID", "AGILENT_CHIP_ID", "AGILENT_ID",
                                 "AGILENT_OLIGO_ID", "APHIDBASE_ID", "BEEBASE_ID", "BEETLEBASE_ID", "BGD_ID", "CGNC_ID",
                                 "CRYPTODB_ID", "DICTYBASE_ID", "ENSEMBL_GENE_ID", "ENSEMBL_TRANSCRIPT_ID", "ENTREZ_GENE_ID",
                                 "FLYBASE_GENE_ID", "GENBANK_ACCESSION", "GENOMIC_GI_ACCESSION", "GENPEPT_ACCESSION",
                                 "LOCUS_TAG", "MGI_ID", "MIRBASE_ID", "MRNA_GI_ACCESSION", "NASONIABASE_ID", "PROTEIN_GI_ACCESSION",
                                 "PSEUDOCAP_ID", "REFSEQ_MRNA", "REFSEQ_PROTEIN", "RGD_ID", "SGD_ID", "TAIR_ID", "UNIGENE",
                                 "UNIPROT_ACCESSION", "UNIPROT_ID", "VECTORBASE_ID", "WORMBASE_GENE_ID", "XENBASE_ID", "ZFIN_ID"],
                        help="Select Identifier")
    parser.add_argument("--listName",
                        metavar="STR",
                        dest="listName",
                        default="make_up",
                        help="the name of input gene list as well as the output file")
    parser.add_argument("--email",
                        metavar="STR",
                        dest="email",
                        required=True,
                        help="David email")
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    else:
        args = parser.parse_args()
        return args


def run_david(inputIds, idType, email, listName, input_dict):

    setup_logging()
    logging.getLogger('suds.client').setLevel(logging.DEBUG)
    url = 'https://david.ncifcrf.gov/webservice/services/DAVIDWebService?wsdl'
    print('url=%s' % url)

    # create a service client using the wsdl.
    client = Client(url)
    client.wsdl.services[0].setlocation(
        'https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/')

    #authenticate user email
    client.service.authenticate(email)

    #add a list
    listType = 0
    client.service.addList(inputIds, idType, listName, listType)

    #print client.service.getDefaultCategoryNames()
    categories = ["BBID", "BIOCARTA", "COG_ONTOLOGY", "GOTERM_BP_FAT", "GOTERM_CC_FAT", "GOTERM_MF_FAT",
                  "INTERPRO", "KEGG_PATHWAY", "OMIM_DISEASE", "PIR_SUPERFAMILY", "SMART", "SP_PIR_KEYWORDS", "UP_SEQ_FEATURE"]
    client.service.setCategories(",".join(categories))

    tableReport = client.service.getTableReport()
    tableRow = len(tableReport)
    print('Total table records:', tableRow)

    res = ['ID\tGene Name\tSpecies\t' + '\t'.join(categories)]
    for tableRecord in tableReport:
        name = tableRecord.name
        species = tableRecord.species
        gene_id = input_dict[tableRecord.values[0].array[0]] if len(
            input_dict) > 0 else tableRecord.values[0].array[0]
        rowList = [gene_id, name, species]
        category_dict = dict.fromkeys(categories, "")
        for annotationRecord in tableRecord.annotationRecords:
            # category_dict[str(annotationRecord.category)] = ",".join([x.split("$")[1] for x in annotationRecord.terms])
            category_dict[str(annotationRecord.category)] = ",".join([x.split("$")[1].replace("\ufffd","") for x in annotationRecord.terms])

        res.append("\t".join(rowList + [category_dict[key] for key in categories]))
            # .encode('utf-8').replace("\xef\xbf\xbd","").decode())

    resF = listName + ".tableReport.txt"
    with open(resF, 'w') as fOut:
        fOut.write("\n".join(res))


def main():
    args = get_args3()
    with open(args.inputIds) as f:
        inputIds = [line.strip() for line in f.readlines()]
        input_dict = {}
    if args.idType == "OFFICIAL_GENE_SYMBOL":
        mg = mygene.MyGeneInfo()
        inputIds = mg.querymany(inputIds, scopes='symbol',
                                species='human', as_dataframe=True)
        inputIds = inputIds["entrezgene"].dropna()
        input_dict = {inputIds[x]: inputIds.index[x]
                      for x in range(len(inputIds))}
        idType = "ENTREZ_GENE_ID"
    else:
        idType = args.idType

    inputIds = ",".join(inputIds)

    run_david(inputIds, idType, args.email, args.listName, input_dict)

if __name__ == '__main__':
    sys.exit(main())
