#!/bin/usr/env python

import sys
import os.path
from argparse import ArgumentParser


def load_maf(maf):
    index = []
    header = ['Chromosome', 'Start_position', 'End_position',
              'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele2',
              'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode']
    maf_dict = {}
    maf_header = []
    for line in open(maf):
        if line.startswith("#"):
            continue
        line = line.rstrip("\n").split("\t")
        if line[0].startswith("Hugo"):
            index = [header.index(ll) for ll in line if ll in header]
            maf_header = line[:]
        else:
            key = "::".join([line[i] for i in index])
            maf_dict.setdefault(key, line)
    return maf_header, maf_dict


def write_xor_maf(keys, maf_dict, header, output_file):
    fo = open(output_file, "w")
    fo.write("\t".join(header)+"\n")
    for key in keys:
        fo.write("\t".join(maf_dict[key])+"\n")
    fo.close()


def main(args):
    maf1_header, maf1_dict = load_maf(args.maf_one)
    maf2_header, maf2_dict = load_maf(args.maf_two)
    
    maf1_keys = set(maf1_dict.keys())
    maf2_keys = set(maf2_dict.keys())

    maf1_output = os.path.basename(args.maf_one).rstrip(".maf") + "_Only.maf"
    maf2_output = os.path.basename(args.maf_two).rstrip(".maf") + "_Only.maf"

    write_xor_maf(list(maf1_keys.difference(maf2_keys)), maf1_dict, maf1_header, maf1_output)
    write_xor_maf(list(maf2_keys.difference(maf1_keys)), maf2_dict, maf2_header, maf2_output)


if __name__ == "__main__":
    parser = ArgumentParser(description="Print out differences in two MAFs as separate MAFs")
    parser.add_argument("maf_one", help="first maf")
    parser.add_argument("maf_two", help="second maf")
    args = parser.parse_args()
    try:
        main(args)
    except KeyboardInterrupt: pass

