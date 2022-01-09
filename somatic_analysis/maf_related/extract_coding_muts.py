#!/usr/bin/env python

import sys
import os.path
from argparse import ArgumentParser

# coding mutations
CODING_MUT_CLASS = [
          'Frame_Shift_Del',
          'Frame_Shift_Ins',
          'In_Frame_Del',
          'In_Frame_Ins',
          'Missense_Mutation',
          'Nonsense_Mutation',
          'Splice_Site',
          'Nonstop_Mutation']

HEADER = "Variant_Classification"


def main(args):
    for line in open(args.maf):
        line = line.rstrip("\n")
        if line.startswith("#"):
            print(line)
            continue
        line = line.split("\t")
        if line[0].startswith("Hugo"):
            variant_class_index = line.index(HEADER)
            print("\t".join(line))
        else:
            if line[variant_class_index] in CODING_MUT_CLASS:
                print("\t".join(line))


if __name__ == "__main__":
    parser = ArgumentParser(description="Extract coding mutations from MAF")
    parser.add_argument("maf", help="Input MAF")
    args = parser.parse_args()
    try:
        main(args)
    except KeyboardInterrupt: pass

