#!/bin/env/python

import sys
import os.path
from argparse import ArgumentParser

MAF_COLUMNS = "/broad/tag_working/jtsuji/utility/src/maf_related/maf_columns.txt"

def main(args):
    index = []
    columns = [line.rstrip("\n") for line in open(MAF_COLUMNS)]
    for line in open(args.maf):
        if line.startswith("#") or line.rstrip("\n") == "":
            continue

        line = line.rstrip("\n").split("\t")
        if line[0].startswith("Hugo"):
            index = [i for i,x in enumerate(line) if x in columns]
        
        print("\t".join([line[i] for i in index]))


if __name__ == "__main__":
    parser = ArgumentParser(description="Create public MAF with selected columns")
    parser.add_argument("maf", help="Input MAF")
    args = parser.parse_args()
    try:
        main(args)
    except KeyboardInterrupt: pass
