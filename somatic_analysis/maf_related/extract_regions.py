#!/usr/bin/env python

import sys
import os.path
from argparse import ArgumentParser

def main(args):
    # load interval list
    regions = {}
    for line in open(args.interval_list):
        if line.startswith("@"):
            continue
        line = line.rstrip().split("\t")
        bed = [line[0], int(line[1]), int(line[2])]
        if bed[1] >= bed[2]:
            continue
        regions.setdefault(bed[0], []).append(bed)
    # extract MAF lines
    for line in open(args.maf):
        line = line.rstrip("\n")
        if line.startswith("#") or line.startswith("Hugo"):
            print(line)
            continue
        line = line.split("\t")
        # this script expect specific column orders
        chrom, beg, end = line[4], int(line[5]), int(line[6])
        if chrom == "M":
            chrom = "MT"
            line[4] = chrom
        # skip any chromosmes that are not in the interval list
        if line[4] not in regions:
            continue
        if line[9] == "DEL":
            beg -= 1

        for pos in regions[chrom]:            
            if (pos[1] <= beg and beg <= pos[2]) or (pos[1] <= end and end <= pos[2]):
                print("\t".join(line))
                break

if __name__ == "__main__":
    description = "Extract MAF entries that overlap with genomic regions in an interval list"
    parser = ArgumentParser(description=description)
    parser.add_argument("interval_list", help="interval list")
    parser.add_argument("maf", help="target maf")
    args = parser.parse_args()
    try:
        main(args)
    except KeyboardInterrupt: pass

