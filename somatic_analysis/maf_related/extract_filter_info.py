#!/usr/bin/env python

import sys
from argparse import ArgumentParser

def main(args):
    chrom_index = -1
    beg_index = -1
    end_index = -1
    mut_index = -1
    filter_index = -1
    alt_index = -1
    ref_index = -1

    white_list = set()
    if args.reference_maf:
        for line in open(args.reference_maf):
            if line.startswith("#"):
                continue
            line = line.rstrip("\n").split("\t")
            if line[0].startswith("Hugo"):
                chrom_index = line.index("Chromosome")
                beg_index = line.index("Start_position")
                end_index = line.index("End_position")
                mut_index = line.index("Variant_Type")
                continue
            else:
                white_list.add(":".join([line[chrom_index], line[beg_index], line[end_index], line[mut_index]]))

    print("\t".join(["Locus", "Filter", "AF"]))
    for line in open(args.union_maf):
        if line.startswith("#"):
            continue
        line = line.rstrip("\n").split("\t")
        if line[0].startswith("Hugo"):
            chrom_index = line.index("Chromosome")
            beg_index = line.index("Start_position")
            end_index = line.index("End_position")
            mut_index = line.index("Variant_Type")
            filter_index = line.index("Failed_Filters")
            alt_index = line.index("t_alt_count")
            ref_index = line.index("t_ref_count")
        else:
            locus = ":".join([line[chrom_index], line[beg_index], line[end_index], line[mut_index]])
            if len(white_list) > 0 and locus not in white_list:
                continue
            af = float(line[alt_index])/(float(line[alt_index]) + float(line[ref_index]))
            if not line[filter_index]:
                print("\t".join([locus, "PASS", str(af)]))
            else:
                for f in line[filter_index].rstrip(",").split(","):
                    print("\t".join([locus, f, str(af)]))


if __name__ == "__main__":
    parser = ArgumentParser(description="List up mutation filters to check the contribution")
    parser.add_argument("--reference-maf", help="MAF to narrow down area to investigate")
    parser.add_argument("union_maf", help="MAF to check filters")
    args = parser.parse_args()
    try:
        main(args)
    except KeyboardInterrupt: pass

