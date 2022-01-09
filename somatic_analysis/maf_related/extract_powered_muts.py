#!/usr/bin/env python

import sys
import os.path
from argparse import ArgumentParser


def main(args):
    fo = None
    if args.output_blacklist:
        fo = open(os.path.basename(args.maf).split(".maf")[0] + ".not_powered.txt", "w")

    HEADER = ["Chromosome", "Start_Position", "End_Position",
              "Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode",
              args.power_column]
    locus_index = []
    for line in open(args.maf):
        line = line.rstrip("\n").split("\t")
        if line[0].startswith("Hugo"):
            index = line.index(args.power_column)
            locus_index = [line.index(h) for h in HEADER]
            print("\t".join(line))
        else:
            try:
                power = float(line[index])
            except ValueError:
                power = 0
            if power >= args.power_cutoff:
                print("\t".join(line))
            else:
                fo.write("\t".join([line[i] for i in locus_index])+"\n")
    fo.close()

               

if __name__ == "__main__":
    description = "Extract powered mutations from the list"
    parser = ArgumentParser(description=description)
    parser.add_argument("--power-cutoff", type=float, default=0.95, help="Cutoff of detection power (default: 0.95)")
    parser.add_argument("--output-blacklist", action="store_true", help="Generate a list of mutations that are not powered")
    parser.add_argument("power_column", help="Column that has power to check (e.g. validation_power_wex, validation_power_targeted)")
    parser.add_argument("maf", help="Validated MAF")
    args = parser.parse_args()
    try:
        main(args)
    except KeyboardInterrupt: pass

