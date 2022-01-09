#!/usr/bin/env python

import sys
from argparse import ArgumentParser


REQ_COL = ["Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Type", "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode"]
#REQ_COL = ["Chromosome", "Start_position", "End_position", "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Type", "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode"]
ref_name = "t_ref_count"
alt_name = "t_alt_count"


def read_tsv(input_file, first_header_elem):
    index = []
    ref_index = -1
    alt_index = -1
    allele_frac = {}
    loci = {}
    for line in open(input_file):
        line = line.rstrip().split("\t")
        if line[0] == first_header_elem:
            index = [line.index(rcol) for rcol in REQ_COL]
            ref_index = line.index(ref_name)
            alt_index = line.index(alt_name)
            continue
        data = [line[i] for i in index]
        af = float(line[alt_index])/(float(line[alt_index]) + float(line[ref_index]))
        allele_frac.setdefault(":".join(data), [str(af), line[alt_index], line[ref_index]])
        sample_name = data[-2] + ":" + data[-1]
        loci.setdefault(sample_name, set()).add(tuple(data[:-2]))
    return allele_frac, loci


def main(args):
    maf1_af, maf1_loci = read_tsv(args.MAF1, "Hugo_Symbol")
    maf2_af, maf2_loci = read_tsv(args.MAF2, "Hugo_Symbol")

    f = open(args.prefix + "-overlap.txt", "w")
    f.write("\t".join(["SAMPLE", "MAF1_ONLY", "OVERLAP", "MAF2_ONLY"]) + "\n")
    for sample_name in maf1_loci:
        common = len(maf1_loci[sample_name].intersection(maf2_loci.get(sample_name, set())))
        maf1_only = len(maf1_loci[sample_name].difference(maf2_loci.get(sample_name, set())))
        maf2_only = len(maf2_loci[sample_name].difference(maf1_loci.get(sample_name, set())))
        f.write("\t".join([sample_name.split(":")[0], str(maf1_only), str(common), str(maf2_only)]) + "\n")
    f.close()

    f = open(args.prefix + "-vaf-count.txt", "w")
    f1 = open(args.prefix + "-maf1-only.txt", "w")
    f.write("\t".join(["sample", "chrom", "start", "end", "ref", "alt", "mut_type",
                       "maf1_vaf", "maf1_alt_count", "maf1_ref_count",
                       "maf2_vaf", "maf2_alt_count", "maf2_ref_count"]) + "\n")
                       
    f1.write("\t".join(["sample", "chrom", "start", "end", "ref", "alt", "mut_type",
                        "maf1_vaf", "maf1_alt_count", "maf1_ref_count"]) + "\n")
    for key in maf1_af:
        elems = key.split(":")[:-2]
        sample_name = key.split(":")[-2]
        if key not in maf2_af:
            f1.write("\t".join([sample_name] + elems + maf1_af[key]) + "\n")
            continue
        f.write("\t".join([sample_name] + elems + maf1_af[key] + maf2_af[key]) + "\n")
    f.close()
    f1.close()

    f2 = open(args.prefix + "-maf2-only.txt", "w")
    f2.write("\t".join(["sample", "chrom", "start", "end", "ref", "alt", "mut_type",
                        "maf2_vaf", "maf2_alt_count", "maf2_ref_count"]) + "\n")
    for key in maf2_af:
        elems = key.split(":")[:-2]
        sample_name = key.split(":")[-2]
        if key not in maf1_af:
            f2.write("\t".join([sample_name] + elems + maf2_af[key]) + "\n")
    f2.close()


if __name__ == "__main__":
    parser = ArgumentParser(description="Output concordance of two MAFs")
    parser.add_argument("--prefix", help="output prefix", default="out")
    parser.add_argument("MAF1", help="MAF input 1")
    parser.add_argument("MAF2", help="MAF input 2")
    args = parser.parse_args()
    try:
        main(args)
    except KeyboardInterrupt: pass

