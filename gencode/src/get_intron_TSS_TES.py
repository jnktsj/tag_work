#!/usr/bin/env python

import sys
import fileinput
from argparse import ArgumentParser


def main(args):
    chrom_len = {}
    for line in open(args.chrom_length):
        line = line.rstrip().split("\t")
        chrom_len.setdefault(line[0], int(line[1]))
    exons = {}
    for line in fileinput.input(args.exon_bed):
        line = line.rstrip().split("\t")
        key = line[3] + ":" + line[0]
        exons.setdefault(key, [])
        exons[key].append(line)

    # print up-/down-stream regions and introns
    fr = open("TSS_TES.bed", "w")
    fi = open("introns.bed", "w")
    for gene_name in exons:
        has_exon = len(exons[gene_name]) > 1
        for i, exon in enumerate(exons[gene_name]):
            chrom = exon[0]
            beg = int(exon[1])
            end = int(exon[2])
            strand = gene_name.split(":")[-2]
            gene = ":".join(gene_name.split(":")[:-1])
            if i == 0 and strand == "+":
                fr.write("\t".join([chrom, str(max(beg-1000, 0)), str(beg), gene+":TSS"]) + "\n")
            if i == 0 and strand == "-":
                fr.write("\t".join([chrom, str(max(beg-1000, 0)), str(beg), gene+":TES"]) + "\n")
            if i == len(exons[gene_name]) - 1 and strand == "+":
                fr.write("\t".join([chrom, str(end), str(min(end+1000, chrom_len[chrom])), gene+":TES"]) + "\n")
            if i == len(exons[gene_name]) - 1 and strand == "-":
                fr.write("\t".join([chrom, str(end), str(min(end+1000, chrom_len[chrom])), gene+":TSS"]) + "\n")
            if (i < len(exons[gene_name]) - 1) and has_exon:
                fi.write("\t".join([chrom, str(end), exons[gene_name][i+1][1], gene+":Intron"]) + "\n")
    fr.close()
    fi.close()


if __name__ == "__main__":
    parser = ArgumentParser(description="Obtain introns and up/down-stream regions of a gene")
    parser.add_argument("chrom_length", help="Chromosome length file")
    parser.add_argument("exon_bed", help="Exon BED (needs to be sorted!)")
    args = parser.parse_args()
    try:
        main(args)
    except KeyboardInterrupt: pass

    
