#!/usr/bin/env python

import sys
import os.path
from argparse import ArgumentParser


def main(args):
    samples = [os.path.basename(f).replace(".segtab.intersected.txt", "") for f in args.intersected_bed]
    genes = [line.rstrip("\n") for line in open(args.gene_list)]
    ploidy_matrix = dict()
    for gene in genes:
        ploidy_matrix.setdefault(gene, dict())
        for sample in samples:
            ploidy_matrix[gene].setdefault(sample, [[], []])
    
    for i, f in enumerate(args.intersected_bed):
        sample = samples[i]
        for line in open(f):
            line = line.rstrip("\n").split("\t")
            gene = line[3]
            allele_1 = -1
            allele_2 = -1
            overlap_len = -1
            try:
                allele_1 = int(line[30])
                allele_2 = int(line[31])
                overlap_len = int(line[49])
                ploidy_matrix[gene][sample][0].append((allele_1, allele_2))
                ploidy_matrix[gene][sample][1].append(overlap_len)
            except ValueError:
                pass

    # AMP >= 4, GAIN = 3, LOSS = 1, DEL = 0
    # NEU-LOH = allele_{1,2} == (2,0) && cn == 2
    fo_matrix = open(args.prefix+".abs_cn_matrix.txt", "w")
    fo_onco = open(args.prefix+".cnv_oncoprint_matrix.txt", "w")
    fo_matrix.write("\t".join(["gene_name"] + samples)+"\n")
    fo_onco.write("\t".join(["gene_name"] + samples)+"\n")
    for gene in genes:
        matrix_line = []
        onco_line = []
        for sample in samples:
            if len(ploidy_matrix[gene][sample][0]) == 0:
                matrix_line.append("NA")
                onco_line.append("NA")
            else:
                allele_1 = 0
                allele_2 = 0
                total_len = sum(ploidy_matrix[gene][sample][1])
                for i in range(len(ploidy_matrix[gene][sample][0])):
                    allele_1 += ploidy_matrix[gene][sample][0][i][0] * ploidy_matrix[gene][sample][1][i]/total_len
                    allele_2 += ploidy_matrix[gene][sample][0][i][1] * ploidy_matrix[gene][sample][1][i]/total_len
                copy_num = round(allele_1 + allele_2)
                label = "NA"
                if copy_num >= 4:
                    label = "amplification"
                elif copy_num == 3:
                    label = "gain"
                elif copy_num == 1:
                    label = "loss"
                elif copy_num == 0:
                    label = "deletion"
                elif copy_num == 2 and (allele_1 == 0 or allele_2 == 0):
                    label = "nLOH"
                matrix_line.append(str(copy_num))
                onco_line.append(label)
        fo_matrix.write("\t".join([gene] + matrix_line)+"\n")
        fo_onco.write("\t".join([gene] + onco_line)+"\n")
    fo_matrix.close()
    fo_onco.close()


if __name__ == "__main__":
    description = "Create gene matrix from ABSOLUTE intersected file"
    parser = ArgumentParser(description=description)
    parser.add_argument("prefix")
    parser.add_argument("gene_list")
    parser.add_argument("intersected_bed", nargs="+", help="ABSOLUTE seg with GENCODE annotation")
    args = parser.parse_args()
    try:
        main(args)
    except KeyboardInterrupt: pass
