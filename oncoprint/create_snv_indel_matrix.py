#!/usr/bin/env python

import sys
import os.path
import argparse

# required header in a maf
HEADER = ['Hugo_Symbol',
          'Reference_Allele',
          'Tumor_Seq_Allele2',
          'Variant_Type',
          'Variant_Classification',
          'Protein_Change',
          'Tumor_Sample_Barcode']
          #'ccf_hat']

# coding mutations
CODING_MUT_CLASS = {
          'Frame_Shift_Del': "INFRAME",
          'Frame_Shift_Ins': "INFRAME",
          'In_Frame_Del': "INFRAME",
          'In_Frame_Ins': "INFRAME",
          'Missense_Mutation': "MISSENSE",
          'Nonsense_Mutation': "TRUNC",
          'Splice_Site': "OTHER",
          'Nonstop_Mutation': "OTHER" }
MUT_CATEG = {
          'Frame_Shift_*': "frame_shift",
          'In_Frame_*': "in_frame_indel",
          'Missense_Mutation': "missense",
          'Nonsense_Mutation': "nonsense",
          'Splice_Site*': "splice_site",
          'Nonstop_Mutation': "other_non_syn",
          'De_novo_Start*': "other_non_syn",
          'Start_*': "other_non_syn",
          'Translation_*': "other_non_syn",
          'Read-through*': "other_non_syn"}


def translate_muts(ref, alt):
    return FOLDED_CONTEXT[ ref + ">" + alt ]

def main(args):
    gene_list = [line.rstrip() for line in open(args.gene_list)]
    sample_names = []

    index = []
    mut_matrix = dict()

    for x in open(args.maf):
        if x.startswith("#"):
            continue
        x = x.split("\t")
        # read header in maf
        if x[0].startswith("Hugo"):
            column_num = len(x)
            try:
                index = [x.index(colname) for colname in HEADER]
            except ValueError as e:
                raise Exception(e.message.split()[0] + " is missing in MAF")
            continue
        # skip blank or aberrent lines
        if len(x) != column_num:
            sys.stderr.write("skipped: " + "\t".join(x) + "\n")
            continue

        sample_name = x[index[HEADER.index('Tumor_Sample_Barcode')]]
        if sample_name not in sample_names:
            sample_names.append(sample_name)
            mut_matrix.setdefault(sample_name, {g:set() for g in gene_list})

        gene = x[index[HEADER.index('Hugo_Symbol')]]
        if gene not in gene_list:
            continue

        added = False
        mut_class = x[index[HEADER.index('Variant_Classification')]]
        for categ in MUT_CATEG:
            if categ == mut_class:
                mut_matrix[sample_name][gene].add(MUT_CATEG[categ])
                added = True
            elif categ.endswith("*"):
                if mut_class.startswith(categ.replace("*","")):
                    mut_matrix[sample_name][gene].add(MUT_CATEG[categ])
                    added = True

        # add clonality
        #ccf = x[index[HEADER.index('ccf_hat')]]
        #if ccf == "NA" or added == False:
        #    continue
        #if float(ccf) >= 0.9:
        #    mut_matrix[sample_name][gene].add("clonal")
        #else:
        #    if "clonal" not in mut_matrix[sample_name][gene]:
        #        mut_matrix[sample_name][gene].add("subclonal")

    print("\t".join(["gene_name"] + sample_names))
    for gene in gene_list:
        line = [gene]
        for sample_name in sample_names:
            flag = ";".join(sorted(mut_matrix[sample_name][gene]))
            if len(flag) == 0:
                flag = "NA"
            else:
                flag += ";"
            line.append(flag)
        print("\t".join(line))



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("gene_list", help="Gene list")
    parser.add_argument("maf", help="input MAF")
    
    args = parser.parse_args()
    try:
        main(args)
    except KeyboardInterrupt: pass
