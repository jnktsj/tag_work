#!/usr/bin/env python

import sys
from argparse import ArgumentParser

prefixes = ('HGVS_', 'CGC_', 'ESP_', 'OREGANNO_', 'UniProt_', 'COSMIC_',
            'TCGAscape_', 'v_1000gp3_', 'x1000gp3_', 'dbNSFP_', 'HGNC_', 'X1000gp3_')

index = []
for line in open(sys.argv[1]):
    line = line.rstrip("\n").split("\t")
    if line[0].startswith("Hugo"):
        for i, x in enumerate(line):
            include = True
            for p in prefixes:
                if x.startswith(p):
                    include = False
                    break
            if include == True:
                index.append(i)
    print("\t".join([line[i] for i in index]))
