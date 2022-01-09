#!/usr/bin/env python

import sys

index = []
for line in open(sys.argv[1]):
    line = line.rstrip("\n").split("\t")
    if line[0].startswith("Hugo"):
        index = [i for i,x in enumerate(line) if "_pon_" not in x]
    print("\t".join([line[i] for i in index]))

