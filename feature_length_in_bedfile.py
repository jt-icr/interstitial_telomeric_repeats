#!/usr/bin/env python3.5
# Finds the motif lengths of features in bed file
# JP Tomkins <jtomkins@icr.org>

import sys

args = sys.argv[1:]
if len(args) != 1 :
    print("Error! Argument needed: input file name")
    sys.exit()
fn = args[0]

L = []

# Put bed file line items into lists in L
with open(fn, 'r') as f:
    for  line in f:
        L.append(line.strip().split())

#Remove empty lists and all that's falsey
L = [x for x in L if x]

# Print out the chr and num bases of motif
for x in L:
    print(x[0] + ": ", int(x[2])-int(x[1]))