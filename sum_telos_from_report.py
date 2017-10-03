#!/usr/bin/env python3.5
# J.P. Tomkins <jtomkins@icr.org>

import sys

args = sys.argv[1:]
if len(args) != 1 :
    print("Error! Argument needed: input file name")
    sys.exit()
filename = args[0]
fi = open(filename, 'r')

exact_fhits_merlist = []
exact_rhits_merlist = []
inexact_fhits_merlist = []
inexact_rhits_merlist = []

for line in fi:
    if line.startswith('<'):
        # Reset all flags and variables
        mer_flag = "none"
        mer_num = 0
        
    elif line.startswith('exact_fhits'):
        mer_flag = "ef"

    elif line.startswith('exact_rhits'):
        mer_flag = "er"

    elif line.startswith('inexact_fhits'):
        mer_flag = "if"

    elif line.startswith('inexact_rhits'):
        mer_flag = "ir"

    else:
        if '_mer' in line:
            # grab the number from the end of the input line
            mer_num = line.strip().rsplit(' ')[2]
            # append to the appropriate list based on mer_flag that is current
            if mer_flag == "ef":
                exact_fhits_merlist.append(int(mer_num))
            elif mer_flag == "er":
                exact_rhits_merlist.append(int(mer_num))
            elif mer_flag == "if":
                inexact_fhits_merlist.append(int(mer_num))
            elif mer_flag == "ir":
                inexact_rhits_merlist.append(int(mer_num))

fi.close()

print("<Summary of telo hits for all chromosomes>")
print("exact_fhits: " + str(sum(exact_fhits_merlist)))
print("exact_rhits: " + str(sum(exact_rhits_merlist)))
print("inexact_fhits: " + str(sum(inexact_fhits_merlist)))
print("inexact_rhits: " + str(sum(inexact_rhits_merlist)))
