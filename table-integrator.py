#!/usr/bin/env python2.7

import argparse as ap
import sys 

# help
if len(sys.argv) < 2:
	sys.argv.append("-h")
if sys.argv[1] in ["-h", "-help", "--help", "usage", "getopt"]:
	sys.exit('''

---
Matteo Schiavinato
December 04, 2017
BOKU - Vienna University of Natural Resources and Life Sciences
---

table-integrator.py [-f|--files] > stdout

-f	--files			Specify as many files as you want, separated by space 
				Files need all to have the following format:
					- header with one name per line including first one 
					- row names on the first column 
				The name of the first column will be discarded (implied: row names) 
				It doesn't matter if a file doesn't contain all the names of another file
				The script will put "NA" where an empty field is found 

-n	--first-line-name	Default="Name"

-m	--missing-values	Default="NA"

''')

p = ap.ArgumentParser()
p.add_argument("-f", "--files", nargs="*", required=True)
p.add_argument("-n", "--first-line-name", default="Name")
p.add_argument("-m", "--missing-values", default="NA")
args = p.parse_args()

header = [args.first_line_name]
Lines = {}

for filename in args.files:
	INPUT = open(filename, "r")
	k=0
	for line in INPUT:
		lst = line.rstrip("\b\r\n").split("\t")
		if k==0:
			invoices = lst[1:]
			for x in invoices:
				header.append(x)
			k+=1
		elif k>0:
			name = lst[0]
			values = lst[1:]
			try:
				Lines[name]
			except KeyError:
				Lines[name] = {}
			#
			for x in invoices:
				idx = invoices.index(x)
				value = values[idx]
				Lines[name][x] = value
			#
			k+=1
	INPUT.close()

for key in Lines:
	for x in header:
		if x not in Lines[key].keys():	
			Lines[key][x] = args.missing_values

sys.stdout.write("\t".join(header) + "\n")
for key in Lines:
	value = []
	for x in header[1:]:
		value.append(Lines[key][x])
	valueFinal = "\t".join(value)
	sys.stdout.write("\t".join([key, valueFinal]) + "\n")
