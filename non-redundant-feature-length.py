#!/usr/bin/env python2.7 

# modules
import sys
import re
import argparse as ap
from time import asctime as at


# help
if len(sys.argv) < 2:
	sys.argv.append("-h")
if sys.argv[1] in ["-h", "-help", "--help", "usage", "getopt"]:
	sys.exit('''

Matteo Schiavinato
September 11, 2017 

--------------------------------------------------------------------------------------------------------

[ input options ]
-g		--gff			Input gff file 						[stdin]

[ criteria ]
-f		--feature		Selected feature					[-]

[ output options ]
-o		--output-file		Output filtered gff file				[stdout]
-v		--verbose		Turn on verbose output					[off]

''')


# parser
p = ap.ArgumentParser()
p.add_argument("-g", "--gff")
p.add_argument("-f", "--feature", required=True, type=str)
p.add_argument("-o", "--output-file")
p.add_argument("-v", "--verbose", action="store_true")
args = p.parse_args()


if args.gff:
	INPUT = open(args.gff, "r")
else:
	INPUT = sys.stdin

if args.output_file:
	OUTPUT = open(args.output_file, "w")
else:
	OUTPUT = sys.stdout

ERROR = sys.stderr


# read input file 
Feature = {}
for line in INPUT:
	lst = line.rstrip("\b\n\r").split("\t")
	if lst[2] == args.feature:
		scaffold = str(lst[0])
		start = int(lst[3])
		end = int(lst[4])
		try:
			Feature[scaffold]
		except KeyError:
			Feature[scaffold] = []
		if len(Feature[scaffold]) > 0:
			for region in Feature[scaffold]:
				if region[0] <= end <= region[1]:
					idx = Feature[scaffold].index(region)
					Feature[scaffold][idx] = [min(start,region[0]), region[1]]
				elif region[0] <= start <= region[1]:
					idx = Feature[scaffold].index(region)
					Feature[scaffold][idx] = [region[0], max(end,region[1])]
				else:
					Feature[scaffold].append([start,end])
		else:
			Feature[scaffold].append([start, end])


# sort each scaffold list in crescent way 


# check for further overlaps
NewFeature = []
counter=1	# fake initialization to introduce the "while"
while (counter > 0):
	counter = 0	# initialize counter
	for scaffold in Feature:
		NewFeature[scaffold] = []
		maxidx = len(Feature[scaffold])-1
		k=0
		while (k < maxidx):
			start = Feature[scaffold][k][0]
			end = Feature[scaffold][k][1]
			nextStart = Feature[scaffold][k+1][0]
			nextEnd = Feature[scaffold][k+1][1]
			if nextStart <= end <= nextEnd:
				newRegion = [min(start, nextStart), max(end, nextEnd)]
				counter += 1
				k += 1	# add +1 to k to skip next element (already merged) 
			

			idx = Feature[scaffold].index(x)
			for region in Feature[scaffold][:idx]	# python list ranges are half-open
				if region[0] <= start <= region[1]:
					
				elif region[0] <= end <= region[1]:
			for region in Feature[scaffold][idx+1:]
			
