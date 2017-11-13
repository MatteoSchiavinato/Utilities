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
September 22, 2017 

--------------------------------------------------------------------------------------------------------

[ input options ]
-a		--gff-a			Input gff file A					[-]
-b		--gff-b			Input gff file B					[-]

[ filtering options ]
-f 		--feature		Feature to be considered in the third field		[exon]
	
[ output options ]
-o		--output-file		Output overlap file					[stdout]
-e		--error-file		Standard error file (log)				[stderr]
-v		--verbose		Turn on verbose output					[off]

''')


# parser
p = ap.ArgumentParser()
p.add_argument("-a", "--gff-a")
p.add_argument("-b", "--gff-b")
p.add_argument("-f", "--feature", default="exon")
p.add_argument("-o", "--output-file")
p.add_argument("-v", "--verbose", action="store_true")
p.add_argument("-e", "--error-file")
args = p.parse_args()

if args.output_file:
	OUTPUT = open(args.output_file, "w")
else:
	OUTPUT = sys.stdout

if args.error_file:
	ERROR = open(args.error_file, "w")
else:
	ERROR = sys.stderr


# reading file A
if args.verbose:
	ERROR.write("[{0}] Reading File A\n".format(at()))

FileA = {}
INPUT = open(args.gff_a, "r")
for line in INPUT:
	lst = line.rstrip("\b\r\n").split("\t")
	if lst[2] == args.feature:
		scaffold = str(lst[0])
		start = int(lst[3])
		end = int(lst[4])
		try:
			if len(FileA[scaffold]) > 0:
				for region in FileA[scaffold]:
					if end == region[0]-1:
						idx = FileA[scaffold].index(region)
						FileA[scaffold][idx] = [start,region[1]]
					elif start == region[1]+1:
						idx = FileA[scaffold].index(region)
						FileA[scaffold][idx] = [region[0],end]
					elif region[0] <= end <= region[1]:
						idx = FileA[scaffold].index(region)
						FileA[scaffold][idx] = [start,region[1]]
					elif region[0] <= start <= region[1]:
						idx = FileA[scaffold].index(region)
						FileA[scaffold][idx] = [region[0],end]
					else:
						FileA[scaffold].append([start,end])
			else:
				FileA[scaffold] = [[start, end]]
		except KeyError:
			FileA[scaffold] = [[start, end]]
INPUT.close()


# reading file B
if args.verbose:
	ERROR.write("[{0}] Reading file B\n".format(at()))

FileB = {}
INPUT = open(args.gff_b, "r")
for line in INPUT:
	lst = line.rstrip("\b\r\n").split("\t")
	if lst[2] == args.feature:
		scaffold = str(lst[0])
		start = int(lst[3])
		end = int(lst[4])
		try:
			if len(FileB[scaffold]) > 0:
				for region in FileB[scaffold]:
					if end == region[0]-1:
						idx = FileB[scaffold].index(region)
						FileB[scaffold][idx] = [start,region[1]]
					elif start == region[1]+1:
						idx = FileB[scaffold].index(region)
						FileB[scaffold][idx] = [region[0],end]
					elif region[0] <= end <= region[1]:
						idx = FileB[scaffold].index(region)
						FileB[scaffold][idx] = [start,region[1]]
					elif region[0] <= start <= region[1]:
						idx = FileB[scaffold].index(region)
						FileB[scaffold][idx] = [region[0],end]
					else:
						FileB[scaffold].append([start,end])
			else:
				FileB[scaffold] = [[start, end]]
		except KeyError:
			FileB[scaffold] = [[start, end]]
INPUT.close()


# calculate overlaps 
if args.verbose:
	ERROR.write("[{0}] Computing overlaps\n".format(at()))

for scaffold in FileA:
	for x in FileA[scaffold]:
		start = x[0]
		end = x[1]
		tot_overlap = 0
		if scaffold in FileB.keys():
			for y in FileB[scaffold]:
				if y[0] <= start <= y[1]:
					overlap = y[1]-start+1
					tot_overlap += overlap
					OUTPUT.write("\t".join([scaffold, str(x[0]), str(x[1]), scaffold, str(y[0]), str(y[1]), str(overlap)]) + "\n")
				elif y[0] <= end <= y[1]:
					overlap = end-y[0]+1
					tot_overlap += overlap 
					OUTPUT.write("\t".join([scaffold, str(x[0]), str(x[1]), scaffold, str(y[0]), str(y[1]), str(overlap)]) + "\n")
				elif y[0] > end:
					break

OUTPUT.close()


# verbose output 
if args.verbose:
	ERROR.write("[{0}] Total overlap: {1}\n".format(at(), str(tot_overlap)))

ERROR.close()
