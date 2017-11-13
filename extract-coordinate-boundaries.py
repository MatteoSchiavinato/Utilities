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
August 24, 2017 

--------------------------------------------------------------------------------------------------------

[ input options ]
-g	--gff		Input GFF/GTF/GFF3 file (only fields 1,4,5 count)	[stdin]

[ output options ]
-o	--output-file	Output file with <TAB> separated Scaffold, start, stop	[stdout]
			representing boundaries of coverage. 
-e	--error-file	Standard error file (log), requires --verbose		[stderr]
-v	--verbose	Turn on verbose standard error 				[off]
''')


# parser
p = ap.ArgumentParser()
p.add_argument("-g", "--gff")
p.add_argument("-o", "--output-file")
p.add_argument("-e", "--error-file")
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

if args.error_file:
	ERROR = open(args.error_file, "w")
else:
	ERROR = sys.stderr


if args.verbose:
	ERROR.write("[{0}] Reading input file, computing boundaries ... \n".format(at()))

dic = {}
for line in INPUT:
	lst = line.rstrip("\r\n\b").split("\t")
	scaffold = lst[0]
	position = int(lst[3])
	#
	try:
		dic[scaffold]
	except KeyError:
		dic[scaffold] = [[position,position]]
	#
	k=0
	lastIndex = len(dic[scaffold])-1
	while ((position > dic[scaffold][k][1]+1) and (k < lastIndex)):
		k+=1
	if (position == dic[scaffold][k][1]+1):
		dic[scaffold][k][1] = position 
	elif (position > dic[scaffold][k][1]+1):
		dic[scaffold].append([position, position])

for scaffold in sorted(dic.keys()):
	for region in dic[scaffold]:
		OUTPUT.write("\t".join([scaffold, str(region[0]), str(region[1])]) + "\n")

INPUT.close()
OUTPUT.close()

if args.verbose:
	ERROR.write("[{0}] Done.\n".format(at()))

ERROR.close()
