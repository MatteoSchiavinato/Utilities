#!/usr/bin/env python

# modules
import argparse as ap
from time import asctime as ts
import sys
def ssw(x):
	sys.stderr.write(x)

# help
if len(sys.argv) < 2:
	sys.argv.append("-h")
if sys.argv[1] in ["-h", "-help", "--help", "usage", "getopt"]:
	sys.exit('''

-in	--input-file		A PSL file									[stdin]
	--remove-header		Remove header from the file (program will not consider the first 5 lines) 	[off]
-out	--output-file		The output PSL file 								[mandatory]
	--dupl-file		A PSL file containing equally good hits that have not been selected.		[off]
				These are hits that were as good as the selected one, but a choice had to	
				be made										
	--no-insertions		Don't keep PSL hits that have insertions in the query or the target sequence	[off]
	--ins-file		Store hits with insertions in this PSL file 					[off]
-v	--verbose		Turn on verbose output (with standard error comments) 				[off]

''')


# parser
p = ap.ArgumentParser()
p.add_argument("-in", "--input-file")
p.add_argument("-out", "--output-file")
p.add_argument("--dupl-file")
p.add_argument("--no-insertions", action="store_true")
p.add_argument("--remove-header", action="store_true")
p.add_argument("--ins-file")
p.add_argument("-v", "--verbose", action="store_true")
args = p.parse_args()

#    0.  matches - Number of matching bases that aren't repeats.
#    1.  misMatches - Number of bases that don't match.
#    2.  repMatches - Number of matching bases that are part of repeats.
#    3.  nCount - Number of 'N' bases.
#    4.  qNumInsert - Number of inserts in query.
#    5.  qBaseInsert - Number of bases inserted into query.
#    6.  tNumInsert - Number of inserts in target.
#    7.  tBaseInsert - Number of bases inserted into target.
#    8.  strand - defined as + (forward) or - (reverse) for query strand. In mouse, a second '+' or '-' indecates genomic strand.
#    9.  qName - Query sequence name.
#    10. qSize - Query sequence size.
#    11. qStart - Alignment start position in query.
#    12. qEnd - Alignment end position in query.
#    13. tName - Target sequence name.
#    14. tSize - Target sequence size.
#    15. tStart - Alignment start position in query.
#    16. tEnd - Alignment end position in query.
#    17. blockCount - Number of blocks in the alignment.
#    18. blockSizes - Comma-separated list of sizes of each block.
#    19. qStarts - Comma-separated list of start position of each block in query.
#    20. tStarts - Comma-separated list of start position of each block in target.

HITS = {}
if args.input_file:
	INPUT = open(args.input_file, "r")
else:
	INPUT = sys.stdin

if args.verbose:
	ssw("[{0}]\tReading input file\n".format(ts()))

k=0
if args.remove_header:
	headerCount = 0
else:
	headerCount = 5
for line in INPUT:
	if headerCount >= 5:
		int(line.rstrip().split()[0])
		k+=1
		linelist = line.rstrip().split("\t")
		ratio = round(float(float(int(linelist[0]) + int(linelist[2])) / int(linelist[10])), 2)
			# store strings at the end, so they don't alter any future numerical sorting of the tuples
			# store the whole line at the end, so it's easy to fetch when printing
		lineTuple = (ratio, linelist[1], linelist[5], linelist[7], linelist[13], line.rstrip())
		try:
			HITS[linelist[9]].append(lineTuple)
		except KeyError:
			HITS[linelist[9]] = [lineTuple]
	else:
		headerCount += 1

INPUT.close()
if args.verbose:
	ssw("[{0}]\tProcessed {1} lines\n".format(ts(), k))

# select best hit
printed = 0
same_match = 0
same_mismatch = 0
same_db_ins = 0
same_query_ins = 0
dupl = 0
ins = 0

if args.output_file:
	OUTPUT = open(args.output_file, "w")
else:
	OUTPUT = sys.stdout

if args.ins_file:
	INS = open(args.ins_file, "w")
if args.dupl_file:
	DUPL = open(args.dupl_file, "w")

if args.verbose:
	ssw("[{0}]\tSelecting a unique hit for each query sequence. Choosing the best hit\n".format(ts()))

for name in HITS:
		# sort reversely (from great to small values) the tuples list
	x = sorted(HITS[name], reverse=True)
		# pick the first tuple which represents the highest number of matching bases
	bestMatch = x[0][0]
	templist = []
	for entry in x:
		if entry[0] == bestMatch:
			templist.append(entry)
	if len(templist) == 1:
		printed += 1
		OUTPUT.write(templist[0][-1] + "\n")
	# if there are more than one result having the same score, go through the mismatches
	else:
		same_match += 1
		x = sorted(templist, reverse=True)
		bestMatch = x[-1][1]
		templist = []
		for entry in x:
			if entry[1] == bestMatch:
				templist.append(entry)
		if len(templist) == 1:
			printed += 1
			OUTPUT.write(templist[-1][-1] + "\n")
		# if still there are more than one, process insertions in the database sequence
		else:
			same_mismatch += 1
			x = sorted(templist, reverse=True)
			bestMatch = x[-1][2]
			templist = []
			for entry in x:
				if entry[2] == bestMatch:
					templist.append(entry)
			if len(templist) == 1:
				printed += 1
				ins += 1
				if not args.no_insertions:
					OUTPUT.write(templist[-1][-1] + "\n")
				elif args.no_insertions:
					if args.ins_file:
						print >>INS, templist[-1][-1]
			# if still > 1, process insertions in the query
			else:
				same_db_ins += 1
				x = sorted(templist, reverse=True)
				bestMatch = x[-1][3]
				templist = []
				for entry in x:
					if entry[3] == bestMatch:
						templist.append(entry)
				ins += 1
				if not args.no_insertions:
					printed += 1
					OUTPUT.write(templist[0][-1] + "\n")
				elif args.no_insertions:
					if args.ins_file:
						print >>INS, templist[0][-1]
				if len(templist) > 1:
					same_query_ins += 1
					del templist[0]
					for element in templist:
						ins += 1
						dupl += 1
						if args.dupl_file:
							print >>DUPL, element[-1]
						if args.ins_file:
							print >>INS, element[-1]


if args.output_file:
	OUTPUT.close()
else:
	pass

if args.dupl_file:
	DUPL.close()
if args.ins_file:
	INS.close()

if args.verbose:
	ssw("[{0}]\tA total of {1} unique entries were printed to output. Of these:\n\t\t\t\t- {2} times a conflict of matching bases was found\n\t\t\t\t- {3} times a conflict of mismatch number was found\n\t\t\t\t- {4} times a conflict of database insertions was found\n\t\t\t\t- {5} times a conflict of query insertions was found\n\t\t\t\t- A total of {6} duplicate lines were found. If --dupl-file is used, these lines are there\n".format(ts(), printed, same_match, same_mismatch, same_db_ins, same_query_ins, dupl))
if args.ins_file:
	if args.verbose:
		ssw("[{0}]\tYour insertion lines file contains {1} lines\n".format(ts(), ins))
