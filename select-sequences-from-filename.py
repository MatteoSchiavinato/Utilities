#!/usr/bin/env python2.7

import sys
import argparse as ap
from time import asctime as at
from Bio import SeqIO

# help
if len(sys.argv) < 2:
	sys.argv.append("-h")
if sys.argv[1] in ["-h", "-help", "--help", "usage", "getopt"]:
	sys.exit('''

select-sequences-from-filename.py
v2.0\tlast update: 23 May 2017

-----------------------------------------------------------------------------------------------------------

-in	--input-file			FASTA input file 					[mandatory]
-out	--output-file			FASTA output file 					[mandatory]
-n	--names				Select only the names specified in this file 		[off]
	--invert			Invert selection (leave -n out)				[off]
-r	--remove-nonseq-chars		Remove non sequence characters (like *) from sequences	[off]

''')

# parser
p = ap.ArgumentParser()
p.add_argument("-in", "--input-file")
p.add_argument("-out", "--output-file")
p.add_argument("--invert", action="store_true")
p.add_argument("-n", "--names", metavar="NAMES")
p.add_argument("-r", "--remove-nonseq-chars", action="store_true")
args = p.parse_args()


# read input file
sys.stderr.write("[{0}]\tImporting list of names\n".format(at()))
FOFN = open(args.names, "r")
NAMES = []
for line in FOFN:
	NAMES.append(line.rstrip().split()[0].strip(">"))
#
FOFN.close()
NAMES = list(set(NAMES))
sys.stderr.write("[{0}]\tStored {1} names\n".format(at(), len(NAMES)))


# parsing input fasta and writing to output
k=0
sys.stderr.write("[{0}]\tParsing input FASTA and writing to output lines matching name list\n".format(at()))

INPUT = open(args.input_file, "r")
OUTPUT = open(args.output_file, "w")
WrittenRecords = []

for record in SeqIO.parse(INPUT, "fasta"):
	name = str(record.id)
	if args.invert:
		if (name not in NAMES):
			if args.remove_nonseq_chars:
				sequence = str(record.seq).strip("*-_[].()~&%$\"?+:,;#\'")
			else:
				sequence = str(record.seq)
			OUTPUT.write(">" + str(record.id) + "\n" + sequence + "\n")
			WrittenRecords.append(str(record.id))
			k+=1
	else:
		if (name in NAMES):
                        if args.remove_nonseq_chars:
                                sequence = str(record.seq).strip("*-_[].()~&%$\"?+:,;#\'")
                        else:
                                sequence = str(record.seq)
                        OUTPUT.write(">" + str(record.id) + "\n" + sequence + "\n")
                        WrittenRecords.append(str(record.id))
                        k+=1

INPUT.close()
OUTPUT.close()
sys.stderr.write("[{0}]\tCreated a subset containing {1} sequences\n".format(at(), k))
